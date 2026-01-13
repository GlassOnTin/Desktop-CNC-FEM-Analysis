"""CLI wrapper to generate combined mesh and run MPC tie solve.

Usage:
  python fem/analysis/run_ttc450_mpc.py --gap 0.2 --mesh-min 8 --mesh-max 25
"""

from __future__ import annotations

import argparse
from pathlib import Path

from fem.geometry.mesh_generator import generate_ttc450_combined_mesh


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate TTC450 combined mesh and run MPC tie solve.")
    parser.add_argument("--gap", type=float, default=0.2, help="Interface gap (mm)")
    parser.add_argument("--mesh-min", type=float, default=8.0, help="Min mesh size (mm)")
    parser.add_argument("--mesh-max", type=float, default=25.0, help="Max mesh size (mm)")
    parser.add_argument("--profile", type=str, default="cshape", choices=["cshape", "hollow", "solid", "cslot"])
    parser.add_argument("--include-base", action="store_true", help="Include base frame")
    parser.add_argument("--mesh-only", action="store_true", help="Generate mesh only; skip MPC solve")
    parser.add_argument("--self-check", action="store_true", help="Validate facet tags after meshing")
    parser.add_argument("--write-xdmf", action="store_true", help="Write XDMF for visualization")
    parser.add_argument("--write-tags-json", action="store_true", help="Write tag IDs + counts to JSON")
    args = parser.parse_args()

    output_msh = Path("fem/results/ttc450_combined.msh")
    generate_ttc450_combined_mesh(
        output_msh,
        mesh_size_min=args.mesh_min,
        mesh_size_max=args.mesh_max,
        profile=args.profile,
        interface_gap=args.gap,
        include_base=args.include_base,
    )

    if args.self_check or args.write_xdmf or args.write_tags_json:
        from fem.analysis.fenicsx_tie_parts import _read_mesh_and_tags
        from dolfinx.io import XDMFFile
        from mpi4py import MPI
        import json
        import gmsh

        domain, cell_tags, facet_tags = _read_mesh_and_tags(output_msh)
        if args.self_check:
            if facet_tags is None:
                raise RuntimeError("No facet tags found in combined mesh.")
            unique = sorted(set(facet_tags.values.tolist()))
            print("Facet tag IDs:", unique)
            for tag_id in unique:
                count = (facet_tags.values == tag_id).sum()
                print(f"  tag {tag_id}: {count} facets")

        if args.write_xdmf:
            xdmf_path = output_msh.with_suffix(".xdmf")
            with XDMFFile(MPI.COMM_WORLD, str(xdmf_path), "w") as xdmf:
                xdmf.write_mesh(domain)
                if cell_tags is not None:
                    xdmf.write_meshtags(cell_tags)
                if facet_tags is not None:
                    xdmf.write_meshtags(facet_tags)
            print(f"Wrote {xdmf_path}")

        if args.write_tags_json:
            if facet_tags is None:
                raise RuntimeError("No facet tags found in combined mesh.")
            gmsh.initialize()
            gmsh.open(str(output_msh))
            field_data = {gmsh.model.getPhysicalName(dim, tag): (tag, dim)
                          for dim, tag in gmsh.model.getPhysicalGroups()}
            gmsh.finalize()

            tag_counts = {}
            for name, (tag_id, dim) in field_data.items():
                if dim != 2:
                    continue
                count = int((facet_tags.values == tag_id).sum())
                tag_counts[name] = {"id": int(tag_id), "dim": int(dim), "count": count}

            json_path = output_msh.with_name(f"{output_msh.stem}_tags.json")
            with open(json_path, "w") as f:
                json.dump(tag_counts, f, indent=2)
            print(f"Wrote {json_path}")

    if args.mesh_only:
        print(f"Mesh written to {output_msh}")
        return

    # Run MPC solve
    from fem.analysis import fenicsx_mpc_tie_combined

    fenicsx_mpc_tie_combined.INTERFACE_GAP = args.gap
    fenicsx_mpc_tie_combined.MESH_PATH = output_msh
    fenicsx_mpc_tie_combined.main()


if __name__ == "__main__":
    main()
