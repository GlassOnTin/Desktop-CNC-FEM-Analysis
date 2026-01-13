"""Load TTC450 parametric part meshes and prepare interface tags for FEniCSx coupling.

This script verifies interface facet tags in the per-part XDMF files and
prints the tag IDs needed for non-conformal coupling (e.g., Nitsche or MPC).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile


PARTS = {
    "x_beam": {
        "file": "x_beam",
        "iface_names": ["x_beam_end_left", "x_beam_end_right"],
    },
    "riser_left": {
        "file": "riser_left",
        "iface_names": ["riser_left_inner"],
    },
    "riser_right": {
        "file": "riser_right",
        "iface_names": ["riser_right_inner"],
    },
}


def _read_mesh_and_tags(msh_path: Path):
    """Load mesh and tags directly from Gmsh .msh using dolfinx."""
    domain, cell_tags, facet_tags = gmshio.read_from_msh(
        str(msh_path), MPI.COMM_WORLD, gdim=3
    )
    return domain, cell_tags, facet_tags


def _write_xdmf(msh_path: Path, xdmf_path: Path) -> None:
    """Write mesh + tags to XDMF for FEniCSx workflows."""
    domain, cell_tags, facet_tags = _read_mesh_and_tags(msh_path)
    with XDMFFile(MPI.COMM_WORLD, str(xdmf_path), "w") as xdmf:
        xdmf.write_mesh(domain)
        if cell_tags is not None:
            xdmf.write_meshtags(cell_tags)
        if facet_tags is not None:
            xdmf.write_meshtags(facet_tags)


def _resolve_tag_ids(field_data: Dict[str, Tuple[int, int]], names: list[str]) -> Dict[str, int]:
    """Resolve interface names to physical IDs for facet tags."""
    ids = {}
    for name in names:
        if name not in field_data:
            raise ValueError(f"Missing physical name '{name}' in field_data.")
        tag_id, dim = field_data[name]
        if dim != 2:
            raise ValueError(f"Physical name '{name}' has dim={dim}, expected 2.")
        ids[name] = tag_id
    return ids


def main(parts_dir: Path, write_xdmf: bool = True) -> None:
    parts_dir = Path(parts_dir)
    for part_name, cfg in PARTS.items():
        msh_path = parts_dir / f"{cfg['file']}.msh"
        xdmf_path = parts_dir / f"{cfg['file']}.xdmf"

        if not msh_path.exists():
            raise FileNotFoundError(f"Missing mesh file for {part_name}")

        domain, cell_tags, facet_tags = _read_mesh_and_tags(msh_path)
        if facet_tags is None:
            raise RuntimeError(f"No facet tags found in {msh_path.name}")

        # Read physical names via gmsh API (msh stores name->id)
        try:
            import gmsh
            gmsh.initialize()
            gmsh.open(str(msh_path))
            field_data = {gmsh.model.getPhysicalName(dim, tag): (tag, dim)
                          for dim, tag in gmsh.model.getPhysicalGroups()}
            gmsh.finalize()
        except Exception:
            field_data = {}

        tag_ids = _resolve_tag_ids(field_data, cfg["iface_names"])

        print(f"{part_name}: mesh={msh_path.name}")
        print(f"  interface tag IDs: {tag_ids}")

        for name, tag_id in tag_ids.items():
            count = (facet_tags.values == tag_id).sum()
            print(f"  {name}: {count} facets")

        if write_xdmf:
            _write_xdmf(msh_path, xdmf_path)
            print(f"  wrote: {xdmf_path.name}")

    print("\nNext step: use these tag IDs for non-conformal coupling in FEniCSx.")
    print("Recommended: Nitsche or MPC-based tying with coincident interfaces.")


if __name__ == "__main__":
    main(Path("fem/results/parts_cshape"), write_xdmf=True)
