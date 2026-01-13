"""Example non-conformal tie coupling in FEniCSx using Nitsche method.

This is a template that loads per-part meshes and demonstrates how to:
  - read facet tags from XDMF
  - identify interface facets on each part
  - assemble a Nitsche-style coupling term

NOTE: This is a starting point. You should tune penalty parameters and
choose a coupling formulation appropriate to your problem.
"""

from __future__ import annotations

from pathlib import Path
from mpi4py import MPI
from dolfinx import fem, io
import ufl


PARTS_DIR = Path("fem/results/parts_cshape")


def _load_mesh(xdmf_path: Path):
    with io.XDMFFile(MPI.COMM_WORLD, str(xdmf_path), "r") as xdmf:
        domain = xdmf.read_mesh(name="Grid")
        try:
            facet_tags = xdmf.read_meshtags(domain, name="Grid")
        except Exception:
            facet_tags = None
    return domain, facet_tags


def _nitsche_coupling(u, v, u2, v2, n, h, gamma=10.0):
    """Basic symmetric Nitsche coupling for displacement continuity."""
    # Jump/average terms
    jump_u = u - u2
    jump_v = v - v2
    return gamma / h * ufl.inner(jump_u, jump_v) * ufl.ds


def main():
    # Load x_beam and riser_left as an example pair
    x_beam, x_facets = _load_mesh(PARTS_DIR / "x_beam.xdmf")
    riser, r_facets = _load_mesh(PARTS_DIR / "riser_left.xdmf")

    if x_facets is None or r_facets is None:
        raise RuntimeError("Facet tags missing. Run fenicsx_tie_parts.py first.")

    # Function spaces
    Vx = fem.functionspace(x_beam, ("Lagrange", 1, (x_beam.geometry.dim,)))
    Vr = fem.functionspace(riser, ("Lagrange", 1, (riser.geometry.dim,)))
    ux = ufl.TrialFunction(Vx)
    vx = ufl.TestFunction(Vx)
    ur = ufl.TrialFunction(Vr)
    vr = ufl.TestFunction(Vr)

    # Extract facet tag IDs (replace with your actual IDs if different)
    # These should match tags in the .msh
    x_iface_id = 2  # x_beam_end_left
    r_iface_id = 2  # riser_left_inner

    # Define measures on each mesh
    ds_x = ufl.Measure("ds", domain=x_beam, subdomain_data=x_facets)
    ds_r = ufl.Measure("ds", domain=riser, subdomain_data=r_facets)

    # Normal vectors
    nx = ufl.FacetNormal(x_beam)
    nr = ufl.FacetNormal(riser)

    # Mesh sizes for penalty scaling
    hx = ufl.CellDiameter(x_beam)
    hr = ufl.CellDiameter(riser)

    # Coupling term (toy example; requires mapping between non-matching meshes)
    # In practice, you will need mortar projection or custom coupling operator.
    # This placeholder shows where the term would go.
    a_couple = (
        _nitsche_coupling(ux, vx, ur, vr, nx, hx, gamma=20.0) * ds_x(x_iface_id)
        + _nitsche_coupling(ur, vr, ux, vx, nr, hr, gamma=20.0) * ds_r(r_iface_id)
    )

    print("Coupling form assembled (placeholder).")
    print("Replace with mortar/Nitsche implementation for non-matching meshes.")


if __name__ == "__main__":
    main()
