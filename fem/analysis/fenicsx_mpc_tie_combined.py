"""MPC-based tie coupling for the combined TTC450 mesh in FEniCSx.

This expects a single .msh containing multiple volumes and facet tags:
  - x_beam_end_left / x_beam_end_right
  - riser_left_inner / riser_right_inner
and volume names for each part.
"""

from __future__ import annotations

from pathlib import Path
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

import ufl
from dolfinx import fem
from dolfinx.io import gmshio

try:
    import dolfinx_mpc
    from dolfinx_mpc import MultiPointConstraint, LinearProblem
    from dolfinx_mpc.utils import create_contact_inelastic_condition
except Exception as exc:  # pragma: no cover
    raise RuntimeError(
        "dolfinx_mpc is required for this script. "
        "Install it and re-run."
    ) from exc


MESH_PATH = Path("fem/results/ttc450_combined.msh")
INTERFACE_GAP = 0.2  # mm, must match mesh generation


def main():
    # Load combined mesh + tags
    domain, cell_tags, facet_tags = gmshio.read_from_msh(
        str(MESH_PATH), MPI.COMM_WORLD, gdim=3
    )

    # Function space
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # MPC constraints
    mpc = MultiPointConstraint(V)

    # Resolve tag IDs from gmsh
    import gmsh
    gmsh.initialize()
    gmsh.open(str(MESH_PATH))
    name_to_id = {gmsh.model.getPhysicalName(dim, tag): tag
                  for dim, tag in gmsh.model.getPhysicalGroups()}
    gmsh.finalize()

    x_left = name_to_id["x_beam_end_left"]
    x_right = name_to_id["x_beam_end_right"]
    r_left = name_to_id["riser_left_inner"]
    r_right = name_to_id["riser_right_inner"]

    # Create MPC ties (allow small separation)
    # eps2 is squared distance tolerance in physical units
    eps2 = (INTERFACE_GAP * 1.5) ** 2
    create_contact_inelastic_condition(mpc, facet_tags, x_left, r_left, eps2=eps2)
    create_contact_inelastic_condition(mpc, facet_tags, x_right, r_right, eps2=eps2)
    mpc.finalize()

    # Simple static elasticity example (no loads by default)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    E = 69.0e3  # N/mm^2
    nu = 0.33
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    def eps(w):
        return ufl.sym(ufl.grad(w))

    def sigma(w):
        return lmbda * ufl.nabla_div(w) * ufl.Identity(len(w)) + 2 * mu * eps(w)

    a = ufl.inner(sigma(u), eps(v)) * ufl.dx
    L = ufl.dot(fem.Constant(domain, PETSc.ScalarType((0.0, 0.0, 0.0))), v) * ufl.dx

    # Example: fix the base (if present) by bounding box
    coords = domain.geometry.x
    zmin = coords[:, 2].min()
    tol = 1e-3

    def on_base(x):
        return x[2] < zmin + tol

    dofs = fem.locate_dofs_geometrical(V, on_base)
    bc = fem.dirichletbc(PETSc.ScalarType((0.0, 0.0, 0.0)), dofs, V)

    problem = LinearProblem(a, L, mpc=mpc, bcs=[bc])
    uh = problem.solve()

    if MPI.COMM_WORLD.rank == 0:
        print("Solved with MPC ties. Displacement norm:", np.linalg.norm(uh.x.array))


if __name__ == "__main__":
    main()
