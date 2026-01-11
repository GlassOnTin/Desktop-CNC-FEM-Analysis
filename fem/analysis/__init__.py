"""Analysis solvers for static and modal FEM.

Two solver backends are available:
- fenicsx: Modern FEniCSx solver (recommended)
- sfepy: Legacy SfePy solver
"""

# FEniCSx solvers (default)
try:
    from .fenicsx_static import solve_static, save_results_xdmf
    from .fenicsx_modal import solve_modal, save_mode_shapes, chatter_analysis
    FENICSX_AVAILABLE = True
except ImportError:
    FENICSX_AVAILABLE = False

# SfePy solvers (legacy)
try:
    from .static_solver import setup_gantry_problem, solve_with_cutting_load
    from .modal_solver import solve_modal as sfepy_solve_modal
    SFEPY_AVAILABLE = True
except ImportError:
    SFEPY_AVAILABLE = False

__all__ = [
    'FENICSX_AVAILABLE',
    'SFEPY_AVAILABLE',
]

if FENICSX_AVAILABLE:
    __all__.extend(['solve_static', 'save_results_xdmf', 'solve_modal', 'save_mode_shapes', 'chatter_analysis'])

if SFEPY_AVAILABLE:
    __all__.extend(['setup_gantry_problem', 'solve_with_cutting_load', 'sfepy_solve_modal'])
