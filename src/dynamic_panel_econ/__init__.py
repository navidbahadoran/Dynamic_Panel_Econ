"""Low-rank dynamic-panel Monte Carlo replication package."""

from .core import Coefficients, Design, adjoint, fitted_values
from .dgp import DGPParameters, PanelData, generate_panel

__all__ = [
    "Coefficients",
    "DGPParameters",
    "Design",
    "PanelData",
    "adjoint",
    "fitted_values",
    "generate_panel",
]

__version__ = "0.1.0"
