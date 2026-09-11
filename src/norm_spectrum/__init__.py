"""Tools for quick normalization of astronomical spectra."""

from .pseudo import DEFAULT_CONTINUUM_WINDOW, estimate_pseudocontinuum, pseudo_norm

__all__ = ["DEFAULT_CONTINUUM_WINDOW", "estimate_pseudocontinuum", "pseudo_norm"]
__version__ = "0.1.0"
