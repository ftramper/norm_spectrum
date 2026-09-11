# norm_spectrum

Lightweight, importable normalization routines for astronomical spectra.

The first routine, `pseudo_norm`, performs a broad running-median normalization. For spectra with a well-defined continuum (for example many O-star spectra), this behaves as a conventional continuum normalization; for spectra dominated by dense absorption structure, it provides a pseudo-continuum normalization. The method is intended for fast analysis tasks such as plotting and radial-velocity measurements rather than detailed atmosphere modelling.

## Installation

From a clone of this repository:

```bash
pip install .
```

## Usage

```python
from norm_spectrum import pseudo_norm

norm_flux = pseudo_norm(wavelength, flux, window_width=180.0)
```

With uncertainties:

```python
norm_flux, norm_err = pseudo_norm(
    wavelength, flux, error, window_width=180.0
)
```

`window_width` is expressed in the same units as the wavelength array. A fixed window in samples can instead be supplied with `window_width=None, window_pixels=101`.

The package is intentionally being structured so additional normalization methods can be added later under the same namespace.
