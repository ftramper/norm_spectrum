# norm_spectrum

Lightweight, importable normalization routines for astronomical spectra.

The first routine, `pseudo_norm`, performs a broad running-median normalization. For spectra with a well-defined continuum (for example many O-star spectra), this behaves as a conventional continuum normalization; for spectra dominated by dense absorption structure, it provides a pseudo-continuum normalization. The method is intended for fast analysis tasks such as plotting and radial-velocity measurements rather than detailed atmosphere modelling.

Additional normalization methods can be added later under the same package namespace (for example the X-shooter NIR normalization routine).

## Method

`pseudo_norm` estimates a broad continuum or pseudo-continuum with a running median and divides the original, unsmoothed spectrum by that estimate.

Invalid flux samples are interpolated only while constructing the continuum. They remain invalid in the returned normalized spectrum. The wavelength array does not need to be sorted beforehand; the routine sorts internally and returns the result in the original input order.

The running-median width can be specified either in wavelength units or directly in pixels.

## Installation

From a clone of this repository:

```bash
pip install .
```

A pip installation is not required. Alternatively, copy the `norm_spectrum` package directory (the directory under `src/`) to a location included in your `PYTHONPATH`, or place it in the directory from which you run your Python code. In either case it can be imported in exactly the same way.

The package can also be installed directly from GitHub:

```bash
pip install git+https://github.com/ftramper/norm_spectrum.git
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

The uncertainty array, when supplied, is divided by the same continuum estimate as the flux.

### Choosing the continuum window

`window_width` is expressed in the same units as the wavelength array. The default value is `180.0`, matching the original prescription when wavelengths are supplied in Angstrom.

For example, the same physical 180 Angstrom window for wavelengths stored in micron is:

```python
norm_flux = pseudo_norm(wavelength_um, flux, window_width=0.018)
```

A fixed width in samples can instead be used with `window_pixels`:

```python
norm_flux = pseudo_norm(
    wavelength, flux,
    window_width=None,
    window_pixels=101,
)
```

The appropriate window depends on the spectral resolution, wavelength range, and scientific purpose. Broad windows preserve narrow and moderately broad spectral features while following slower continuum structure.

## Public API

Currently available:

```python
from norm_spectrum import pseudo_norm
```

A future release is intended to add the X-shooter normalization routine under the same namespace, e.g.:

```python
from norm_spectrum import RJ_norm
```
