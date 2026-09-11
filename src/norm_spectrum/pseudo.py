"""Quick continuum and pseudo-continuum normalization for 1-D spectra."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import median_filter

DEFAULT_CONTINUUM_WINDOW = 180.0


def _odd_size(value: int) -> int:
    """Return an odd filter size, with a minimum of five samples."""
    value = max(5, int(value))
    return value if value % 2 else value + 1


def _validate_inputs(wavelength, flux, error=None):
    """Convert inputs to 1-D float arrays and check compatible shapes."""
    wave = np.asarray(wavelength, dtype=float)
    flx = np.asarray(flux, dtype=float)

    if wave.ndim != 1 or flx.ndim != 1:
        raise ValueError("wavelength and flux must be one-dimensional arrays")
    if wave.shape != flx.shape:
        raise ValueError("wavelength and flux must have identical shapes")

    err = None
    if error is not None:
        err = np.asarray(error, dtype=float)
        if err.ndim != 1:
            raise ValueError("error must be a one-dimensional array")
        if err.shape != wave.shape:
            raise ValueError("error must have the same shape as wavelength and flux")

    return wave, flx, err


def _resolve_window_pixels(
    wavelength_sorted,
    flux_sorted,
    *,
    window_width,
    window_pixels,
):
    """Convert the requested continuum window to an odd pixel width."""
    if window_width is not None and window_pixels is not None:
        raise ValueError("Specify either window_width or window_pixels, not both")

    if window_pixels is not None:
        if (
            not np.isfinite(window_pixels)
            or float(window_pixels) <= 0
            or not float(window_pixels).is_integer()
        ):
            raise ValueError("window_pixels must be a positive integer")
        return _odd_size(int(window_pixels))

    if window_width is None:
        window_width = DEFAULT_CONTINUUM_WINDOW

    if not np.isfinite(window_width) or float(window_width) <= 0:
        raise ValueError("window_width must be finite and > 0")

    good = np.isfinite(wavelength_sorted) & np.isfinite(flux_sorted)
    if np.count_nonzero(good) < 2:
        raise ValueError("At least two finite wavelength/flux samples are required")

    good_wave = wavelength_sorted[good]
    steps = np.diff(good_wave)
    steps = steps[np.isfinite(steps) & (steps > 0)]
    if steps.size == 0:
        raise ValueError("Cannot determine a positive wavelength sampling step")

    step = float(np.nanmedian(steps))
    return _odd_size(round(float(window_width) / step))


def estimate_pseudocontinuum(
    wavelength,
    flux,
    *,
    window_width=DEFAULT_CONTINUUM_WINDOW,
    window_pixels=None,
):
    """Estimate a broad running-median continuum or pseudo-continuum.

    ``window_width`` is expressed in the same units as ``wavelength``. Invalid
    flux samples are interpolated only while constructing the continuum.
    """
    wave, flx, _ = _validate_inputs(wavelength, flux)

    finite_wave = np.isfinite(wave)
    good = finite_wave & np.isfinite(flx)
    if np.count_nonzero(good) < 2:
        raise ValueError("At least two finite wavelength/flux samples are required")

    finite_indices = np.flatnonzero(finite_wave)
    order = np.argsort(wave[finite_indices])
    sorted_indices = finite_indices[order]
    wave_sorted = wave[sorted_indices]
    flux_sorted = flx[sorted_indices]

    good_sorted = np.isfinite(flux_sorted)
    good_wave = wave_sorted[good_sorted]
    good_flux = flux_sorted[good_sorted]

    size = _resolve_window_pixels(
        wave_sorted,
        flux_sorted,
        window_width=window_width,
        window_pixels=window_pixels,
    )

    flux_interp = np.interp(wave_sorted, good_wave, good_flux)
    continuum_sorted = median_filter(flux_interp, size=size, mode="nearest")

    finite_cont = np.isfinite(continuum_sorted) & (continuum_sorted != 0.0)
    if np.any(finite_cont):
        replacement = float(np.nanmedian(continuum_sorted[finite_cont]))
    else:
        replacement = 1.0
    continuum_sorted[~finite_cont] = replacement

    continuum = np.full(wave.shape, np.nan, dtype=float)
    continuum[sorted_indices] = continuum_sorted
    return continuum


def pseudo_norm(
    wavelength,
    flux,
    error=None,
    *,
    window_width=DEFAULT_CONTINUUM_WINDOW,
    window_pixels=None,
):
    """Normalize a 1-D spectrum to a broad running-median continuum.

    For spectra with a true continuum (for example many O-star spectra), the
    result is a conventional continuum normalization. For strongly line-blanketed
    or structured spectra, the result is a pseudo-continuum normalization.

    Parameters
    ----------
    wavelength : array-like
        One-dimensional wavelength array in any unit.
    flux : array-like
        One-dimensional flux array.
    error : array-like, optional
        Formal uncertainty array. If supplied, it is divided by the absolute
        continuum and returned as a second array.
    window_width : float, optional
        Running-median width in the same units as ``wavelength``. Default 180.0.
    window_pixels : int, optional
        Running-median width in pixels. Set ``window_width=None`` when using it.

    Returns
    -------
    normalized_flux : numpy.ndarray
        Normalized flux, retaining NaNs/non-finite input flux samples.
    normalized_error : numpy.ndarray, optional
        Returned only if ``error`` is supplied.
    """
    wave, flx, err = _validate_inputs(wavelength, flux, error)
    continuum = estimate_pseudocontinuum(
        wave,
        flx,
        window_width=window_width,
        window_pixels=window_pixels,
    )

    normalized_flux = np.full(flx.shape, np.nan, dtype=float)
    usable = np.isfinite(flx) & np.isfinite(continuum) & (continuum != 0.0)
    normalized_flux[usable] = flx[usable] / continuum[usable]

    if err is None:
        return normalized_flux

    normalized_error = np.full(err.shape, np.nan, dtype=float)
    usable_err = usable & np.isfinite(err)
    normalized_error[usable_err] = err[usable_err] / np.abs(continuum[usable_err])
    return normalized_flux, normalized_error


__all__ = [
    "DEFAULT_CONTINUUM_WINDOW",
    "estimate_pseudocontinuum",
    "pseudo_norm",
]
