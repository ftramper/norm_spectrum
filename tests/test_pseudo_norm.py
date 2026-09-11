import numpy as np

from norm_spectrum import pseudo_norm


def test_flat_spectrum_normalizes_to_one():
    wave = np.arange(100, dtype=float)
    flux = np.full(100, 7.5)
    norm = pseudo_norm(wave, flux, window_width=11.0)
    assert np.allclose(norm, 1.0)


def test_error_is_scaled_by_same_continuum():
    wave = np.arange(100, dtype=float)
    flux = np.full(100, 10.0)
    err = np.full(100, 2.0)
    norm, norm_err = pseudo_norm(wave, flux, err, window_width=11.0)
    assert np.allclose(norm, 1.0)
    assert np.allclose(norm_err, 0.2)


def test_invalid_flux_is_not_filled_in_output():
    wave = np.arange(100, dtype=float)
    flux = np.full(100, 5.0)
    flux[50] = np.nan
    norm = pseudo_norm(wave, flux, window_width=11.0)
    assert np.isnan(norm[50])
    assert np.allclose(norm[np.isfinite(norm)], 1.0)
