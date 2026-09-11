import numpy as np
import pytest

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


def test_unsorted_wavelength_preserves_input_order():
    rng = np.random.default_rng(1234)
    wave = np.linspace(4000.0, 5000.0, 501)
    flux = 1.0 + 0.05 * np.sin(wave / 35.0)

    reference = pseudo_norm(wave, flux, window_width=180.0)

    order = rng.permutation(wave.size)
    shuffled = pseudo_norm(wave[order], flux[order], window_width=180.0)

    restored = np.empty_like(shuffled)
    restored[order] = shuffled
    assert np.allclose(restored, reference)


def test_window_pixels_mode():
    wave = np.linspace(4000.0, 5000.0, 501)
    flux = 1.0 + 0.05 * np.sin(wave / 35.0)

    norm_pixels = pseudo_norm(
        wave, flux, window_width=None, window_pixels=91
    )

    step = np.median(np.diff(wave))
    norm_width = pseudo_norm(wave, flux, window_width=91 * step)

    assert np.allclose(norm_pixels, norm_width)


def test_equivalent_wavelength_units_give_same_result():
    wave_angstrom = np.linspace(4000.0, 5000.0, 501)
    wave_micron = wave_angstrom / 10000.0
    flux = 1.0 + 0.05 * np.sin(wave_angstrom / 35.0)

    norm_angstrom = pseudo_norm(
        wave_angstrom, flux, window_width=180.0
    )
    norm_micron = pseudo_norm(
        wave_micron, flux, window_width=0.018
    )

    assert np.allclose(norm_angstrom, norm_micron)


def test_both_window_definitions_raise():
    wave = np.arange(20, dtype=float)
    flux = np.ones(20)

    with pytest.raises(ValueError, match="either window_width or window_pixels"):
        pseudo_norm(wave, flux, window_width=5.0, window_pixels=5)


def test_non_integer_window_pixels_raise():
    wave = np.arange(20, dtype=float)
    flux = np.ones(20)

    with pytest.raises(ValueError, match="positive integer"):
        pseudo_norm(wave, flux, window_width=None, window_pixels=5.5)


def test_mismatched_shapes_raise():
    wave = np.arange(20, dtype=float)
    flux = np.ones(19)

    with pytest.raises(ValueError, match="identical shapes"):
        pseudo_norm(wave, flux)
