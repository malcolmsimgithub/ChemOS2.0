from typing import Optional, Tuple
import numpy as np
from scipy.constants import c, nano


def calculate_gain_factor(
        emission_spectrum: np.ndarray,
        spectral_interval: Optional[Tuple[float, float]] = None,
        gain_interval: Optional[Tuple[float, float]] = None,
        maximum: bool = True
) -> Tuple[float, float]:
    """
    Calculates the maximum gain factor for a given spectrum.

    Args:
        emission_spectrum: 2D Numpy array (2, N) containing the wavelength (in nm) and intensity values.
        spectral_interval: Wavelength interval to take into account for calculating the gain spectrum (min, max).
        gain_interval: Wavelength interval to take into account for identifying the maximum gain factor (should not
                       overlap with the absorption spectrum of the sample).
        maximum: If True, the maximum gain factor is returned. If False, the gain factor at the maximum emission
                 wavelength is returned.

    Returns:
        float: Maximum gain factor (in m^2 s).
        float: Wavelength at which the maximum gain factor occurs (in nm).
    """
    # mask the emission spectrum to the spectral interval
    if spectral_interval is not None:
        emission_spectrum = emission_spectrum[:, np.logical_and(emission_spectrum[0] >= spectral_interval[0], emission_spectrum[0] <= spectral_interval[1])]

    # convert to SI units
    emission_spectrum[0] *= nano

    # calculate the gain spectrum as lambda ** 4 * emission_spectrum / (8 * pi * c * int(emission_spectrum))
    gain_spectrum = emission_spectrum[0] ** 4 * emission_spectrum[1] / (8 * np.pi * c * np.trapz(emission_spectrum[1], emission_spectrum[0]))

    # mask the gain spectrum to the gain interval
    if gain_interval is not None:
        gain_spectrum = gain_spectrum[np.logical_and(gain_spectrum[0] >= gain_interval[0], gain_spectrum[0] <= gain_interval[1])]

    # calculate the maximum gain factor
    if maximum:
        max_idx = np.argmax(gain_spectrum[1])
    else:
        max_idx = np.argmax(emission_spectrum[1])

    return gain_spectrum[1][max_idx], gain_spectrum[0][max_idx] / nano












