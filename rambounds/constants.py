"""Physical constants and unit conventions used throughout the package.

Units: frequencies in GHz, thicknesses in mm, magnetic masses in kg/m^2,
saturation magnetization 4*pi*Ms in kG, densities in g/cm^3.
Attenuations are given in dB by the user and converted to nepers internally.
"""
import numpy as np

C = 299.792458                      # speed of light, mm.GHz
DB_PER_NEPER = 20.0 / np.log(10.0)  # 8.686 dB per neper (amplitude)
K_DB = 2.0 * np.pi**2 * DB_PER_NEPER  # 171.45 dB, so that c/K_DB = 1.75 mm.GHz/dB
GAMMA_BAR = 2.8                     # reduced gyromagnetic ratio: 2.8 MHz/Oe = 2.8 GHz per kG


def db_to_np(x_db):
    """Attenuation in dB -> nepers."""
    return np.asarray(x_db, dtype=float) / DB_PER_NEPER


def np_to_db(x_np):
    """Attenuation in nepers -> dB."""
    return np.asarray(x_np, dtype=float) * DB_PER_NEPER
