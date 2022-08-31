# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi

# get current loading from Essons Number and Air Gap Flux Density Amplitude
def K_from(C_e, B_airgap_amplitude, kw=1):
    nom = np.sqrt(2) * C_e
    den = pi**2 * kw * B_airgap_amplitude
    return nom / den


# get Air Gap Flux Density Amplitude from Essons Number and current loading
def K_from(C_e, B_airgap_amplitude, kw=1):
    nom = np.sqrt(2) * C_e
    den = pi**2 * kw * B_airgap_amplitude
    return nom / den