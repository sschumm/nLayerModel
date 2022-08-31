# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi

# get current loading from Essons Number and Air Gap Flux Density Amplitude
def K_from(B_airgap_amplitude, C_e, kw=1):
    nom = np.sqrt(2) * C_e
    den = pi**2 * kw * B_airgap_amplitude
    return nom / den


# get Air Gap Flux Density Amplitude from Essons Number and current loading
def B_from(K, C_e, kw=1):
    nom = np.sqrt(2) * C_e
    den = pi**2 * kw * K
    return nom / den


# get main flux per pole from Air Gap Flux Density Amplitude and dimensions
def Phi_from(B_airgap_amplitude, taup, l_e):
    nom = 2 * taup * l_e * B_airgap_amplitude
    den = pi    
    return nom / den