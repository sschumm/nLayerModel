# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi


def K_ph_amplitude(T_ph, r, I, n=1, o=1, m=3):
    b = n * o/2
    k = np.sin(b * pi) / b
    return I * (m/pi) * (T_ph * k)/r 


def K_n_amplitude(T, r, I, n=1, o=1):
    b = n * o/2
    k = np.sin(b * pi) / b
    return I * (2 * T * k)/(pi * r)


def A(z, i, d):
    return (z * i) / (d * pi)