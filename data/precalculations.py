# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi


def kd(m, q):
    nom = np.sin(pi/(2*m))
    den = q * np.sin(pi/(2 * m*q))
    return nom / den


def kw(n, o):
    num = np.sin(n * o/2)
    den = n * o/2
    return num / den


def q(Q, m, p):
    num = Q
    den = 2 * p * m    
    return num / den


def Ns(p, q, Nc, a):
    num = p * q * Nc
    den = a
    return num / den


def K(m, I, d, N):
    num = 2 * m * N * I
    den = pi * d
    return num / den


def _K(m, I, d, p, q, Nc, a):
    num = 2 * m * Ns(p, q, Nc, a) * I
    den = pi * d
    return num / den


def Kamp(kw, K):
    return np.sqrt(2) * kw * K





def K_ph_amplitude(T_ph, r, I, n=1, o=1, m=3):
    b = n * pi * o/2
    k = np.sin(b) / b
    return I * (m/pi) * (T_ph * k)/r 


def K_n_amplitude(T, r, I, n=1, o=1):
    b = n * pi * o/2
    k = np.sin(b) / b
    return I * (2 * T * k)/(pi * r)


def A(z, i, d):
    return (z * i) / (d * pi)