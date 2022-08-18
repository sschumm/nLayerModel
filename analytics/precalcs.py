# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi


# zone factor
def kd(m, q):
    nom = np.sin(pi/(2*m))
    den = q * np.sin(pi/(2 * m*q))
    return nom / den

# breadth factor
def kb(n, sigma):
    return np.sinc(n * sigma/2)
    # num = np.sin(n * sigma/2)
    # den = n * sigma/2
    # return num / den

# current loading
def K(m, I, d, N):
    num = 2 * m * N * I
    den = pi * d
    return num / den


def q(Q, m, p):
    num = Q
    den = 2 * p * m    
    return num / den


def Ns(p, q, Nc, a):
    num = p * q * Nc
    den = a
    return num / den


