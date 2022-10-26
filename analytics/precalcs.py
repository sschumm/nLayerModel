# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi


# zone factor
def kd(m, q):
    nom = np.sin(pi/(2*m))
    den = q * np.sin(pi/(2 * m*q))
    return nom / den

# pitching factor
def kp(w, tp, n=1):
    return np.sin(n * (pi/2) * (w/tp))

# breadth factor
def kb(n, sigma):
    return np.sinc(n * sigma/2)

# current loading
def K(m, I, d, N):
    num = 2 * m * N * I
    den = pi * d
    return num / den

# pole pitch
def taup(d_si, p):
    num = d_si * pi
    den = 2 * p
    return num / den

def q(Q, m, p):
    num = Q
    den = 2 * p * m    
    return num / den


def Ns(p, q, Nc, a):
    num = p * q * Nc
    den = a
    return num / den

def yoke_height(p, r_si, B, B_max):
    num = r_si * B
    den = p * B_max
    return num / den
