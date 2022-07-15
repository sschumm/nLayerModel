# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:38:57 2022

@author: svens
"""
import numpy as np



def r0_to_x(r, theta):
    x = r * np.cos(theta)
    return x


def r0_to_y(r, theta):
    y = r * np.sin(theta)
    return y


def r0_to_xy(r, theta):
    x = r0_to_x(r, theta)
    y = r0_to_y(r, theta)    
    return x, y


# -----------------------

def fa(r, p):
    return r**p


def fb(r, p):
    return r**-p
    

def dfa(r, p):
    return p * r**(p-1)


def dfb(r, p):
    return -p * r**-(p+1)


# -----------------------

def solution_for_A(r, p, aj, bj, fa, fb):
    return aj * fa(r, p) + bj * fb(r, p)  


def Ar_no_k(r, p, aj, bj):
    return solution_for_A(r, p, aj, bj, fa(r, p), fb(r, p))


def dAr_no_k(r, p, aj, bj):
    return solution_for_A(r, p, aj, bj, dfa(r, p), dfb(r, p))


# --- Testbench ------------------------------------------------------


theta = np.arange(0, 2 * np.pi, np.pi/4)
r = np.arange(theta.shape[0])
p = 2

print("r0_to_xy(r, theta): \n", r0_to_xy(r, theta), "\n")
print("fa(r, p): \n", fa(r, p), "\n")


















