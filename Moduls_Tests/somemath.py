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
    if np.all(r):
        return r**-p
    else:
        raise ZeroDivisionError(r)
    

def dfa(r, p):
    return p * r**(p-1)


def dfb(r, p):
    if np.all(r):
        return -p * r**-(p+1)
    else:
        raise ZeroDivisionError(r)


# -----------------------

def solution_for_A(r, p, aj, bj, fa, fb):
    # Note: Huge bug potential when trying to use arrays for r and for aj/bj
    if type(aj) is not type(bj):
        raise Exception("(sschumm): aj and bj have to be of the same datatype")
    if type(aj) is np.ndarray and type(r) is np.ndarray:
        raise Exception("(sschumm): this function is only designed to accept multiple r's or multiple solutions aj and bj")
    
    return np.add(aj * fa(r, p), bj * fb(r, p))
    


def Ar_no_k(r, p, aj, bj):
    return solution_for_A(r, p, aj, bj, fa, fb)


def dAr_no_k(r, p, aj, bj):
    return solution_for_A(r, p, aj, bj, dfa, dfb)


# --- Testbench ------------------------------------------------------


theta = 2 # np.arange(0, 2 * np.pi, np.pi/4, dtype=np.longdouble)
r = 3 # np.linspace(1, 8, theta.shape[0], dtype=np.longdouble)
p = 2
aj = 3 # np.linspace(1,5, 5) # 4
bj = 3  #np.linspace(5,1, 5) # 4

print("r0_to_xy(r, theta): \n", r0_to_xy(r, theta), "\n")
#
print("fa(r, p): \n", fa(r, p), "\n")
print("fb(r, p): \n", fb(r, p), "\n")
print("dfa(r, p): \n", dfa(r, p), "\n")
print("dfb(r, p): \n", dfb(r, p), "\n")
#
#print("solution_for_A(r, p, aj, bj, fa, fb): \n", solution_for_A(r, p, aj, bj, fa, fb), "\n")
print("Ar_no_k(r, p, aj, bj): \n", Ar_no_k(r, p, aj, bj), "\n")
print("dAr_no_k(r, p, aj, bj): \n", dAr_no_k(r, p, aj, bj), "\n")




















