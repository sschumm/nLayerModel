# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:38:57 2022

@author: svens
"""
import numpy as np



def r0_to_x(r, theta):
    return r * np.cos(theta)


def r0_to_y(r, theta):
    return r * np.sin(theta)


def r0_to_xy(r, theta):
    return r0_to_x(r, theta), r0_to_y(r, theta)    
    

# -----------------------

def fa(p, r):
    return r**p


def fb(p, r):
    if np.all(r):
        return r**-p
    else:
        raise ZeroDivisionError(r)
    

def dfa(p, r):
    return p * r**(p-1)


def dfb(p, r):
    if np.all(r):
        return -p * r**-(p+1)
    else:
        raise ZeroDivisionError(r)


# -----------------------

def solution_for_A(p, r, aj, bj, fa, fb):
    # Note: Huge bug potential when trying to use arrays for r and for aj/bj
    if type(aj) is not type(bj):
        raise Exception("(sschumm): aj and bj have to be of the same datatype")
    if type(aj) is np.ndarray and type(r) is np.ndarray:
        raise Exception("(sschumm): this function is only designed to accept multiple r's or multiple solutions aj and bj")
    else:
        return np.add(aj * fa(p, r), bj * fb(p, r))
    


def Ar_no_k(p, r, aj, bj):
    return solution_for_A(p, r, aj, bj, fa, fb)


def dAr_no_k(p, r, aj, bj):
    return solution_for_A(p, r, aj, bj, dfa, dfb)


# --------------------------------------------------

def Az_no_k(p, r, theta, aj, bj):
    if np.all(r):
        return np.sin(p * theta) * Ar_no_k(p, r, aj, bj)
    else:
        raise ZeroDivisionError(r)    
    

def Br_no_k(p, r, theta, aj, bj):
    if np.all(r):
        return (p / r) * np.cos(p * theta) * Ar_no_k(p, r, aj, bj)
    else:
        raise ZeroDivisionError(r)
        
        
def B0_no_k(p, r, theta, aj, bj):
    if np.all(r):
        return - np.sin(p * theta) * dAr_no_k(p, r, aj, bj)
    else:
        raise ZeroDivisionError(r)
        
        
def Hr_no_k(p, r, theta, aj, bj, mu):
    if type(mu) is np.ndarray and (type(r) or type(theta) or type(aj) or type(bj)) is np.ndarray:
        raise Exception("(sschumm): this function is only designed to vary either aj/bj or r/theta or mu")
    if np.all(mu):
        return (p / (r * mu)) * np.cos(p * theta) * Ar_no_k(p, r, aj, bj)
    else:
        raise ZeroDivisionError(r)


def H0_no_k(p, r, theta, aj, bj, mu):
    if type(mu) is np.ndarray and (type(r) or type(theta) or type(aj) or type(bj)) is np.ndarray:
        raise Exception("(sschumm): this function is only designed to vary either aj/bj or r/theta or mu")
    if np.all(mu):
        return (1 / mu) * np.sin(p * theta) * dAr_no_k(p, r, aj, bj)
    else:
        raise ZeroDivisionError(r)
    


# --- Testbench ------------------------------------------------------


# =============================================================================
# p = 2
# theta = 2 * np.arange(0, 2 * np.pi, np.pi/9, dtype=np.longdouble)
# r     = 3 #* np.linspace(1, 8, theta.shape[0], dtype=np.longdouble)
# aj    = 3 #* np.linspace(1,5, 5) # 4
# bj    = 3 #* np.linspace(5,1, 5) # 4
# mu    = 4 * np.pi * 10e-6  #* np.linspace(1,10,10)
# 
# print("r0_to_xy(r, theta): \n", r0_to_xy(r, theta), "\n")
# 
# print("fa(p, r): \n", fa(p, r), "\n")
# print("fb(p, r): \n", fb(p, r), "\n")
# print("dfa(p, r): \n", dfa(p, r), "\n")
# print("dfb(p, r): \n", dfb(p, r), "\n")
# 
# print("Ar_no_k(p, r, aj, bj): \n", Ar_no_k(p, r, aj, bj), "\n")
# print("dAr_no_k(p, r, aj, bj): \n", dAr_no_k(p, r, aj, bj), "\n")
# 
# print("Az_no_k(p, r, theta, aj, bj): \n", Az_no_k(p, r, theta, aj, bj), "\n")
# print("Br_no_k(p, r, theta, aj, bj): \n", Br_no_k(p, r, theta, aj, bj), "\n")
# print("B0_no_k(p, r, theta, aj, bj): \n", B0_no_k(p, r, theta, aj, bj), "\n")
# print("Hr_no_k(p, r, theta, aj, bj, mu): \n", Hr_no_k(p, r, theta, aj, bj, mu), "\n")
# print("H0_no_k(p, r, theta, aj, bj, mu): \n", H0_no_k(p, r, theta, aj, bj, mu), "\n")
# 
# 
# =============================================================================














