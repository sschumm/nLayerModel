# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 08:06:20 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0

# -- iron -- | | -- air -- |k| -- air --


a_1 = lambda p, r1, r2, K, mu_1: (mu_0*mu_1) / (mu_0+mu_1) * (K/p) * r2**(-p+1)
b_1 = 0.
a_2 = lambda p, r1, r2, K, mu_1: (mu_0*K) / (2*p) * r2**(-p+1)
b_2 = lambda p, r1, r2, K, mu_1: (a_1(p, r1, r2, K, mu_1) - a_2(p, r1, r2, K, mu_1)) * r1**(2*p)
a_3 = 0.
b_3 = lambda p, r1, r2, K, mu_1: a_2(p, r1, r2, K, mu_1) * r2**(2*p) + b_2(p, r1, r2, K, mu_1)
# =============================================================================
# b_2 = lambda p, r1, r2, K, mu_1: (K/p) * ((mu_0*mu_1)/(mu_0+mu_1) - (mu_0/2)) * r1**(2*p) * r2**(-p+1)
# a_3 = 0.
# b_3 = lambda p, r1, r2, K, mu_1: (K/p) * ((mu_0/2) * r2**(p+1) + ((mu_0*mu_1)/(mu_0+mu_1) - (mu_0/2)) * r1**(2*p) * r2**(-p+1))
# =============================================================================


def analytic_solution(p, r1, r2, K, mu_r):
    mu_1 = mu_r * mu_0
    solution = [a_1(p, r1, r2, K, mu_1), 
                b_1, 
                a_2(p, r1, r2, K, mu_1), 
                b_2(p, r1, r2, K, mu_1), 
                a_3, 
                b_3(p, r1, r2, K, mu_1)]
    return np.asarray(solution)
    