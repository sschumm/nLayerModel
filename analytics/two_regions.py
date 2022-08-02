# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 06:59:55 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0

# -- air -- |k| -- air --


a_1 = lambda p, r, K: 0.5 * mu_0 * (K / p) * r**-(p-1)
b_1 = 0.
a_2 = 0.
b_2 = lambda p, r, K: 0.5 * mu_0 * (K / p) * r**(p+1)


def analytic_solution(p, r, K):
    solution = [a_1(p, r, K), b_1, a_2, b_2(p, r, K)]
    return np.asarray(solution)
    