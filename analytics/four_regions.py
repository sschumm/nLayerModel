# -*- coding: utf-8 -*-
import numpy as np

from scipy.constants import mu_0

# -- iron -- |k| -- air -- | | -- iron -- | | -- air --



def analytic_solution(p, K1, r1, r2, r3, mu_r1, mu_r3):
    eta_r3 = (1-mu_r3) / (1+mu_r3)
    eta_r1 = (mu_r1-1) / (mu_r1+1)
    lambda1= (mu_r1 / (mu_r1+1)) * ((mu_0*K1)/p)
    gamma = mu_r3 * (r2**(2*p)*eta_r3 + r3**(2*p)) / (r3**(2*p) - r2**(2*p) * eta_r3)
    
    # 1.
    a2 = lambda1 * r1**(p+1) * (1 - gamma) / (gamma*(eta_r1*r1**(2*p) - r2**(2*p)) - (eta_r1*r1**(2*p) + r2**(2*p)))
    
    # III.
    b2 = eta_r1 * r1**(2*p) * a2 + lambda1*r1**(p+1)
    
    # D''
    a3 = ((mu_r3 * eta_r3) / (r3**(2*p) - eta_r3*r2**(2*p))) * (a2 * (eta_r1*r1**(2*p) - r2**(2*p)) + lambda1*r1**(p+1))
    
    # II.
    b3 = r3**(2*p) * ((mu_r3) / (r3**(2*p) - eta_r3*r2**(2*p))) * (a2 * (eta_r1*r1**(2*p) - r2**(2*p)) + lambda1*r1**(p+1))
    
    # B'
    a1 = (mu_r1/(mu_r1+1)) * 2*a2 + lambda1 * r1**(-p+1)
    
    # E
    b4 = a3 * r3**(2*p) + b3
    
    # zeros
    b1 = 0
    a4 = 0
    
    solution = [a1, 
                b1, 
                a2, 
                b2, 
                a3, 
                b3,
                a4,
                b4]
    return np.asarray(solution)
    
