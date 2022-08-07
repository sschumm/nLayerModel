# -*- coding: utf-8 -*-
import numpy as np

from scipy.constants import mu_0, pi

# -- iron -- |k| -- air -- | | -- iron -- | | -- air --
def analytic_solution_K1(p, K1, r1, r2, r3, mu_r1, mu_r3):
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
    


# -- iron -- | | -- air -- |k| -- iron -- | | -- air --
def analytic_solution_K2(p, K2, r1, r2, r3, mu_r1, mu_r3):
    eta_r1 = (mu_r1 - 1) / (mu_r1 + 1)
    eta_r3 = (1 - mu_r3) / (1 + mu_r3)
    
    x = mu_r3 * (r2**(2*p) + r1**(2*p) * eta_r1)
    y = (r2**(2*p) - r1**(2*p) * eta_r1) * (r3**(2*p) + r2**(2*p) * eta_r3) * mu_r3 
    z = (r3**(2*p) - r2**(2*p) * eta_r3) * (r2**(2*p) + r1**(2*p) * eta_r1)
    
    # D
    a3 = ((eta_r3 * x) / (y + z)) * ((mu_0 * K2)/p) * r2**(p+1)
    
    # C
    a2 = (r3**(2*p) + eta_r3 * r2**(2*p)) / (eta_r3 * (r2**(2*p) + r1**(2*p) * eta_r1)) * a3
    
    # E
    b3 = (x / (y + z)) * ((mu_0 * K2)/p) * r2**(p+1) * r3**(2*p)
    
    # F
    b4 = a3 * r3**(2*p) + b3
    
    # A
    a1 = (mu_r1 / (1 + mu_r1)) * 2 * a2
    
    # B
    b2 = eta_r1 * r1**(2*p) * a2
    
    # zeros
    a4 = 0
    b1 = 0
    
    solution = [a1, 
                b1, 
                a2, 
                b2, 
                a3, 
                b3,
                a4,
                b4]
    return np.asarray(solution)
    
    
def analytic_torque_on_K1(p, l, K1, r1, a1_K2, alpha1, alpha2):
    const = l * K1 * p * a1_K2 * r1**(p+1)
    x = pi * np.sin(alpha1 - alpha2)
    y = (1 / (4*p)) * np.cos(4*p*pi + alpha1 + alpha2)
    z = (1 / (4*p)) * np.cos(alpha1 + alpha2)
    return np.abs(const * (x - y + z))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    