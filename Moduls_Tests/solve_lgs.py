# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 08:39:42 2022

@author: svens
"""

import numpy as np


#%%

A = np.array([
    [0., 1.],
    [1., 0.]
    ])

b = np.array([0., 0.])

#%%

no_of_layers = 1
no_of_variables = (no_of_layers * 2 + 2) 

M = np.zeros((no_of_variables, no_of_variables))
N = np.zeros(no_of_variables)

M[ 0,  1] = 1
M[-1, -2] = 1


i = 1
r = 3
p = 2
inv_mu1 = 1 / (4 * np.pi * 10**(-7))
inv_mu2 = 1 / (4 * np.pi * 10**(-7))

# current loading
# B_r
M[  i, i-1:i+3] = np.array([-r**p, -r**-p, 
                             r**p,  r**-p])
# H_0
M[i+1, i-1:i+3] = np.array([-inv_mu1 * r**(p-1),  inv_mu1 * r**-(p+1),
                             inv_mu2 * r**(p-1), -inv_mu2 * r**-(p+1)]) * p
N[i+1] = 3


print("M:\n", M, "\n")
print("N:", N)




#%%
print("A:\n", A, "\n")
print("b:", b)

x = np.linalg.solve(A, b)
print("x:", x)