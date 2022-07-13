# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 10:13:17 2022

@author: svens
"""

import numpy as np

# =============================================================================
# A = np.array([
#     [3., 1, 5],
#     [7., 8, 3],
#     [4., 1, 2]
#     ])
# =============================================================================

n = 20
i = 10

A = np.random.random((n, n)) * 10

b = np.concatenate([np.array([5.]), np.zeros(n-1)]) # np.array([5., 0, 0])


B = np.stack([b for i in range(i)], axis=0)

# =============================================================================
# M = np.array([A, A, A])
# 
# 
# x = np.linalg.solve(M, B)
# =============================================================================


# =============================================================================
# A = np.array([
#     [np.arange(1,11), np.ones(10)],
#     [np.ones(10), np.ones(10)]
#     ])
# 
# B = np.array([
#     [4, 1],
#     [3, 2]
#     ])
# =============================================================================


E = np.ones((i, n, n))
E[:,0,0] = np.arange(1,i+1)

F = E * A


x = np.linalg.solve(F, B)



for iterate in range(i):
    print(np.allclose(np.dot(F[iterate], x[iterate]), b))



#%%

import numpy as np
import time

# =============================================================================
# i = 3
# A = np.array([
#     [3., 1, 5],
#     [7., 8, 3],
#     [4., 1, 2]
#     ])
# 
# b = np.array([0., 30, 5])
# =============================================================================

i = 10
A = np.random.randn(i,i)
b = np.random.randn(i)

# print("sol: \n", np.linalg.solve(A, b), "\n")

r = 1000
(rx, ry) = 0, 2

R = np.ones((r, i, i))
R[:, rx, ry] = np.linspace(1, 2, r)

A1 = R * A
B1 = np.stack([b for j in range(r)], axis=0)

# print("sol1: \n", np.linalg.solve(A1, B1), "\n")

mu = 1000
(mur, mux, muy) = 3, 1, 1

MU = np.ones((mu, r, i, i))
MU[:, mur, mux, muy] = np.linspace(1,4,mu)

A2 = MU * A1
B2 = np.stack([B1 for j in range(mu)], axis=0)

# print("sol2: \n", np.linalg.solve(A2, B2), "\n")

#%%
start = time.process_time_ns()

x = np.linalg.solve(A2, B2)

end = time.process_time_ns()
runtime = end-start
print("Runtime:", runtime/1000, "microseconds")





























