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

A = np.array([
    [3., 1, 5],
    [7., 8, 3],
    [4., 1, 2]
    ])

b = np.array([0., 30, 5])#.reshape(3,1)
# b = np.array([[0.], [30], [5]])

print("sol: \n", np.linalg.solve(A, b), "\n")


r = 5
(rx, ry) = 0, 2

R = np.ones((r, 3, 3))
R[:, rx, ry] = np.linspace(1, 2, r)

A1 = R * A
B1 = np.stack([b for i in range(r)], axis=0)

print("sol1: \n", np.linalg.solve(A1, B1), "\n")

mu = 4
(mur, mux, muy) = 3, 1, 1

MU = np.ones((mu, r, 3, 3))
MU[:, mur, mux, muy] = np.linspace(1,4,mu)

A2 = MU * A1
B2 = np.stack([B1 for i in range(mu)], axis=0)

# =============================================================================
# print(A.shape)
# print(b.shape)
# 
# print(A1.shape)
# print(B1.shape)
# 
# print(A2.shape)
# print(B2.shape)
# =============================================================================

print("sol2: \n", np.linalg.solve(A2, B2), "\n")


































