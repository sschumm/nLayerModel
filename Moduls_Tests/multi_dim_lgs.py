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

i = 2
A = np.random.randn(i,i)
b = np.random.randn(i)

# print("sol: \n", np.linalg.solve(A, b), "\n")

r = 4
(rx, ry) = 0, 1

R = np.ones((r, i, i))
R[:, rx, ry] = np.linspace(1, 2, r)

A1 = R * A
B1 = np.stack([b for j in range(r)], axis=0)

# print("sol1: \n", np.linalg.solve(A1, B1), "\n")

mu = 3
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


#%%
# deprecated
# =============================================================================
# i = 5
# A = np.random.randn(i,i)
# B = np.random.randn(i)
# var1 = 4
# rng = np.linspace(1,13,var1)
# coorlst = [get_coor(0, 1), get_coor(2, 2)]
# 
# R = r_matrix(rng, A.shape, A.shape, coorlst)
# C, D = find_matrices(R, A, B)
# 
# for j in range(C.shape[0]):
#     print(C[j] - A)
#     
# var2 = 3    
# rng2 = np.linspace(1,21,var2)
# coorlst2 = [all_slice(get_coor(2, 0))] # [add_coor(0, get_coor(2, 0))] 
# R2 = r_matrix(rng2, C.shape, A.shape, coorlst2)
# E, F = find_matrices(R2, C, D)
# 
# 
# for j in range(E.shape[0]):
#     print(E[j] - A)
# 
# np.set_printoptions(suppress=True, linewidth=200, precision=4)
# x = np.linalg.solve(E, F)
# print("x: \n", x, "\n")
# print("x: \n", x[0], "\n")
# print("x: \n", x[0, 0], "\n")
# =============================================================================



#%%

import numpy as np

def get_slice(y, x):
    return (slice(y, y+1), slice(x, x+1))

def add_slice(z, slc):
    return (slice(z, z+1),) + slc

def all_slice(slc):
    if isinstance(slc, list):
        return [(slice(0, None),) + j for j in slc]
    return (slice(0, None),) + slc

def get_coor(y, x):
    return (y, x)

def add_coor(z, coor):
    return (z,) + coor

def find_matrices(R, A, B):
    A1 = R * A
    B1 = np.stack([B for j in range(R.shape[0])], axis=0)
    return A1, B1


def r_matrix(var_arrays, shp_2, shp_1, idx):
    arr = var_arrays[idx]["arr"]
    
    R = np.ones(arr.shape + shp_2)
    
    if idx != 0:
        
        arr1 = None
        for j in range(idx):
            
            
            
            mult = R.shape[-(j+4)]
            arr1 = np.stack([var_arrays[j]["arr"] for i in range(mult)], axis=-1)
        ##
        # mult = R.shape[diff]
        # arr2 = np.stack([arr for i in range(mult)], axis=1)
    else:
        arr1 = arr
    
    for coor in var_arrays[idx]["loc"]:
        this_slc = all_slice(coor)
        R[this_slc] = arr1
    return R


def vary(A, B, var_arrays):
    dims = A.shape
    # variations = len(var_arrays)
    A1 = A
    B1 = B
    
    ##
    for idx, dic in enumerate(var_arrays):
        #coorlst = dic["loc"]
        for i in range(idx):
            dic["loc"] = all_slice(dic["loc"])
    ##
    
    
    for idx, dic in enumerate(var_arrays):
        #rng = dic["arr"]
        #coorlst = dic["loc"]
        
        ##
        # for i in range(idx):
        #     coorlst = all_slice(coorlst)
        ##
        R = r_matrix(var_arrays, A1.shape, dims, idx)
        # print(R)
        # print(R.shape)
        A1, B1 = find_matrices(R, A1, B1)
    x = np.linalg.solve(A1, B1)
    return x, A1, B1
        
    
i = 3
A = np.random.randn(i,i)
B = np.random.randn(i)

var1, var2, var3 = 2, 2, 3
rng1, rng2, rng3 = np.linspace(1,13,var1), np.linspace(1,21,var2), np.ones(var3) * 999

var_arrays = []
var_arrays.append({
    "arr": rng1, 
    "loc": [get_coor(0,1), get_coor(2,2)]
    })
var_arrays.append({
    "arr": rng2,
    "loc": [get_coor(2, 0)]
    })
# var_arrays.append({
#     "arr": rng3,
#     "loc": [get_coor(0, 0)]
#     })

x, At, Bt = vary(A, B, var_arrays)



AA = np.copy(A)
x_test = np.empty((var2,var1,i))        
for par2 in range(var2):
    AA = np.copy(A)
    AA[get_coor(2,0)] *= rng2[par2]
    
    for par in range(var1):
        AAA = np.copy(AA)
        AAA[var_arrays[0]["loc"][0]] *= rng1[par]
        AAA[var_arrays[0]["loc"][1]] *= rng1[par]
        
        x1 = np.linalg.solve(AAA, B)
        x_test[par2, par] = x1
        #print(x1)
    #print("")

np.set_printoptions(suppress=True, linewidth=200, precision=4)

## Goal
print("x: \n", x, "\n")

## Review
print("x1: \n")
print(x_test)
print("")

## Compare
print("Comparison:")
print(np.round(x - x_test, 5))

    
# print("At: \n", At, "\n")
# print("Bt: \n", Bt, "\n")




