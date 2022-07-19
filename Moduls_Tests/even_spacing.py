# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:08:06 2022

@author: svens
"""

import numpy as np
import matplotlib.pyplot as plt 
import somemath as sm
from scipy import interpolate

np.set_printoptions(suppress=True, linewidth=200, precision=2)

r_i = 0
r_o = 7.
lim = r_o * 1.1
detail = 8
much_detail = detail * 10
r = np.linspace(r_i, r_o, detail)
t = np.linspace(0, 2 * np.pi, detail)
R, T = np.meshgrid(r, t)

f1 = lambda r, t: r * t**2
f2 = lambda r, t: -r * (t- 3)


Br = f1(R, T)
Bt = f2(R, T)

X, Y = sm.r0_to_xy(R, T)
U, V = sm.r0_to_xy(Br, Bt)

space = np.linspace(-5., 5, 11)


def interpol(array):
    f = None
    f = interpolate.interp2d(R, T, array, kind="cubic")
    arr = np.round(f(np.linspace(r_i, r_o, much_detail), np.linspace(0, 2*np.pi, much_detail)), 1)
    idx = []
    for j in range(arr.shape[1]):
        # print("iter:", j)
        lst = []
        arr_l = arr[j].tolist()
        for i in space:
            if i in arr_l:
                lst.append(arr_l.index(i))
            else:
                # lst.append(0)
                pass
        while len(lst) < 6:
            lst.append(0)
        idx.append(lst)
    return np.asarray(idx)


def uv(array):
    g = None
    g = interpolate.interp2d(R, T, array, kind="cubic")
    arr = g(np.linspace(r_i, r_o, much_detail), np.linspace(0, 2*np.pi, much_detail))
    return arr 

Xi = interpol(X)
Yi = interpol(Y)
Ui = uv(U)
Vi = uv(V)

newU = np.zeros(Ui.shape)
for j in range(Xi.shape[0]):
    for i in range(Xi.shape[1]):
        newU[j, i] = Ui[Yi[j, i], Xi[j, i]]



# =============================================================================
# fy = None
# fx = interpolate.interp2d(R, T, X, kind="linear")
# X1 = np.round(fx(np.linspace(r_i, r_o, much_detail), np.linspace(0, 2*np.pi, much_detail)), 1)
# X_idx =[]
# for j in range(X1.shape[1]):
#     # print("iter:", j)
#     test = []
#     X1_l = X1[j].tolist()
#     for i in np.linspace(-5., 5., 11):
#         if i in X1_l:
#             test.append(X1_l.index(i))
#         else:
#             test.append(0)
#     X_idx.append(test)
# 
# fx = None
# fy = interpolate.interp2d(R, T, Y, kind="linear")
# Y1 = np.round(fy(np.linspace(r_i, r_o, much_detail), np.linspace(0, 2*np.pi, much_detail)), 1)
# 
# Y_idx =[]
# for j in range(Y1.shape[1]):
#     # print("iter:", j)
#     test1 = []
#     Y1_l = Y1[j].tolist()
#     for i in np.linspace(-5., 5., 11):
#         if i in Y1_l:
#             test1.append(Y1_l.index(i))
#         else:
#             test1.append(0)
#     Y_idx.append(test1)
# 
# =============================================================================
#%%
plt.figure(figsize=(10, 10))
ax = plt.subplot()
ax.set_aspect( 1 )
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

ax.quiver(X, Y, U, V, color="b")
print("")

# mx = np.max(X)
# mn = np.min(X)


# ax.streamplot(X, Y, U, V, density=1.4, linewidth=None, color="black")



#%%
lim = 7.
detail = 7
c = np.linspace(-lim, lim, detail)
X, Y = np.meshgrid(c, c)

getR = lambda x, y: np.sqrt(x**2, y**2)
get0 = lambda x, y: np.arctan2(y, x)

R = getR(X, Y)
THETA = get0(X, Y)

# Rb = R[0, 4:6]
# Rc = np.linspace(Rb[0], Rb[-1], detail)
# Rd = np.tile(Rc, (detail, 1))

R1, T1 = np.meshgrid(np.linspace(0.01,lim,detail), np.linspace(0, 2*np.pi, detail*2))

getX = lambda r, t: r * np.cos(t)
getY = lambda r, t: r * np.sin(t)

X1 = getX(R1, T1)
Y1 = getY(R1, T1)


"""
Try out:
    X, Y meshgrid
    get r
    get 0
    
    create new r matrix with linspace -> np.tile(arr, (Theta.shape[0], 1))
    compute Br, B0
    compute X, Y back
"""
fu = lambda x, y: x**2 + y
fv = lambda x, y: - x + y**2

U = fu(X, Y)
V = fv(X, Y)



plt.figure(figsize=(10, 10))
ax = plt.subplot()
ax.set_aspect( 1 )
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

ax.streamplot(X, Y, U, V, density=1.4, linewidth=None, color="black")
# ax.quiver(X, Y, U, V, color="b")
print("")