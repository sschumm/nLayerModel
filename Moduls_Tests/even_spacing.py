# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:08:06 2022

@author: svens
"""

import numpy as np
import matplotlib.pyplot as plt 

np.set_printoptions(suppress=True, linewidth=200, precision=2)

lim = 7.
detail = 7
c = np.linspace(-lim, lim, detail)
X, Y = np.meshgrid(c, c)

getR = lambda x, y: np.sqrt(x**2, y**2)
get0 = lambda x, y: np.arctan2(y, x)

R = getR(X, Y)
THETA = get0(X, Y)

#%%
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



#%%
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