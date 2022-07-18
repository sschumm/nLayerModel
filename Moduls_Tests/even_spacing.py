# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 18:08:06 2022

@author: svens
"""

import numpy as np
import matplotlib.pyplot as plt 


lim = 7.
c = np.linspace(-lim, lim, 16)
X, Y = np.meshgrid(c, c)

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