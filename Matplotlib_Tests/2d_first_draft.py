# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 22:46:45 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 

r = np.array([1, 2, 3, 4])
lim = np.max(r) * 1.2

fig, ax = plt.subplots(1)
ax.set_aspect( 1 )
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

for radius in r:
    ax.add_artist(plt.Circle((0, 0), 
                             radius = radius, 
                             fill = False, 
                             color = "black", 
                             linestyle = "--"))
 
