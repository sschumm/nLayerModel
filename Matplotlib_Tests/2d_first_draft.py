# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 22:46:45 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 

r = np.array([1, 2, 3, 4])
c = ["black", "blue", "red", "green"]
lim = np.max(r) * 1.2

plt.figure(figsize=(10, 10))
ax = plt.subplot()
ax.set_aspect( 1 )
#ax.set_xlim(-lim, lim)
#ax.set_ylim(-lim, lim)

for i in range(3, -1, -1):
    ax.add_artist(plt.Circle((0, 0), 
                             radius = r[i], 
                             fill = True, 
                             facecolor = c[i], 
                             linestyle = "--",
                             alpha=0.2))
 
   
    
""" Adding some streamplot stuff """
# 1D arrays
x = np.arange(-5,5,0.1)
y = np.arange(-5,5,0.1)
  
# Meshgrid
X,Y = np.meshgrid(x,y)
  
# Assign vector directions
Ex = (X + 1)/((X+1)**2 + Y**2) - (X - 1)/((X-1)**2 + Y**2)
Ey = Y/((X+1)**2 + Y**2) - Y/((X-1)**2 + Y**2)

ax.streamplot(X,Y,Ex,Ey, density=1.4, linewidth=None, color="black")



""" Adding some quiver stuff """
x,y = np.meshgrid(np.linspace(-5,5,10),np.linspace(-5,5,10))

u = x/np.sqrt(x**2 + y**2)
v = y/np.sqrt(x**2 + y**2)
plt.quiver(x,y,u,v, color="red")


print("")














