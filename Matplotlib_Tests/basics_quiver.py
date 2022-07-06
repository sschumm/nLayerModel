# -*- coding: utf-8 -*-
"""
Basics on plotting vector field with quiver()

Created on Wed Jul  6 14:53:27 2022

@author: svens
"""


#%%
import numpy as np
import matplotlib.pyplot as plt

#%% Single Arrow

X = [0]
Y = [0]
U = [2]  
V = [1]  

plt.title("some Title")
plt.xlabel("x - Label")
plt.ylabel("y - Label")
plt.grid(True)  
plt.xlim(-2, 5)
plt.ylim(-2, 2.5)


plt.quiver(X, Y, U, V, color='green', units='xy', scale=1)
plt.show()


#%% Basic Field

x,y = np.meshgrid(np.linspace(-5,5,10),np.linspace(-5,5,10))

u = 1
v = -1

plt.title("some Title")
plt.xlabel("x - Label")
plt.ylabel("y - Label")
plt.grid(True)
plt.xlim(-6, 6)
plt.ylim(-6, 6)

# =============================================================================
# # use ctrl+4 to comment out a whole section
# plt.annotate('center', xy=(0, 0), xytext=(2.8, 1.5),
#              arrowprops=dict(facecolor='red', shrink=0.05),
#              )
# =============================================================================

plt.quiver(x,y,u,v, color="blue")
plt.show()


#%% Field with Function

x,y = np.meshgrid(np.linspace(-5,5,10),np.linspace(-5,5,10))

u = x/np.sqrt(x**2 + y**2)
v = y/np.sqrt(x**2 + y**2)

plt.title("some Title")
plt.xlabel("x - Label")
plt.ylabel("y - Label")
plt.grid(True)
plt.xlim(-6, 6)
plt.ylim(-6, 6)

plt.quiver(x,y,u,v, color="red")
plt.show()


#%% Field with varying vector lengths

x, y = np.meshgrid(np.linspace(-5, 5, 10), 
                   np.linspace(-5, 5, 10))
  
u = -y/np.sqrt(x**2 + y**2)
v = x/(x**2 + y**2)
  
plt.title("some Title")
plt.xlabel("x - Label")
plt.ylabel("y - Label")
plt.grid(True)
lim = 7
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)

plt.quiver(x, y, u, v, color='g')
plt.show()

