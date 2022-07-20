# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 12:26:16 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import griddata
from .layer import CurrentLoading, MagneticLayer, AirLayer
from matplotlib import cm
from matplotlib.colors import ListedColormap

Reds = cm.get_cmap("Reds")
sampled_Reds = Reds(np.linspace(0,1,320))
myReds = ListedColormap(sampled_Reds[64:, :])

Winter = cm.get_cmap("winter_r")
sampled_Winter = Winter(np.linspace(0,1,360))
myWinter = ListedColormap(sampled_Winter[104:, :])


def plot_contour(X, Y, Z, radii):
    
    lim = max(radii)*1.1
    
    plt.figure(figsize=(10, 10))
    ax = plt.subplot()
    ax.set_aspect( 1 )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    # ax = plt.subplots()
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, inline=True, fontsize=10)
    ax.set_title('Simplest default with labels')
    



def plot(X, Y, U, V, radii, layers, style="quiver"):
    
    # r = np.array([1, 2, 3, 4])
    # c = ["black", "blue", "red", "green"]
    # lim = np.max(r) * 1.2
    lim = max(radii)*1.1
    
    plt.figure(figsize=(10, 10))
    ax = plt.subplot()
    ax.set_aspect( 1 )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    
    for lay in reversed(layers.values()):
        edgecolor = "black"
        facecolor = "#d9d9d9" # "white"
        alpha = 0.7
        fill = True
        if isinstance(lay, CurrentLoading):
            edgecolor = "red"
        elif isinstance(lay, MagneticLayer):
            facecolor = "#a6a6a6"
        else: 
            pass
        ax.add_artist(plt.Circle((0, 0), 
                                 lay.r, 
                                 fill = fill, 
                                 edgecolor = edgecolor,
                                 facecolor = facecolor, 
                                 linestyle = "-",
                                 linewidth = 3, 
                                 alpha=alpha))
                                
    
# =============================================================================
#     for r in radii:
#         ax.add_artist(plt.Circle((0, 0), 
#                                  r, 
#                                  fill = False, 
#                                  edgecolor = "r",
#                                  facecolor = "b", 
#                                  linestyle = "-",
#                                  linewidth = 3, 
#                                  alpha=1))
# =============================================================================
       
        
    """ Adding some streamplot stuff """
    x = np.linspace(X.min(), X.max(), 500)
    y = np.linspace(Y.min(), Y.max(), 500)
    xi, yi = np.meshgrid(x,y)
    intensity = np.sqrt(U**2, V**2)
    
    px = X.flatten()
    py = Y.flatten()
    pu = U.flatten()
    pv = V.flatten()
    pi = intensity.flatten()
    gu = griddata((px,py), pu, (xi,yi))
    gv = griddata((px,py), pv, (xi,yi))
    gi = griddata((px,py), pi, (xi, yi))
    lw = (3.5 / np.nanmax(gi)) * gi + 0.5
    
    
    
    # ax.streamplot(X, Y, U, V, density=1.4, linewidth=None, color="black")
    if style == "quiver":
        ax.quiver(X, Y, U, V, color="b")
    elif style == "streamplot":
        ax.streamplot(x,y,gu,gv, density=2, linewidth=lw, color=gi, cmap=myWinter)
    elif style == "all":
        ax.quiver(X, Y, U, V, color="b")
        ax.streamplot(x,y,gu,gv, density=3, linewidth=1, color="b", cmap=plt.cm.jet)
    else:
        print("...")
