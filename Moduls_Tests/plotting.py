# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 12:26:16 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import griddata
from .layer import CurrentLoading, MagneticLayer, AirLayer


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
        facecolor = "white"
        alpha = 0.7
        fill = True
        if isinstance(lay, CurrentLoading):
            edgecolor = "blue"
        elif isinstance(lay, MagneticLayer):
            facecolor = "#cccccc"
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
    
    
    
    # ax.streamplot(X, Y, U, V, density=1.4, linewidth=None, color="black")
    if style == "quiver":
        ax.quiver(X, Y, U, V, color="b")
    elif style == "streamplot":
        ax.streamplot(x,y,gu,gv, density=3, linewidth=1, color=gi, cmap=plt.cm.Reds)
    elif style == "all":
        ax.quiver(X, Y, U, V, color="b")
        ax.streamplot(x,y,gu,gv, density=3, linewidth=1, color="b", cmap=plt.cm.jet)
    else:
        print("...")
