# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 12:26:16 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import griddata
from matplotlib import cm
from matplotlib.colors import ListedColormap
from .layer import CurrentLoading, MagneticLayer, AirLayer

Reds = cm.get_cmap("Reds")
sampled_Reds = Reds(np.linspace(0,1,320))
myReds = ListedColormap(sampled_Reds[64:, :])

Winter = cm.get_cmap("winter_r")
sampled_Winter = Winter(np.linspace(0,1,360))
myWinter = ListedColormap(sampled_Winter[104:, :])

fgsz = 20


def streamplot_B(X, Y, U, V, radii, layers, axis=None):
    if axis is None:
        plt.figure(figsize=(fgsz, fgsz))
        ax = plt.subplot()
    else:
        ax = axis
        
    add_plot_dimensions(ax, radii)
    
    add_machine_dimensions(ax, layers)
                                
    
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
    lw = ( (0.35 * fgsz) / np.nanmax(gi) ) * gi + 0.5
    
    ax.streamplot(x,y,gu,gv, density=2, 
                  linewidth=lw, arrowsize= 0.1*fgsz, 
                  color=gi, cmap=myWinter)


def quiver_B(X, Y, U, V, radii, layers, axis=None):
    if axis is None:
        plt.figure(figsize=(fgsz, fgsz))
        ax = plt.subplot()
    else:
        ax = axis
        
    add_plot_dimensions(ax, radii)
    
    add_machine_dimensions(ax, layers)                           
    
    ax.quiver(X, Y, U, V, color="b")


def contour_A(X, Y, Z, radii, layers, axis=None):
    if axis is None:
        plt.figure(figsize=(fgsz, fgsz))
        ax = plt.subplot()
    else:
        ax = axis
        
    add_plot_dimensions(ax, radii)
    
    add_machine_dimensions(ax, layers)
    
    lvls = np.linspace(np.min(Z), np.max(Z), 10)
    CS = ax.contour(X, Y, Z, levels=lvls, linewidths=0.2*fgsz, cmap=myWinter)
    ax.clabel(CS, inline=True, fontsize=fgsz)
    print("")
    


def multi_figure(nx, ny, data):
    # plt.figure(figsize=(fgsz, fgsz))
    fig, axs = plt.subplots(ny, nx)
    fig.set_figheight(ny * fgsz)
    fig.set_figwidth(nx * fgsz)
    
    def ax_fct(ix, iy):
        if nx == 1:
            return axs[iy]
        elif ny == 1:
            return axs[ix]
        else:
            return axs[iy, ix]
    
    if nx*ny == len(data):
        for iy in range(ny):
            for ix in range(nx):
                fct = list(data)[iy+ix]
                if fct == "streamplot":
                    streamplot_B(*data[fct], ax_fct(ix, iy))
                elif fct == "quiver":
                    quiver_B(*data[fct], ax_fct(ix, iy))
                elif fct == "contour":
                    contour_A(*data[fct], ax_fct(ix, iy))
                else:
                    raise Exception("fct not implemented")
    else:
        raise Exception("data does not fit the number of subplots")
    


def add_machine_dimensions(ax, layers):
    for lay in reversed(layers.values()):
        edgecolor = "black"
        facecolor = "#e6e6e6" # "white"
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
                                 linewidth = 0.25 * fgsz, 
                                 alpha=alpha))
        

def add_plot_dimensions(ax, radii):
    
    lim = max(radii)*1.1

    ax.set_aspect( 1 )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    ax.tick_params(axis='both', which='major', labelsize=fgsz)
    ax.tick_params(axis='both', which='minor', labelsize=int(0.8 * fgsz))

 


# deprecated
# =============================================================================
# def plot(X, Y, U, V, radii, layers, style="quiver"):
#     
#     lim = max(radii)*1.1
#     
#     plt.figure(figsize=(fgsz, fgsz))
#     ax = plt.subplot()
# 
#     ax.set_aspect( 1 )
#     ax.set_xlim(-lim, lim)
#     ax.set_ylim(-lim, lim)
#     
#     
#     add_machine_dimensions(ax, layers)
#                                 
#     
#     x = np.linspace(X.min(), X.max(), 500)
#     y = np.linspace(Y.min(), Y.max(), 500)
#     xi, yi = np.meshgrid(x,y)
#     intensity = np.sqrt(U**2, V**2)
#     
#     px = X.flatten()
#     py = Y.flatten()
#     pu = U.flatten()
#     pv = V.flatten()
#     pi = intensity.flatten()
#     gu = griddata((px,py), pu, (xi,yi))
#     gv = griddata((px,py), pv, (xi,yi))
#     gi = griddata((px,py), pi, (xi, yi))
#     lw = (3.5 / np.nanmax(gi)) * gi + 0.5
#     
#     
#     
#     # ax.streamplot(X, Y, U, V, density=1.4, linewidth=None, color="black")
#     if style == "quiver":
#         ax.quiver(X, Y, U, V, color="b")
#     elif style == "streamplot":
#         ax.streamplot(x,y,gu,gv, density=2, linewidth=lw, color=gi, cmap=myWinter)
#     elif style == "all":
#         ax.quiver(X, Y, U, V, color="b")
#         ax.streamplot(x,y,gu,gv, density=3, linewidth=1, color="b", cmap=plt.cm.jet)
#     else:
#         print("...")
# 
# =============================================================================
