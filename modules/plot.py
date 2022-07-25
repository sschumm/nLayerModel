# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 11:56:57 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 

from matplotlib import cm
from matplotlib.colors import ListedColormap

from scipy.interpolate import griddata
from scipy.constants import pi
from .layer import MagneticLayer, CurrentLoading
from .model import Model

Winter = cm.get_cmap("winter_r")
sampled_Winter = Winter(np.linspace(0,1,360))
myWinter = ListedColormap(sampled_Winter[104:, :])

res = 256
n = 1
sampled_Winter = Winter(np.linspace(0,1,2*res))
# sampled_Winter[:res, -1] = np.flip(np.logspace(n, 0, res))
# sampled_Winter[res:, -1] = np.logspace(n, 0, res)

sampled_Winter[:res, -1] = np.flip(np.linspace(0, 1, res))**n
sampled_Winter[res:, -1] = np.linspace(0, 1, res)**n
contourWinter = ListedColormap(sampled_Winter)


class Plot():
    
    def __init__(self, model: Model, fgsz = 20):
        self.m = model
        
        self.fgsz = fgsz
        self.lim = max(model.radii) * 1.2
        self.r_max = np.sqrt(2 * self.lim**2)
        
        
    def _set_plot_dims(self, ax):
        ax.set_aspect(1)
        ax.set_xlim(-self.lim, self.lim)
        ax.set_ylim(-self.lim, self.lim)
    
        ax.tick_params(axis='both', which='major', labelsize=self.fgsz)
        ax.tick_params(axis='both', which='minor', labelsize=int(0.8 * self.fgsz))
        
    
    def _set_machine_dims(self, ax):
        for layer in reversed(self.m.layers):
            edgecolor = "black"
            facecolor = "#e6e6e6" # "white"
            if isinstance(layer, CurrentLoading):
                edgecolor = "red"
            elif isinstance(layer, MagneticLayer):
                facecolor = "#a6a6a6"
            else: 
                pass
            ax.add_artist(plt.Circle((0, 0), 
                                     layer.r, 
                                     fill = True, 
                                     edgecolor = edgecolor,
                                     facecolor = facecolor, 
                                     linestyle = "-",
                                     linewidth = 0.25 * self.fgsz, 
                                     alpha=0.7))

    
        
    def streamplot(self, dr, dt):
        plt.figure(figsize=(self.fgsz, self.fgsz))
        ax = plt.subplot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt)
        
        print("INFO: computing streamplot...")
        
        X, Y, U, V = self.m.get_B_data(r, t)
        
        x = np.linspace(X.min(), X.max(), self.fgsz * 50)
        y = np.linspace(Y.min(), Y.max(), self.fgsz * 50)
        xi, yi = np.meshgrid(x,y)
        intensity = np.sqrt(U**2 + V**2)
        
        px = X.flatten()
        py = Y.flatten()
        pu = U.flatten()
        pv = V.flatten()
        i = intensity.flatten()
        gu = griddata((px,py), pu, (xi,yi))
        gv = griddata((px,py), pv, (xi,yi))
        gi = griddata((px,py), i, (xi, yi))
        lw = ( (0.38 * self.lim * self.fgsz) / np.nanmax(gi) ) * gi + 0.2
        
        ax.streamplot(x,y,gu,gv,  
                      linewidth=lw, 
                      density= 3 * self.lim, 
                      arrowsize= 0.1*self.fgsz, 
                      color=gi, cmap=myWinter)
        print("INFO: finished streamplot.")
    
    
    def quiver():
        
        pass
        
    
    def contour(self, dr, dt, **kwargs):
        fig = plt.figure(figsize=(self.fgsz, self.fgsz))
        ax = plt.subplot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        lvls = kwargs.get("lvls", 50)
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt) 
        
        print("INFO: computing contour...")
        
        X, Y, Az = self.m.get_A_xy_data(r, t)
        # Az = Az / np.nanmax(Az)

        
        cs = ax.contourf(X, Y, Az,
                         levels=lvls,
                         vmin=np.nanmin(Az),
                         vmax=np.nanmax(Az),
                         cmap=contourWinter)
        #ax.clabel(cs, inline=True, fontsize=self.fgsz)
        fig.colorbar(cs, shrink = 0.8)
        
        print("INFO: finished contour.")
        
    
    
    
        
    
        