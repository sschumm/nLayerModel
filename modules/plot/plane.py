# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:32:31 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 

from matplotlib import cm, patches
from matplotlib.colors import ListedColormap

from scipy.interpolate import griddata
from scipy.constants import pi, mu_0
from ..layer import MagneticLayer, CurrentLoading
from ..model import Model

Winter = cm.get_cmap("winter_r")
sampled_Winter = Winter(np.linspace(0,1,360))
myWinter = ListedColormap(sampled_Winter[104:, :])

res = 256
n = 0.65
sampled_Winter = Winter(np.linspace(0,1,2*res))
sampled_Winter[:res, -1] = np.flip(np.linspace(0., 1, res))**n
sampled_Winter[res:, -1] = np.linspace(0., 1, res)**n
contourWinter = ListedColormap(sampled_Winter)


class PlanePlot():
    
    def __init__(self, model: Model, fgsz = 20):
        self.m = model
        
        self.fgsz = fgsz
        self.lim = max(model.radii) * 1.2
        self.r_max = np.sqrt(2 * self.lim**2)
        
        
    def _set_up_plot(self):
        fig = plt.figure(figsize=(self.fgsz, self.fgsz))
        ax = plt.subplot()        
        return fig, ax
        
        
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
                if not layer.mu == mu_0:
                    facecolor = "#a6a6a6"
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
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt)
        
        print("INFO: computing streamplot...")
        
        X, Y, U, V, _, _, _, _ = self.m.get_B_data(r, t)
        
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
        lw = ( (0.38 * self.fgsz) / np.nanmax(gi) ) * gi + 0.2
        
        ax.streamplot(x,y,gu,gv,  
                      linewidth=lw, 
                      density= 3 * self.lim, 
                      arrowsize= 0.1*self.fgsz, 
                      color=gi, cmap=myWinter)
        print("INFO: finished streamplot.")
           
    
    def contour(self, dr, dt, **kwargs):
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        lvls = kwargs.get("lvls", 50)
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt) 
        
        print("INFO: computing contour...")
        
        X, Y, Az, _, _ = self.m.get_A_data(r, t)
                
        cs = ax.contourf(X, Y, Az,
                         levels=lvls,
                         vmin=np.nanmin(Az),
                         vmax=np.nanmax(Az),
                         cmap=contourWinter)
        #ax.clabel(cs, inline=True, fontsize=self.fgsz)
        # fig.colorbar(cs, shrink = 0.8)
        
        print("INFO: finished contour.")
        
    
    def quiver(self, dr, dt):
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt)
        
        print("INFO: computing quiver...")
        
        X, Y, U, V, _, _, _, _ = self.m.get_B_data(r, t)
        intensity = np.sqrt(U**2 + V**2)
    
        ax.quiver(X, Y, U, V, intensity,
                  scale=3*dt,
                  cmap=myWinter) # color="b")
        print("INFO: finished quiver.")
    

class PlaneDoublePlot(PlanePlot):
    
    def __init__(self, model: Model, fgsz = 20):
        super().__init__(model, fgsz)
        
        self.ny = 1
        self.nx = 2
        self.fig = None
        self.axs = None
        self.count_x = 0
        
        
    def _set_up_plot(self):
        if self.count_x == 0:
            self.fig, self.axs = plt.subplots(self.ny, self.nx)
            self.fig.set_figheight(self.ny * self.fgsz)
            self.fig.set_figwidth(self.nx * self.fgsz)
        
        this_count_x = self.count_x
        self.count_x += 1
        return self.fig, self.axs[this_count_x]
    
    
    def plot_BandA(self, dr, dt):
        self.streamplot(dr, dt)
        self.contour(dr, dt)
    
      
   
    
    