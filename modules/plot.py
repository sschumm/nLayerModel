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
from .layer import MagneticLayer, CurrentLoading
from .model import Model

Winter = cm.get_cmap("winter_r")
sampled_Winter = Winter(np.linspace(0,1,360))
myWinter = ListedColormap(sampled_Winter[104:, :])


class Plot():
    
    def __init__(self, model: Model, fgsz = 20):
        self.m = model
        
        self.fgsz = fgsz
        self.lim = max(model.radii) * 1.2
        
        
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

    
        
    def streamplot(self):
        plt.figure(figsize=(self.fgsz, self.fgsz))
        ax = plt.subplot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
    
    
    def quiver():
        
        pass
        
    
    def contour():
        pass
    
    
    
        
    
        