# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:36:32 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 

from matplotlib import patches

from scipy.constants import pi, mu_0
from ..layer import MagneticLayer, CurrentLoading
from ..model import Model



class RadialPlot():
    
    
    def __init__(self, model: Model, fgsz = 20):
        self.m = model
        
        self.fgsz_x = fgsz
        self.fgsz_y = fgsz * 0.6
        self.machinesz = self.fgsz_y * 0.1
        self.lim = max(model.radii) * 1.2
        self.ymax = 1
        self.ymin = -1
        
        
    def _set_up_plot(self):
        fig = plt.figure(figsize=(self.fgsz_x, self.fgsz_y))
        # plt.style.use('_mpl-gallery')
        ax = plt.subplot()        
        return fig, ax
        
    
    def _set_plot_dims(self, ax):
        ax.set_xlim(0, self.lim)
        ax.set_ylim(self.ymin * 1.1, self.ymax * 1.1)
        
        ax.tick_params(axis='both', which='major', 
                       labelsize=self.fgsz_x, pad=self.fgsz_x)
        ax.tick_params(axis='both', which='minor', 
                       labelsize=int(0.8 * self.fgsz_x), pad=self.fgsz_x)
        ax.set_ylabel
        
        
    def _set_machine_dims(self, ax):
        for layer in reversed(self.m.layers):
            v_color = "black"
            if isinstance(layer, CurrentLoading):
                v_color = "red"

            ax.axvline(layer.r, color=v_color, linewidth=4)
                    
            
    def _unpack_kwargs(self, kwargs):
        dr = kwargs.get("dr", self.fgsz_x * 50)
        dt = kwargs.get("dt", self.fgsz_x * 50)
        a = kwargs.get("angle", 0)
        r = np.linspace(0, self.lim, dr)
        t = np.linspace(0, 2*pi, dt)
        return r, t, a
                    
        
    def plot_radial_Az(self, **kwargs):
        fig, ax = self._set_up_plot()
        r, t, a = self._unpack_kwargs(kwargs)
        _, _, Az, R, T = self.m.get_A_data(r, t)
        
        
        this_Az = Az[int(len(t) * a)]
        self.ymax = np.nanmax(this_Az)
        self.ymin = np.nanmin(this_Az)
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        ax.grid(True)
        ax.plot(r, this_Az)        
        

    def plot_radial_Br(self, **kwargs):
        fig, ax = self._set_up_plot()
        r, t, a = self._unpack_kwargs(kwargs)
        _, _, _, _, _, _, Br, _ = self.m.get_B_data(r, t)
        
        this_Br = Br[int(len(t) * a)]
        self.ymax = np.nanmax(this_Br)
        self.ymin = np.nanmin(this_Br)        
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        ax.plot(r, this_Br) 
    
    
    def plot_radial_Ht(self, **kwargs):
        fig, ax = self._set_up_plot()
        r, t, a = self._unpack_kwargs(kwargs)
        _, _, _, _, _, _, _, Ht = self.m.get_H_data(r, t)
        

        this_Ht = Ht[int(len(t) * a)]
        self.ymax = np.nanmax(this_Ht)
        self.ymin = np.nanmin(this_Ht)
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        ax.plot(r, this_Ht) 
        
        

    
class RadialMsizePlot(RadialPlot):
    def __init__(self, model: Model, fgsz = 20):
        super().__init__(model, fgsz)
               
        
    def _set_up_plot(self):
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [9, 1]})
        fig.set_figheight(self.fgsz_y + self.machinesz)
        fig.set_figwidth(self.fgsz_x)
        
        self._add_machine_dims(axs[1])
        
        return fig, axs[0]
    
    
    def _add_machine_dims(self, ax):
        
        ax.set_xlim(0, self.lim)
        ax.set_ylim(0,1)
        ax.grid(False)
        ax.tick_params(axis='both', which='both', 
                       left=False,
                       labelleft=False,
                       bottom=False,
                       labelbottom=False)
        
        for layer in reversed(self.m.layers):
            v_color = "black"
            h_color = "#e6e6e6" # "white"
            if isinstance(layer, CurrentLoading):
                if not layer.mu == mu_0:
                    h_color = "#a6a6a6"
                v_color = "red"
            if isinstance(layer, MagneticLayer):
                h_color = "#a6a6a6"
                
            ax.add_patch(patches.Rectangle((layer.r, 0), 
                                           width=1,
                                           height=layer.r,
                                           angle=90,
                                           facecolor=h_color
                                           )
                         )
            ax.axvline(layer.r, color=v_color, linewidth=4)

    
    
    
    
    
    
    