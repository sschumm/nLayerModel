# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:36:32 2022

@author: svens
"""

import numpy as np 
import tikzplotlib as tp
import matplotlib.pyplot as plt 

from matplotlib import patches
from scipy.constants import pi, mu_0
from .colors import border_default, border_current, layer_default, layer_magnetic
from .colors import colorAccent, colorRed
from ..layer import MagneticLayer, CurrentLoading
from ..model import Model



class RadialPlot():
    
    
    def __init__(self, model: Model, fgsz = 20):
        self.m = model
        
        self.fgsz_x = fgsz
        self.fgsz_y = fgsz
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
                v_color = colorRed

            ax.axvline(layer.r, color=v_color, linewidth=2)
                    
            
    def _unpack_kwargs(self, ax, **kwargs):
        dr = kwargs.get("dr", self.fgsz_x * 50)
        r = np.linspace(0, self.lim, dr)
        t = kwargs.get("t", 0)/self.m.p
        
        ax.grid(True)
        title = kwargs.get("title", "")
        ax.set_title(title, fontsize=self.fgsz_x*1.4)
        tikz = kwargs.get("tikz", False)

        return r, t, tikz
                    
        
    def plot_radial_Az(self, **kwargs):
        fig, ax = self._set_up_plot()
        r, t, tikz = self._unpack_kwargs(ax, **kwargs)
        d = self.m.get_A_data(r, t)
        
        
        this_Az = d.Az.flatten() *1e-3
        self.ymax = np.nanmax(this_Az)
        self.ymin = np.nanmin(this_Az)
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        fig.tight_layout()
        ax.plot(r, this_Az, linewidth=3, color=colorAccent) 
        
        if tikz:
            tp.clean_figure()
            tp.save("msm_radial_Ht.tex")
        

    def plot_radial_Br(self, **kwargs):
        fig, ax = self._set_up_plot()
        r, t, tikz = self._unpack_kwargs(ax, **kwargs)
        d = self.m.get_B_data(r, t)
        
        this_Br = d.Br.flatten()
        self.ymax = np.nanmax(this_Br)
        self.ymin = np.nanmin(this_Br)        
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        fig.tight_layout()
        ax.plot(r, this_Br, linewidth=3, color=colorAccent) 
        
        if tikz:
            tp.clean_figure()
            tp.save("msm_radial_Ht.tex") 
    
    
    def plot_radial_Ht(self, **kwargs):
        fig, ax = self._set_up_plot()
        r, t, tikz = self._unpack_kwargs(ax, **kwargs)
        d = self.m.get_H_data(r, t)
        

        this_Ht = d.Ht.flatten() *1e-3
        self.ymax = np.nanmax(this_Ht)
        self.ymin = np.nanmin(this_Ht)
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        fig.tight_layout()
        ax.plot(r, this_Ht, linewidth=3, color=colorAccent) 
        
        if tikz:
            tp.clean_figure()
            tp.save("msm_radial_Ht.tex")
        


    

class RadialMsizePlot(RadialPlot):
    def __init__(self, model: Model, fgsz = 20):
        super().__init__(model, fgsz)
               
        
    def _set_up_plot(self):
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [9, 1]})
        fig.set_figheight(self.fgsz_y + self.machinesz)
        fig.set_figwidth(self.fgsz_x)
        fig.tight_layout() 
        
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
            v_color = border_default
            h_color = layer_default
            if isinstance(layer, CurrentLoading):
                if not layer.mu == mu_0:
                    h_color = layer_magnetic
                v_color = border_current
            if isinstance(layer, MagneticLayer):
                h_color = layer_magnetic
                
            ax.add_patch(patches.Rectangle((layer.r, 0), 
                                           width=1,
                                           height=layer.r,
                                           angle=90,
                                           facecolor=h_color
                                           )
                         )
            ax.axvline(layer.r, color=v_color, linewidth=4)



  
  
class RadialMultiPlot(RadialMsizePlot):
    
    def __init__(self, model: Model, fgsz = 20):
        super().__init__(model, fgsz)
        
        self.fig = None
        self.axs = None
        self.count_y = 0
        self.multi = False
        
        self.Az_details = dict()
        self.Br_details = dict()
        self.Ht_details = dict()
        
        
    def _set_up_plot(self):
        if not self.multi:
            fig = plt.figure(figsize=(self.fgsz_x, self.fgsz_y))
            ax = plt.subplot()        
            return fig, ax
        else:           
            return self.fig, self.axs[self.count_y]
    
    
    def set_Az_details(self, **kwargs):
        self.Az_details = self.Az_details | kwargs
        
    
    def set_Br_details(self, **kwargs):
        self.Br_details = self.Br_details | kwargs
        
        
    def set_Ht_details(self, **kwargs):
        self.Ht_details = self.Ht_details | kwargs
    
    
    def multiplot(self, quantities: list = ["Br", "Ht"], **kwargs):
        self.multi = True
        self.count_y = 0
        
        # ------- get information on what to draw
        ny = len(quantities)
        self.fgsz_y = self.fgsz_x * 0.8/ny
        machine_dims = kwargs.get("machine_dims", True)
        gspecs = {"height_ratios": [self.fgsz_y] * ny}
        height = self.fgsz_y * ny        
        if machine_dims:
            ny += 1
            height += 1
            gspecs["height_ratios"].append(1)
        
        
        # ------- build figure -------
        self.fig, self.axs = plt.subplots(ny, 1, gridspec_kw=gspecs)
        self.fig.set_figheight(height)
        self.fig.set_figwidth(self.fgsz_x)
        if machine_dims:
            self._add_machine_dims(self.axs[-1])
        
        
        # ------- handle kwargs ------- 
        a = kwargs.get("angle", 0)
        if "title" in kwargs:
            self.fig.suptitle(kwargs["title"], fontsize=self.fgsz_x * 2)
        
        
        
        # ------- draw the different plots -------
        for q in quantities:
            if q == "Az":
                self.Az_details["angle"] = a
                self.plot_radial_Az(**self.Az_details)
            elif q == "Br":
                self.Br_details["angle"] = a
                self.plot_radial_Br(**self.Br_details)
            elif q == "Ht":
                self.Ht_details["angle"] = a
                self.plot_radial_Ht(**self.Ht_details)
            else:
                print(f"INFO: Plot Quantity {q} is unknown.")
                
            self.count_y += 1

        self.fig.tight_layout()
    
    
    
    