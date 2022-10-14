# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:32:31 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 

from scipy.interpolate import griddata
from scipy.constants import pi, mu_0
from .colors import strp_winter, cntr_winter, cntr_jet, cntr_rb, flux, flux_cutoff
from .colors import border_default, border_current, layer_default, layer_magnetic
from ..layer import MagneticLayer, CurrentLoading
from ..model import Model


class PlanePlot():
    
    def __init__(self, model: Model, fgsz = 50):
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
            edgecolor = border_default
            facecolor = layer_default
            if isinstance(layer, CurrentLoading):
                if not layer.mu == mu_0:
                    facecolor = layer_magnetic
                edgecolor = border_current
            elif isinstance(layer, MagneticLayer):
                facecolor = layer_magnetic
            else: 
                pass
            ax.add_artist(plt.Circle((0, 0), 
                                     layer.r, 
                                     fill = True, 
                                     edgecolor = edgecolor,
                                     facecolor = facecolor, 
                                     linestyle = "-",
                                     linewidth = 0.1 * self.fgsz, 
                                     alpha=0.7))

    
        
    def streamplot(self, dr, dt):
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt)
        
        print("INFO: computing streamplot...")
        
        d = self.m.get_B_data(r, t)
        
        x = np.linspace(d.X.min(), d.X.max(), self.fgsz * 50)
        y = np.linspace(d.Y.min(), d.Y.max(), self.fgsz * 50)
        xi, yi = np.meshgrid(x,y)
        intensity = np.sqrt(d.U**2 + d.V**2)
        
        px = d.X.flatten()
        py = d.Y.flatten()
        pu = d.U.flatten()
        pv = d.V.flatten()
        i = intensity.flatten()
        gu = griddata((px,py), pu, (xi,yi))
        gv = griddata((px,py), pv, (xi,yi))
        gi = griddata((px,py), i, (xi, yi))
        lw = ( (0.38 * self.fgsz) / np.nanmax(gi) ) * gi + 0.2
        
        ax.streamplot(x,y,gu,gv,  
                      linewidth=lw, 
                      density= 3 * self.lim, 
                      arrowsize= 0.1*self.fgsz, 
                      color=gi, cmap=strp_winter)
        print("INFO: finished streamplot.")
           
    
    def contour(self, dr, dt, **kwargs):
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        lvls = kwargs.get("lvls", 50)
        style = kwargs.get("style", "winter")
        pdf = kwargs.get("pdf", False)
        pdf_dpi = kwargs.get("pdf_dpi", 300)
        pdf_name = kwargs.get("pdf_name", "pdf_contour.pdf")
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt) 
        
        print("INFO: computing contour...")
        
        d = self.m.get_A_data(r, t)
        if style in ["winter", "jet", "rb"]:
            if style == "winter":
                cmap=cntr_winter
            if style == "jet":
                cmap=cntr_jet
            if style == "rb":
                cmap = cntr_rb
            cs = ax.contourf(d.X, d.Y, d.Az,
                             levels=lvls,
                             vmin=np.nanmin(d.Az),
                             vmax=np.nanmax(d.Az),
                             cmap=cmap)
            
            if pdf: # only when printing pdfs for better performance otherwise
                # This is the fix for the white lines between contour levels
                for c in cs.collections:
                    c.set_edgecolor("face")
            
        else: 
            cs = ax.contour(d.X, d.Y, d.Az,
                             levels=lvls,
                             vmin=np.nanmin(d.Az),
                             vmax=np.nanmax(d.Az),
                             colors=style)
            
        #ax.clabel(cs, inline=True, fontsize=self.fgsz)
        # fig.colorbar(cs, shrink = 0.8)
        if pdf:
            plt.savefig(pdf_name, dpi=pdf_dpi, bbox_inches="tight")
        
        print("INFO: finished contour.")
        
    
    def quiver(self, dr, dt, **kwargs):
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        self._set_machine_dims(ax)
        
        scale = kwargs.get("scale", None)
        width = kwargs.get("width", None)
        pdf = kwargs.get("pdf", False)
        pdf_dpi = kwargs.get("pdf_dpi", 300)
        pdf_name = kwargs.get("pdf_name", "pdf_quiver.pdf")
        
        # r = np.linspace(0, self.r_max, dr)
        r = np.linspace(min(self.m.radii), max(self.m.radii), dr)
        t = np.linspace(0, 2*pi, dt)
        
        print("INFO: computing quiver...")
        
        d = self.m.get_B_data(r, t)
        intensity = np.sqrt(d.U**2 + d.V**2)
    
        ax.quiver(d.X, d.Y, d.U, d.V, intensity,
                  cmap=flux, 
                  scale=scale, 
                  width=width) # color="b")
        
        if pdf:
            plt.savefig(pdf_name, dpi=pdf_dpi, bbox_inches="tight")
        
        print("INFO: finished quiver.")
        
        
    def fluxplot(self, dr, dt, **kwargs):
        fig, ax = self._set_up_plot()
        self._set_plot_dims(ax)
        # self._set_machine_dims(ax)
        
        # --- detail ---
        lvls = kwargs.get("lvls", 50)
        lw = kwargs.get("lw", None)
        
        # --- file export ---
        pdf = kwargs.get("pdf", False)
        pdf_dpi = kwargs.get("pdf_dpi", 300)
        pdf_name = kwargs.get("pdf_name", "pdf_fluxplot.pdf")
        svg = kwargs.get("svg", False)
        svg_dpi = kwargs.get("svg_dpi", 300)
        svg_name = kwargs.get("svg_name", "svg_fluxplot.svg")
        
        r = np.linspace(0, self.r_max, dr)
        t = np.linspace(0, 2*pi, dt)
        
        print("INFO: computing fluxplot...")
        
        dA = self.m.get_A_data(r, t)
        dB = self.m.get_B_data(r, t)
        A = np.abs(dA.Az)
        B = np.sqrt(dB.Br**2 + dB.Bt**2)
        # print(np.nanmin(B))
        # print(np.nanmax(B))
        
        cs_cf = ax.contourf(dB.X, dB.Y, B,
                            levels=200,
                            vmin=np.nanmin(B),
                            vmax=np.nanmax(B),
                            cmap=flux)
        
        if pdf or svg: # only when printing pdfs for better performance otherwise
            # This is the fix for the white lines between contour levels
            for c in cs_cf.collections:
                c.set_edgecolor("face")
        
        cs_c = ax.contour(dB.X, dB.Y, A,
                          levels=lvls,
                          vmin=np.nanmin(A),
                          vmax=np.nanmax(A),
                          colors="black",
                          linewidths=lw)
        cbar = fig.colorbar(cs_cf, shrink = 0.8)
        cbar.ax.tick_params(labelsize=self.fgsz)
        
        if pdf:
            plt.savefig(pdf_name, dpi=pdf_dpi, bbox_inches="tight")
        
        if svg:
            plt.savefig(svg_name, dpi=svg_dpi, bbox_inches="tight")

        print("INFO: finished fluxplot.")
    

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
    
      
   
    
    