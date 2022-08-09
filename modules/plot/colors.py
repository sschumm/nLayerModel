# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 07:17:22 2022

@author: svens
"""

import numpy as np
from matplotlib import cm#, patches
from matplotlib.colors import ListedColormap


winter = cm.get_cmap("winter_r")
jet = cm.get_cmap("jet_r")
rb = cm.get_cmap("gist_rainbow")


# build winter colormap for streamplot and quiver
strp_sample = winter(np.linspace(0,1,360))
strp_winter = ListedColormap(strp_sample[104:, :]) # output


#build winter colormap for contour
res, n = 256, 0.65
cntr_sampled = winter(np.linspace(0,1,2*res))
cntr_sampled[:res, -1] = np.flip(np.linspace(0., 1, res))**n
cntr_sampled[res:, -1] = np.linspace(0., 1, res)**n
cntr_winter = ListedColormap(cntr_sampled) # output

# build jet colormap for contour
res, n = 256, 0.5
jet_sampled = jet(np.linspace(0,1,2*res))
jet_sampled[:res, -1] = np.flip(np.linspace(0., 1, res))**n
jet_sampled[res:, -1] = np.linspace(0., 1, res)**n
cntr_jet = ListedColormap(jet_sampled) # output

# build rainbow colormap for contour
front, back, res, n = 20, 172,256, 0.7
rb_sampled = rb(np.linspace(0,1,front + 2*res + back))
rb_sampled = rb_sampled[front:front + 2*res, :]
rb_sampled[:res, -1] = np.flip(np.linspace(0., 1, res))**n
rb_sampled[res:, -1] = np.linspace(0., 1, res)**n
cntr_rb = ListedColormap(rb_sampled) # output



# machine dimensions
border_default = "black"
border_current = "red"
layer_default = "white"
layer_magnetic = "#a6a6a6"