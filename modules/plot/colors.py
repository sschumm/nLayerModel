# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 07:17:22 2022

@author: svens
"""

import numpy as np
from matplotlib import cm#, patches
from matplotlib.colors import ListedColormap


winter = cm.get_cmap("winter_r")
strp_sample = winter(np.linspace(0,1,360))
strp_winter = ListedColormap(strp_sample[104:, :])

res, n = 256, 0.65
cntr_sampled = winter(np.linspace(0,1,2*res))
cntr_sampled[:res, -1] = np.flip(np.linspace(0., 1, res))**n
cntr_sampled[res:, -1] = np.linspace(0., 1, res)**n
cntr_winter = ListedColormap(cntr_sampled)