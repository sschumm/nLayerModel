# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 11:08:20 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot, PlaneDoublePlot, RadialPlot, RadialMsizePlot

p = 3

cl1 = CurrentLoading(K=5e6, r=0.5, mu_r=1)
al1 = AirLayer(r=0.7)
ml1 = MagneticLayer(r=0.8, mu_r=1e5)

model = Model(p)

model.add_layer(cl1)
model.add_layer(al1)
model.add_layer(ml1)

model.build()
model.solve()

# =============================================================================
# print("x= \n", model.x, "\n")
# print("sysA= \n", model.sysA, "\n")
# print("sysb= \n", model.sysb, "\n")
# 
# print("")
# 
# i = 1
# r = np.array([4., 3, 2, 1])
# 
# print("A= \n", model.layers[i].A(p, R=r, a_j=model.x[i], b_j=model.x[i+1]), "\n")
# print("dA= \n", model.layers[i].dA(p, R=r, a_j=model.x[i], b_j=model.x[i+1]), "\n")
# print("Az= \n", model.layers[i].Az(p, R=r, T=pi, a_j=model.x[i], b_j=model.x[i+1]), "\n")
# print("Br= \n", model.layers[i].Br(p, R=r, T=pi, a_j=model.x[i], b_j=model.x[i+1]), "\n")
# print("Bt= \n", model.layers[i].Bt(p, R=r, T=pi, a_j=model.x[i], b_j=model.x[i+1]), "\n")
# =============================================================================

#%%
p_plot = PlanePlot(model)
# p_plot.streamplot(400, 400)

#%%
p_plot.contour(400, 400)

#%%
pD_plot = PlaneDoublePlot(model)
# pD_plot.plot_BandA(100, 100)

#%%
# r_plot = RadialPlot(model)
# r_plot.plot_radial_Az(angle=0.25)
# r_plot.plot_radial_Br(angle=0.25)
# r_plot.plot_radial_Ht(angle=0.25)

#%%
rM_plot = RadialMsizePlot(model)
angle=0.
rM_plot.plot_radial_Az(angle=angle)
rM_plot.plot_radial_Br(angle=angle)
# rM_plot.plot_radial_Ht(angle=0.25)


