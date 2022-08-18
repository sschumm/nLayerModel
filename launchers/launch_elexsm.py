# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:42:44 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot, PlaneDoublePlot


# -------- init params --------
p = 3
wd = 0.1

r_ri = 0.4
r_ro = 0.7
m_rr = 8 * 1e3
K_r  = 4 * 1e5

r_si = 1.4
r_so = 2.0
m_rs = 5 * 1e3
K_s  = 2 * 1e5


# -------- create model --------
model = Model(p)

model.add_layer(AirLayer(r=r_ri))
model.add_layer(MagneticLayer(r=r_ro, mu_r=m_rr))
model.add_layer(CurrentLoading(K=K_r, r=r_ro+wd, alpha=pi * 0.0))
model.add_layer(CurrentLoading(K=K_s, r=r_si-wd, alpha=pi * 0.5))
model.add_layer(AirLayer(r_si))
model.add_layer(MagneticLayer(r_so, mu_r=m_rs))

model.build()
model.solve()
model.total_torque()

# -------- output -------- 
print("x =", model.x, "\n")
print("M =", model.M, "\n")

rM_plot = RadialMultiPlot(model)
rM_plot.set_Az_details(title="Az")
rM_plot.set_Br_details(title="Br")
rM_plot.set_Ht_details(title="Ht")
# rM_plot.multiplot(["Az"], angle = 45)
# rM_plot.multiplot(["Az", "Br"])

p_plot = PlanePlot(model)
# p_plot.contour(dr=400, dt=200, style="black")
p_plot.contour(dr=400, dt=200, style="jet")
p_plot.fluxplot(dr=400, dt=200)
# p_plot.streamplot(dr=400, dt=200)

pd_plot = PlaneDoublePlot(model)
# pd_plot.plot_BandA(dr=400, dt=200)



