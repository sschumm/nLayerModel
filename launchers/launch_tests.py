# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 11:08:20 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=3)

from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.plane import PlanePlot, PlaneDoublePlot
from modules.plot.radial import RadialPlot, RadialMsizePlot

p = 2


model = Model(p)

# model.add_layer(CurrentLoading(K=5e6, r=0.4, mu_r=1e5))

model.add_layer(AirLayer(r=0.2))
model.add_layer(MagneticLayer(r=0.3, mu_r=8))
model.add_layer(AirLayer(r=0.4))

model.add_layer(CurrentLoading(K=5e6, r=0.5, mu_r=1, alpha=(np.pi/2) * 0.))

model.add_layer(MagneticLayer(r=0.6, mu_r=7))
model.add_layer(AirLayer(r=0.7))

model.add_layer(MagneticLayer(r=0.8, mu_r=5))
# model.add_layer(CurrentLoading(K=3e6, r=0.85))


model.build()
model.solve()

M = model.total_torque()
# print(M)


p_plot = PlanePlot(model)
# pD_plot = PlaneDoublePlot(model)
# r_plot = RadialPlot(model)
rM_plot = RadialMsizePlot(model)

detail = 400
# p_plot.contour(detail, detail)
# p_plot.streamplot(detail, detail)
# p_plot.quiver(40, 40)
angle=80
rM_plot.plot_radial_Az(angle=angle)
rM_plot.plot_radial_Br(angle=angle)
rM_plot.plot_radial_Ht(angle=angle)


