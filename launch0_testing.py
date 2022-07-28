# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 11:08:20 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=3)

from modules.model2 import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.plane import PlanePlot, PlaneDoublePlot
from modules.plot.radial import RadialPlot, RadialMsizePlot

p = 4

cl1 = CurrentLoading(K=5e6, r=0.5, mu_r=1)
al1 = AirLayer(r=0.7)
ml1 = MagneticLayer(r=0.8, mu_r=1e5)

model = Model(p)

model.add_layer(cl1)
model.add_layer(al1)
model.add_layer(ml1)

# =============================================================================
# model.add_layer(CurrentLoading(K=3e6, r=0.85))
# model.add_layer(AirLayer(r=0.2))
# model.add_layer(CurrentLoading(K=2e6, r=0.4, mu_r=1e5))
# model.add_layer(MagneticLayer(r=0.45, mu_r=1e5))
# =============================================================================

model.build()
model.solve()

# =============================================================================
# p_plot = PlanePlot(model)
# pD_plot = PlaneDoublePlot(model)
# r_plot = RadialPlot(model)
# rM_plot = RadialMsizePlot(model)
# =============================================================================

# =============================================================================
# p_plot.contour(100, 100)
# angle=0.17
# rM_plot.plot_radial_Az(angle=angle)
# rM_plot.plot_radial_Br(angle=angle)
# rM_plot.plot_radial_Ht(angle=angle)
# =============================================================================


