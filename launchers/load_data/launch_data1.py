# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot, PlaneDoublePlot


# -------- init params --------
from data.load import info, load_json, data

data, layers = data("kostopoulos", verbose=(True))

# -------- create model --------

model = Model(data["p"])
model.add_layer(CurrentLoading(K=layers[0]["K"], r=layers[0]["r"], alpha=pi*0.5))
model.add_layer(CurrentLoading(K=layers[1]["K"], r=layers[1]["r"], alpha=pi*0.0))


model.build()
model.solve()
model.total_torque()

# -------- output -------- 
print("x =", model.x, "\n")
print("M =", model.M, "\n")

rM_plot = RadialMultiPlot(model)
rM_plot.set_Br_details(title="Br")
rM_plot.set_Ht_details(title="Ht")
# rM_plot.multiplot(angle=90)

p_plot = PlanePlot(model)
p_plot.contour(dr=400, dt=200)
# p_plot.streamplot(dr=400, dt=200)

pd_plot = PlaneDoublePlot(model)
# pd_plot.plot_BandA(dr=400, dt=200)




