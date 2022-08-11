# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot, PlaneDoublePlot


# -------- init params --------
from data.load import data
from data.precalculations import K_ph_amplitude, K_n_amplitude

data, layers = data("song", verbose=(False))
Kr = K_n_amplitude(T=24*1400, r=1.756, I = np.sqrt(2) * 680)
# print(Kr)

# -------- create model --------

model = Model(data["p"], data["l"])

loadings = []
for layer in layers:
    if layer["type"] == "AirLayer":
        model.add_layer(AirLayer(r=layer["r"]))
    elif layer["type"] == "MagneticLayer":
        model.add_layer(MagneticLayer(r=layer["r"], mu_r=layer["mu_r"]))
    else:
        loadings.append(layer)

    
# model.add_layer(CurrentLoading(K=loadings[0]["K"], r=loadings[0]["r"], alpha=pi*0.0))
model.add_layer(CurrentLoading(K=Kr, r=loadings[0]["r"], alpha=pi*0.5))
model.add_layer(CurrentLoading(K=loadings[1]["K"], r=loadings[1]["r"], alpha=pi*0.))
# model.add_layer(CurrentLoading(K=Ks, r=loadings[1]["r"], alpha=pi*0.5))



model.build()
model.solve()
model.total_torque()

# -------- output -------- 
# print("x =", model.x, "\n")
print("M =", model.M, "[Nm] \n")

print("P_max =", 2*pi* (data["n_rated"]/60) * (model.M / 1e6), "[MW] \n")
print("P_out =", (data["P_out"])/1e6, "[MW] \n")


# rM_plot = RadialMultiPlot(model)
# rM_plot.set_Br_details(title="Br")
# rM_plot.set_Ht_details(title="Ht")
# rM_plot.multiplot(angle=90)

# p_plot = PlanePlot(model, fgsz=100)
# p_plot.contour(dr=400, dt=200)
# p_plot.streamplot(dr=400, dt=200)

# pd_plot = PlaneDoublePlot(model)
# pd_plot.plot_BandA(dr=400, dt=200)




