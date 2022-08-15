# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot


# -------- init params --------
from data.load import data
from data.precalculations import kw, K, Kamp

data, layers = data("song_gen1", verbose=(False))
# data, layers = data("song_gen2", verbose=(True))


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

# Kr = K_n_amplitude(T=data["T_p_pole"] * 2*data["p"], 
#                    r=loadings[0]["r"], 
#                    I=data["I_r"] * np.sqrt(2), o=0.5)   
# Kr = 5.447 * 1e6
k_w = kw(1, pi * 0.5)
print(k_w)
K_r = K(m=1, I=680, d=2*1.756, Ns=2*12*1400)
print(K_r)
K_r_amp = Kamp(k_w, K_r)
print(K_r_amp)

K_s = loadings[1]["K"]
print(K_s)
K_s_amp = Kamp(0.9577, K_s)
print(K_s_amp)

model.add_layer(CurrentLoading(K=K_r_amp, r=loadings[0]["r"], 
                               alpha=pi*0.5))
model.add_layer(CurrentLoading(K=K_s_amp, r=loadings[1]["r"], 
                               alpha=pi*0.))

# -------- compute model --------
model.build()
model.solve()
model.total_torque()

# -------- output -------- 
# print("x =", model.x, "\n")
# print("M =", model.M, "[Nm] \n")

print("P_max =", 2*pi* (data["n_rated"]/60) * (model.M / 1e6), "[MW] \n")
print("P_out =", (data["P_out"])/1e6, "[MW] \n")


rM_plot = RadialMultiPlot(model)
rM_plot.set_Br_details(title="Br")
rM_plot.set_Ht_details(title="Ht")
# for angle in np.linspace(0,20,20):  
#     rM_plot.set_Br_details(title=f"Br at {angle}°.")
#     rM_plot.set_Ht_details(title=f"Ht at {angle}°.")
#     rM_plot.multiplot(dr=1000, dt=5000, angle=angle)

p_plot = PlanePlot(model, fgsz=70)
p_plot.contour(dr=400, dt=200, style="rb")
# p_plot.streamplot(dr=400, dt=200)


# d= model.get_B_data(r=np.linspace(0,2,100), 
#                     t=np.linspace(0,2*pi,100))






