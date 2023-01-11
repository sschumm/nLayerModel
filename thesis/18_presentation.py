# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi
import tikzplotlib as tkz
import matplotlib.pyplot as plt_mpl
from data import Generator as gn
from modules import Model, AirLayer, MagneticLayer, CurrentLoading
from modules.plot import PlanePlot

n_syn = gn.n_syn


p, l_e = 4, 0.3

r1, r2, r3, r4, r5, r6 = 1.4, 2, 2.3, 2.8, 3.2, 3.6 # [m]
mu_r = 1e5 # [-]

alpha_r = 0.0*pi # [rad] =  0 [°]
alpha_s = 0.5*pi # [rad] = 90 [°]

K_r_hat = 8e5 # [A/m]
K_s_hat = 4e5 # [A/m]


model = Model(p, l_e)

model.add_layer(AirLayer(r=r1))
model.add_layer(MagneticLayer(r=r2, mu_r=mu_r))
model.add_layer(CurrentLoading(K=K_r_hat, r=r3, alpha=alpha_r))
model.add_layer(CurrentLoading(K=K_s_hat, r=r4, alpha=alpha_s))
model.add_layer(AirLayer(r=r5))
model.add_layer(MagneticLayer(r=r6, mu_r=mu_r))

model.build()
model.solve()
# model.tangential_forces()

plt = PlanePlot(model)
plt.fluxplot(dr=500, dt=500, lvls=11)
# plt.just_dims(dr=500, dt=500)