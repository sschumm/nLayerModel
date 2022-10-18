# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi
import tikzplotlib as tkz
import matplotlib.pyplot as plt_mpl

from modules import Model, AirLayer, MagneticLayer, CurrentLoading
from modules.plot import PlanePlot



p, l_e = 4, 0.3

r1, r2, r3, r4, r5, r6 = 1.4, 2, 2.3, 2.8, 3.2, 3.6 # [m]
mu_r = 1e5 # [-]

alpha_r = 0.0*pi # [rad] =  0 [°]
alpha_s = 0.5*pi # [rad] = 90 [°]

K_r_hat = 1e3 # [A/m]
K_s_hat = 4e2 # [A/m]


model = Model(p, l_e)

model.add_layer(AirLayer(r=r1))
model.add_layer(MagneticLayer(r=r2, mu_r=mu_r))
model.add_layer(CurrentLoading(K=K_r_hat, r=r3, alpha=alpha_r))
model.add_layer(CurrentLoading(K=K_s_hat, r=r4, alpha=alpha_s))
model.add_layer(AirLayer(r=r5))
model.add_layer(MagneticLayer(r=r6, mu_r=mu_r))

model.build()
x = model.solve()
# model.total_torque()
#%%
# plt_plane = PlanePlot(model, fgsz=50)
# plt_plane.fluxplot(dr=1000, dt=1000, lvls=10)

plt_plane = PlanePlot(model, fgsz=50)
margin = 0.05
# plt_plane.fluxplot(dr=200, dt=200, lvls=10, show_borders=True)
plt_plane.fluxplot(dr=100, dt=100, lvls=10, lw=None,
                    r_min=1.4-margin, r_max=3.6+margin,
                    t_min=0, t_max=np.pi/2,
                    custom_dims=True,
                    y0=0, x0=0,
                    show_cbar=True, 
                    show_axis=False,
                    show_borders=True,
                    transparent=True,
                    padding=0.1,
                    pdf=False, pdf_dpi=300, 
                    svg=False, svg_dpi=300)