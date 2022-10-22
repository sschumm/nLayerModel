# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi
import tikzplotlib as tkz
import matplotlib.pyplot as plt_mpl
from data import Generator as gn
from modules import Model, AirLayer, MagneticLayer, CurrentLoading

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



x1, x2 = model.solve()
x1.shape  # (14,)
x2.shape  # (14,)

A = model.get_A_data(r=np.linspace(r1, r6, 1000), t=np.linspace(0, 2*pi, 1000))
A.Az.shape # (1000, 1000)

B = model.get_B_data(r=np.array([r3]), t=np.linspace(0, pi, 1000))
B.Br.shape # (1000, 1)

H = model.get_H_data(r=np.linspace(r1, r6, 1000), t=np.array([0.25*pi]))
H.Ht.shape # (1, 1000)

M = model.total_torque()
P = 2*pi*n_syn*M

#%% angle sweep for torque demo

def angle_sweep(alpha_s):
    model = Model(p, l_e)

    model.add_layer(AirLayer(r=r1))
    model.add_layer(MagneticLayer(r=r2, mu_r=mu_r))
    model.add_layer(CurrentLoading(K=K_r_hat, r=r3, alpha=alpha_r))
    model.add_layer(CurrentLoading(K=K_s_hat, r=r4, alpha=alpha_s))
    model.add_layer(AirLayer(r=r5))
    model.add_layer(MagneticLayer(r=r6, mu_r=mu_r))

    model.build()
    model.solve()
    model.total_torque()
    return model.Mpos


x = np.linspace(0, pi, 1000)
y = []
for alpha_s in x:
    y.append(angle_sweep(alpha_s))

# import matplotlib.pyplot as plt
# plt.figure(dpi=300)
# plt.plot(x, y)

# tkz.clean_figure()
# tkz.save("msm_torque_over_angle.tex")


#%%
from modules.plot import PlanePlot
# PlanePlot(model).fluxplot(dr=1000, dt=1000, lvls=10)
# PlanePlot(model).quiver(dr=20, dt=50)
# PlanePlot(model).streamplot(dr=20, dt=20)
# PlanePlot(model).contour(dr=100, dt=100)

from modules.plot import RadialPlot
# RadialPlot(model).plot_radial_Ht(t=0*pi, tikz=True)
# RadialPlot(model).plot_radial_Br(t=0.5*pi)





# model.total_torque()
#%%
# plt_plane = PlanePlot(model, fgsz=50)
# plt_plane.fluxplot(dr=1000, dt=1000, lvls=10)

# plt_plane = PlanePlot(model, fgsz=50)
# margin = 0.05
# # plt_plane.fluxplot(dr=200, dt=200, lvls=10, show_borders=True)
# plt_plane.fluxplot(dr=70, dt=80, lvls=10, lw=None,
#                     r_min=1.4-margin, r_max=3.6+margin,
#                     t_min=0, t_max=np.pi/2,
#                     custom_dims=True,
#                     y0=0, x0=0,
#                     show_cbar=True, 
#                     show_axis=False,
#                     show_borders=True,
#                     transparent=True,
#                     padding=0,
#                     pdf=False, pdf_dpi=300, 
#                     svg=True, svg_dpi=300)