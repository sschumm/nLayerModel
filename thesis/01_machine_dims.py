# -*- coding: utf-8 -*-
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as plt_mpl

from modules import Model, AirLayer, MagneticLayer, CurrentLoading
from modules.plot import PlanePlot


r_ri = 0.05
r_ro = 0.1
r_rF = 0.108
r_sA = 0.135
r_si = 0.145
r_so = 0.17

model = Model(p=4, l=0.5)

model.add_layer(AirLayer(r=r_ri))
model.add_layer(MagneticLayer(r=r_ro, mu_r=1e5))
# model.add_layer(CurrentLoading(K=1e4, r=r_rF, alpha=0.0))
# model.add_layer(CurrentLoading(K=1e3, r=r_sA, alpha=0.5))
model.add_layer(AirLayer(r=r_si))
model.add_layer(MagneticLayer(r=r_so, mu_r=1e5))

model.build()
model.solve()
model.total_torque()


#%%

# x1 = np.linspace(-4, 4, 5000)
# y1 = x1**2

# fig = plt_mpl.figure(dpi=1000)
# ax = plt_mpl.subplot()

# ax.plot(x1, y1)


# tkz.clean_figure()
# tkz.save("test1.tex")



#%%
plt_plane = PlanePlot(model, fgsz=50)
# plt_plane.fluxplot(dr=1000, dt=1000, lvls=10)
# plt_plane.contour(dr=1000, dt=1000, pdf=True, pdf_dpi=300)

#%%
# plt_plane.quiver(dr=20, dt=200, scale=300, width=0.001, pdf=True, pdf_dpi=300)

#%%
plt_plane = PlanePlot(model, fgsz=50)
margin = 0.003
# plt_plane.fluxplot(dr=200, dt=200, lvls=10, show_borders=True)
# plt_plane.fluxplot(dr=100, dt=100, lvls=10, lw=None,
#                     r_min=r_ri-margin, r_max=r_so+margin,
#                     t_min=0, t_max=np.pi/2,
#                     custom_dims=True,
#                     y0=0, x0=0,
#                     show_cbar=True, 
#                     show_axis=False,
#                     show_borders=True,
#                     transparent=True,
#                     padding=0.0,
#                     pdf=False, pdf_dpi=300, 
#                     svg=False, svg_dpi=300)

plt_plane.machineplot(r_min=r_ri-margin, r_max=r_so+margin,
                      t_min=0, t_max=np.pi/2,
                      custom_dims=True,
                      y0=0, x0=0,
                      show_cbar=True, 
                      show_axis=False,
                      show_borders=False,
                      transparent=True,
                      border_width=4,
                      padding=0.0,
                      pdf=False, pdf_dpi=300, 
                      svg=True, svg_dpi=300)
