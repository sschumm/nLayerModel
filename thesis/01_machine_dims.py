# -*- coding: utf-8 -*-
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as plt_mpl

from modules import Model, AirLayer, MagneticLayer, CurrentLoading
from modules.plot import PlanePlot


model = Model(p=4, l=0.5)

model.add_layer(AirLayer(r=0.4))
model.add_layer(MagneticLayer(r=0.5, mu_r=1e5))
model.add_layer(CurrentLoading(K=1e4, r=0.6, alpha=0.0))
model.add_layer(CurrentLoading(K=1e3, r=0.8, alpha=0.5))
model.add_layer(AirLayer(r=0.9))
model.add_layer(MagneticLayer(r=1.0, mu_r=1e5))

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
# plt_plane.quiver(dr=20, dt=200, scale=250, width=0.001, pdf=True, pdf_dpi=300)

#%%
plt_plane = PlanePlot(model, fgsz=20)
# plt_plane.fluxplot(dr=200, dt=200, lvls=10, pdf=True, pdf_dpi=300)
plt_plane.fluxplot(dr=150, dt=150, lvls=10, svg=True, svg_dpi=300)