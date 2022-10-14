# -*- coding: utf-8 -*-
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


plt_plane = PlanePlot(model, fgsz=100)
plt_plane.fluxplot(dr=1000, dt=1000, lvls=10)