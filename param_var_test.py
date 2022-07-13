# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 07:24:32 2022

@author: svens
"""

import numpy as np

from Moduls_Tests.layer import CurrentLoading, MagneticLayer, AirLayer
from Moduls_Tests.model import Model


l1 = AirLayer(r = 0.3)
l2 = MagneticLayer(r = 0.5, mu_r = 5000)
l3 = CurrentLoading(r = 0.55, K = 300000)
l4 = AirLayer(r = 0.65)
l5 = CurrentLoading(r = 0.7, K = 100000)
l6 = AirLayer(r = 0.7)
l7 = MagneticLayer(r = 1.0, mu_r = 5000)

#%%
m1 = Model(p = 2)

m1.add_layer(l1)
m1.add_layer(l2)
m1.add_layer(l3)
m1.add_layer(l4)
m1.add_layer(l5)
m1.add_layer(l6)
m1.add_layer(l7)

m1.build()

np.set_printoptions(suppress=True, linewidth=250, precision=0)
print("M=\n", m1.M, "\n")
np.set_printoptions(suppress=True, linewidth=200, precision=8)
print("y=\n", m1.y, "\n")

x = m1.solve()
print("x=\n", x, "\n")

#%%
for i in range(0, 101):
    print("Progress: ", i, "/ 100")
    m1 = Model(p = 2)
    for j in range(0, 100):
        m1.add_layer(l1)
        m1.add_layer(l2)
        m1.add_layer(l3)
        m1.add_layer(l4)
        m1.add_layer(CurrentLoading(r = 0.7, K = 100000 + i * 10))
        m1.add_layer(l6)
        m1.add_layer(MagneticLayer(r = 1.0, mu_r = 5000 + i * 30))
        
        m1.build()
        x = m1.solve()
        

for j in range(0, 100):
    for i in range(0, 100):
        m1.vary_param(l2, "r", 0.4 + i * 0.001)
















