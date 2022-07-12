# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:36:48 2022

@author: svens
"""
import numpy as np

from Moduls_Tests.layer import CurrentLoading, MagneticLayer, AirLayer
from Moduls_Tests.model import Model


l1 = AirLayer(r = 0.3)                          # r: [m]
l2 = CurrentLoading(r = 0.425, K = 332000)      # K: [A/m]
l3 = AirLayer(r = 1.4)                          # r: [m]
l4 = MagneticLayer(r = 2.0, mu_r = 5000)        # mu_r: 300 - 10.000 [-]



#%%

m1 = Model(p = 2)

m1.add_layer(l1)
m1.add_layer(l2)
m1.add_layer(l3)
m1.add_layer(l4)


#%%

m1.build()


np.set_printoptions(suppress=True, linewidth=200, precision=5)
print("M=\n", m1.M, "\n")
np.set_printoptions(suppress=True, linewidth=200, precision=8)
print("y=\n", m1.y, "\n")

#%%

x = m1.solve()
print("x=\n", x, "\n")
