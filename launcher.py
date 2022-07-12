# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:36:48 2022

@author: svens
"""
import numpy as np

from Moduls_Tests.layer import Layer, CurrentLoading
from Moduls_Tests.model import Model


l1 = Layer(r = 3)
l2 = Layer(r = 3, mu_r = 2)
l3 = CurrentLoading(K = 7777, r = 4)

#%%

m1 = Model(p = 3)

#m1.add_layer(l1)
#m1.add_layer(l2)
m1.add_layer(l3)


#%%

m1.build()


np.set_printoptions(suppress=True, linewidth=150, precision=8)
print("M=\n", m1.M, "\n")
print("y=\n", m1.y, "\n")

#%%

x = m1.solve()
print("x=\n", x, "\n")
