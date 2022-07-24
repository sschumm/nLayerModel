# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 11:08:20 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from modules.model import Model
from modules.layer import Layer, AirLayer, CurrentLoading
from scipy.constants import pi

p = 2

cl1 = CurrentLoading(K=5e6, r=0.5, mu_r=1)
al1 = AirLayer(r=0.7)
ml1 = Layer(r=0.8, mu_r=10e3)

model = Model(p)

model.add_layer(cl1)
model.add_layer(al1)
model.add_layer(ml1)

model.build()
model.solve()

print("x= \n", model.x, "\n")
print("sysA= \n", model.sysA, "\n")
print("sysb= \n", model.sysb, "\n")

print("")

print("A= \n", cl1.A(p, R=0.2, a_j=model.x[0], b_j=model.x[1]), "\n")
print("dA= \n", cl1.dA(p, R=0.2, a_j=model.x[0], b_j=model.x[1]), "\n")
print("Az= \n", cl1.Az(p, R=0.2, T=pi, a_j=model.x[0], b_j=model.x[1]), "\n")
print("Br= \n", cl1.Br(p, R=0.2, T=pi, a_j=model.x[0], b_j=model.x[1]), "\n")
print("Bt= \n", cl1.Bt(p, R=0.2, T=pi, a_j=model.x[0], b_j=model.x[1]), "\n")