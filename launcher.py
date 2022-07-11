# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:36:48 2022

@author: svens
"""

from Moduls_Tests.layer import Layer
from Moduls_Tests.model import Model


l1 = Layer(r = 4)
l2 = Layer(r = 3, mu_r = 2)

m1 = Model()

m1.add_layer(l1)
m1.add_layer(l2)

