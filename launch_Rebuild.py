# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:32:04 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=3)

from modules.model import Model
from modules.layer import Layer, CurrentLoading
from modules.plot import Plot


# parameters from bumby appendix
p       = 1
K_f0    = 5.31e6
r_f     = 0.425
# random mu_r for iron
r_1     = 0.2
mu_1    = 10e3

m1 = Model(p)

m1.add_layer(Layer(r=r_1, mu_r=mu_1))
m1.add_layer(CurrentLoading(K=K_f0, r=r_f, mu_r=1))

m1.build()
m1.solve()
m1.x

a1, b1, a2, b2, a3, b3 = m1.x 

print("Solution: ")
print(f"a1 - analytic: 6.672743 - numeric: {a1}")
print(f"b1 - analytic:        0 - numeric: {b1}")
print(f"a2 - analytic: 3.336371 - numeric: {a2}")
print(f"b2 - analytic: 0.133455 - numeric: {b2}")
print(f"a3 - analytic:        0 - numeric: {a3}")
print(f"b3 - analytic: 0.736087 - numeric: {b3}")

#%%
p = Plot(m1)
p.streamplot()