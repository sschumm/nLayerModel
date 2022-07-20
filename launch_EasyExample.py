# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 15:26:23 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=200, precision=8)

from Moduls_Tests.model import Model
from Moduls_Tests.layer import CurrentLoading
from Moduls_Tests import plotting


# parameters from bumby appendix
p       = 1
K_f0    = 5.31e6
r_f     = 0.425


m1 = Model(p)
m1.add_layer(CurrentLoading(r=r_f, K=K_f0))
m1.build()

x = m1.solve()

a1, b1, a2, b2 = x 

print("Solution: ")
print(f"a1 - analytic: 3.3364 - numeric: {a1}")
print(f"b1 - analytic:      0 - numeric: {b1}")
print(f"a2 - analytic:      0 - numeric: {a2}")
print(f"b2 - analytic: 0.6026 - numeric: {b2}")



X1, Y1, U1, V1 = m1.get_B_plot()
X2, Y2, Z2 = m1.get_A_plot()
radii = m1.get_radii_data()


#%%
data = {"streamplot": (X1, Y1, U1, V1, radii, m1.layers),
        "contour": (X2, Y2, Z2, radii, m1.layers)}
plotting.multi_figure(2, 1, data=data)


#%%
plotting.streamplot_B(X1, Y1, U1, V1, radii, m1.layers)
#%%
plotting.quiver_B(X1, Y1, U1, V1, radii, m1.layers)
#%%
plotting.contour_A(X2, Y2, Z2, radii, m1.layers)