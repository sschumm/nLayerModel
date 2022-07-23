# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 12:05:15 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=200, precision=8)

from Moduls_Tests.model import Model
from Moduls_Tests.layer import CurrentLoading, MagneticLayer
from Moduls_Tests import plotting


# parameters from bumby appendix
p       = 1
K_f0    = 5.31e6
r_f     = 0.425
# random mu_r for iron
r_1     = 0.2
mu_1    = 10e3


m1 = Model(p)
m1.add_layer(MagneticLayer(r=r_1, mu_r=mu_1))
m1.add_layer(CurrentLoading(r=r_f, K=K_f0))
m1.build()

x = m1.solve()

a1, b1, a2, b2, a3, b3 = x 

print("Solution: ")
print(f"a1 - analytic: 6.672743 - numeric: {a1}")
print(f"b1 - analytic:        0 - numeric: {b1}")
print(f"a2 - analytic: 3.336371 - numeric: {a2}")
print(f"b2 - analytic: 0.133455 - numeric: {b2}")
print(f"a3 - analytic:        0 - numeric: {a3}")
print(f"b3 - analytic: 0.736087 - numeric: {b3}")



X1, Y1, U1, V1 = m1.get_B_plot()
X2, Y2, Z2 = m1.get_A_plot()
radii = m1.get_radii_data()


#%%
# =============================================================================
# data = {"streamplot": (X1, Y1, U1, V1, radii, m1.layers),
#         "contour": (X2, Y2, Z2, radii, m1.layers)}
# plotting.multi_figure(2, 1, data=data)
# =============================================================================


#%%
plotting.streamplot_B(X1, Y1, U1, V1, radii, m1.layers)
#%%
plotting.quiver_B(X1, Y1, U1, V1, radii, m1.layers)
#%%
# plotting.contour_A(X2, Y2, Z2, radii, m1.layers)