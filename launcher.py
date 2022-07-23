# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 15:36:48 2022

@author: svens
"""
import numpy as np
from Moduls_Tests import plotting

from Moduls_Tests.layer import CurrentLoading, MagneticLayer, AirLayer
from Moduls_Tests.model import Model

m1 = Model(p = 4)


# =============================================================================
# l1 = AirLayer(r = 2)                          # r: [m]
# l2 = CurrentLoading(r = 0.8, K = 332000)      # K: [A/m]
# l3 = AirLayer(r = 1.4)                          # r: [m]
# l4 = MagneticLayer(r = 1.4, mu_r = 10000)        # mu_r: 300 - 10.000 [-]
# 
# m1.add_layer(l2)
# m1.add_layer(l3)
# # m1.add_layer(l4)
# m1.add_layer(l1)
# =============================================================================


m1.add_layer(AirLayer(r = 0.4))
m1.add_layer(MagneticLayer(r = 0.7, mu_r=10000))
m1.add_layer(CurrentLoading(r= 0.8, K = 300000))
m1.add_layer(AirLayer(r = 1))
m1.add_layer(CurrentLoading(r= 1.3, K = 200000))
m1.add_layer(AirLayer(r = 1.4))
m1.add_layer(MagneticLayer(r = 2.0, mu_r=10000))
# m1.add_layer(AirLayer(r = 2.5))


#%%
m1.build()

# =============================================================================
# np.set_printoptions(suppress=True, linewidth=200, precision=5)
# print("M=\n", m1.M, "\n")
# np.set_printoptions(suppress=True, linewidth=200, precision=8)
# print("y=\n", m1.y, "\n")
# =============================================================================

#%%
x = m1.solve()

# print("x=\n", x, "\n")
# print("M=\n", m1.M, "\n")

#%%
X1, Y1, U1, V1 = m1.get_B_plot()
X2, Y2, Z2 = m1.get_A_plot()
radii = m1.get_radii_data()

# =============================================================================
# np.set_printoptions(suppress=True, linewidth=200, precision=2)
# print("R=\n", R, "\n")
# print("THETA=\n", THETA, "\n")
# print("Br=\n", Br, "\n")
# print("B0=\n", B0, "\n")
# print("X=\n", X, "\n")
# print("Y=\n", Y, "\n")
# print("U=\n", U, "\n")
# print("V=\n", V, "\n")
# =============================================================================

#%%
plotting.streamplot_B(X1, Y1, U1, V1, radii, m1.layers)
#%%
# plotting.quiver_B(X1, Y1, U1, V1, radii, m1.layers)
#%%
# plotting.contour_A(X2, Y2, Z2, radii, m1.layers)
#%%
# =============================================================================
# data = {"streamplot": (X1, Y1, U1, V1, radii, m1.layers),
#         "contour": (X2, Y2, Z2, radii, m1.layers)}
# plotting.multi_figure(2, 1, data=data)
# =============================================================================


































