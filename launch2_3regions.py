# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 08:22:53 2022

@author: svens
"""

import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=3)

from modules.model import Model
from modules.layer import CurrentLoading, MagneticLayer
from modules.plot import PlanePlot

from analytics.three_regions import analytic_solution


# random parameters (scale approx to bumby)
p       = 4
r_1     = 0.7 # m
r_2     = 1.2
K_f     = 5 * 1e6 # A/m
mu_1    = 10 * 1e3

model = Model(p)
model.add_layer(MagneticLayer(r=r_1, mu_r=mu_1))
model.add_layer(CurrentLoading(K=K_f, r=r_2, mu_r=1))
model.build()

x_numeric = model.solve()
x_analytic= analytic_solution(p, r1=r_1, r2=r_2, K=K_f, mu_1=mu_1)

a1_n, b1_n, a2_n, b2_n, a3_n, b3_n = x_numeric
a1_a, b1_a, a2_a, b2_a, a3_a, b3_a = x_analytic

print("Solution: ")
print(f"a1 - analytic: {a1_a}{(24 - len(str(a1_a))) * ' '}- numeric: {a1_n}")
print(f"b1 - analytic: {b1_a}{(24 - len(str(b1_a))) * ' '}- numeric: {b1_n}")
print(f"a2 - analytic: {a2_a}{(24 - len(str(a2_a))) * ' '}- numeric: {a2_n}")
print(f"b2 - analytic: {b2_a}{(24 - len(str(b2_a))) * ' '}- numeric: {b2_n}")
print(f"a3 - analytic: {a3_a}{(24 - len(str(a3_a))) * ' '}- numeric: {a3_n}")
print(f"b3 - analytic: {b3_a}{(24 - len(str(b3_a))) * ' '}- numeric: {b3_n}")
print(f"\nNumeric solution is {np.allclose(x_numeric, x_analytic)}. \n")

p_plot = PlanePlot(model)
# p_plot.streamplot(dr=100, dt=100)
# p_plot.contour(dr=100, dt=100)

#%%

# =============================================================================
# for i in range(1, 1001):
#     p_i = np.random.randint(1, 20)
#     K_i = K_f0 * np.random.rand()
#     r_i = np.random.rand() * np.random.randint(1, 10)
#     
#     model = Model(p_i)
#     model.add_layer(CurrentLoading(K=K_i, r=r_i, mu_r=1))
#     model.build()
# 
#     x_numeric = model.solve()
#     x_analytic= analytic_solution(p_i, r_i, K_i)
#     close = np.allclose(x_numeric, x_analytic)
#     if close:
#         print("=", end="")
#     else: 
#         print("\n ! \n")
#     if i % 100 == 0:
#         print("")
# =============================================================================
    
