# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=3)

from modules.model import Model
from modules.layer import CurrentLoading, MagneticLayer
from modules.plot.plane import PlanePlot

from analytics.four_regions import analytic_solution_K2


# random parameters (scale approx to bumby)
p       = 2
r_1     = 0.6 # m
r_2     = 0.8
r_3     = 1.0
K_f     = 5 * 1e6 # A/m
mu_r1   = 1 * 1e4
mu_r3   = 1 * 1e3

model = Model(p)
model.add_layer(CurrentLoading(K=0, r=r_1, mu_r=mu_r1))
model.add_layer(CurrentLoading(K=K_f, r=r_2, mu_r=1))
model.add_layer(MagneticLayer(r=r_3, mu_r=mu_r3))
model.build()

x_numeric = model.solve()
x_analytic= analytic_solution_K2(p, K2=K_f, r1=r_1, r2=r_2, r3=r_3,  
                              mu_r1=mu_r1, mu_r3=mu_r3)

a1_n, b1_n, a2_n, b2_n, a3_n, b3_n, a4_n, b4_n = x_numeric[1]
a1_a, b1_a, a2_a, b2_a, a3_a, b3_a, a4_a, b4_a = x_analytic

print("Solution: ")
print(f"a1 - analytic: {a1_a}{(24 - len(str(a1_a))) * ' '}- numeric: {a1_n}")
print(f"b1 - analytic: {b1_a}{(24 - len(str(b1_a))) * ' '}- numeric: {b1_n}")
print(f"a2 - analytic: {a2_a}{(24 - len(str(a2_a))) * ' '}- numeric: {a2_n}")
print(f"b2 - analytic: {b2_a}{(24 - len(str(b2_a))) * ' '}- numeric: {b2_n}")
print(f"a3 - analytic: {a3_a}{(24 - len(str(a3_a))) * ' '}- numeric: {a3_n}")
print(f"b3 - analytic: {b3_a}{(24 - len(str(b3_a))) * ' '}- numeric: {b3_n}")
print(f"a4 - analytic: {a4_a}{(24 - len(str(a4_a))) * ' '}- numeric: {a4_n}")
print(f"b4 - analytic: {b4_a}{(24 - len(str(b4_a))) * ' '}- numeric: {b4_n}")

atol = 1e-8
close_default = np.allclose(x_numeric[1], x_analytic)
close_edit = np.allclose(x_numeric[1], x_analytic, atol=atol, rtol=0)
print(f"\nNumeric solution is {close_default}. ")
print(f"\nNumeric solution is {close_edit} with a tolerance of {atol}. \n")

p_plot = PlanePlot(model)
# p_plot.streamplot(dr=1000, dt=200)
# p_plot.quiver(50, 40)
# p_plot.contour(dr=100, dt=100)

#%%

# =============================================================================
# for i in range(1, 1001):
#     p_i = np.random.randint(1, 8)
#     K_i = K_f * np.random.rand()
#     r_i1 = np.random.rand() * np.random.randint(1, 10) + 0.1
#     r_i2 = r_i1 * np.random.randint(2, 10)
#     r_i3 = r_i2 * np.random.randint(2, 10)
#     mu_ir1 = 10000 + 1000 * np.random.randint(-5, 5)
#     mu_ir3 = 10000 + 1000 * np.random.randint(-5, 5)
#     
#     model = Model(p_i)
#     model.add_layer(CurrentLoading(K=0, r=r_i1, mu_r=mu_ir1))
#     model.add_layer(CurrentLoading(K=K_i, r=r_i2, mu_r=1))
#     model.add_layer(MagneticLayer(r=r_i3, mu_r=mu_ir3))
#     model.build()
# 
#     x_numeric = model.solve()
#     x_analytic= analytic_solution_K2(p_i, K2=K_i, r1=r_i1, r2=r_i2, r3=r_i3,  
#                                   mu_r1=mu_ir1, mu_r3=mu_ir3)
#     
#     atol = 1e-5
#     close_default = np.allclose(x_numeric[0], x_analytic)
#     close_edit = np.allclose(x_numeric[0], x_analytic, atol=atol, rtol=0)
#     
#     if close_default: #  and close_edit:
#         print("=", end="")
#     else: 
#         print("\n ! \n")
#         a1_n, b1_n, a2_n, b2_n, a3_n, b3_n, a4_n, b4_n = x_numeric[0]
#         a1_a, b1_a, a2_a, b2_a, a3_a, b3_a, a4_a, b4_a = x_analytic
#         print("Solution: ")
#         print(f"a1 - analytic: {a1_a}{(24 - len(str(a1_a))) * ' '}- numeric: {a1_n}")
#         print(f"b1 - analytic: {b1_a}{(24 - len(str(b1_a))) * ' '}- numeric: {b1_n}")
#         print(f"a2 - analytic: {a2_a}{(24 - len(str(a2_a))) * ' '}- numeric: {a2_n}")
#         print(f"b2 - analytic: {b2_a}{(24 - len(str(b2_a))) * ' '}- numeric: {b2_n}")
#         print(f"a3 - analytic: {a3_a}{(24 - len(str(a3_a))) * ' '}- numeric: {a3_n}")
#         print(f"b3 - analytic: {b3_a}{(24 - len(str(b3_a))) * ' '}- numeric: {b3_n}")
#         print(f"a4 - analytic: {a4_a}{(24 - len(str(a4_a))) * ' '}- numeric: {a4_n}")
#         print(f"b4 - analytic: {b4_a}{(24 - len(str(b4_a))) * ' '}- numeric: {b4_n}")
#         
#         print(p_i)
#         print(K_i)
#         print(r_i1)
#         print(r_i2)
#         break
#     if i % 100 == 0:
#         print("")
# =============================================================================
    
