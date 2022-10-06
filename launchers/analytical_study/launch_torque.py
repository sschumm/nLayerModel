# -*- coding: utf-8 -*-
import math
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=6)

from modules.model import Model
from modules.layer import CurrentLoading, MagneticLayer
from modules.plot.plane import PlanePlot
from modules.plot.radial import RadialMultiPlot

from analytics.four_regions import analytic_solution_K1, analytic_solution_K2
from analytics.four_regions import analytic_torque_on_K1


# random parameters (scale approx to bumby)
p       = 2
l       = 1.0
r_1     = 0.6 # m
r_2     = 0.8
r_3     = 1.0
K_1     = 5 * 1e6 # A/m
K_2     = 3 * 1e6
mu_r1   = 1 * 1e4
mu_r3   = 1 * 1e3

alpha1 = (np.pi/2) * 0.0
alpha2 = (np.pi/2) * 0.5

model = Model(p, l)
model.add_layer(CurrentLoading(K=K_1, r=r_1, mu_r=mu_r1, alpha=alpha1))
model.add_layer(CurrentLoading(K=K_2, r=r_2, mu_r=1, alpha=alpha2))
model.add_layer(MagneticLayer(r=r_3, mu_r=mu_r3))
model.build()

x_numeric = model.solve()
x_analytic= [
    analytic_solution_K1(p, K1=K_1, r1=r_1, r2=r_2, r3=r_3,  
                                  mu_r1=mu_r1, mu_r3=mu_r3),
    analytic_solution_K2(p, K2=K_2, r1=r_1, r2=r_2, r3=r_3,  
                              mu_r1=mu_r1, mu_r3=mu_r3)    
    ]


close_default = np.allclose(x_numeric, x_analytic)
print(f"\nNumeric solution is {close_default}. ")

M_numeric = model.total_torque()
M_analytic = analytic_torque_on_K1(p, l, K1=K_1, r1=r_1, 
                                   a1_K2=x_analytic[1][0],
                                   alpha1=alpha1, alpha2=alpha2)

print("M_numeric = ", M_numeric)
print("M_analytic =", M_analytic)


p_plot = PlanePlot(model)
# p_plot.streamplot(dr=1000, dt=200)
# p_plot.quiver(50, 40)
p_plot.contour(dr=100, dt=100)
p_plot.contour(dr=100, dt=100, style="jet")
p_plot.contour(dr=100, dt=100, style="rb")
p_plot.contour(dr=400, dt=200, style="black")

rM_plot = RadialMultiPlot(model)
rM_plot.set_Az_details(title="Az")
rM_plot.set_Br_details(title="Br")
rM_plot.set_Ht_details(title="Ht")
rM_plot.multiplot(["Az", "Br", "Ht"], angle = 45)


#%% angle sweep
if True:
    print("INFO: computing load angle sweep...")
    x_angle = (np.pi/2) * np.linspace(0, 2, 20)
    y_torque = list()
    alpha2 = (np.pi/2) * 0.0

    for i in x_angle:
        alpha1 = i
    
        model = Model(p, l)
        model.add_layer(CurrentLoading(K=K_1, r=r_1, mu_r=mu_r1, alpha=alpha1))
        model.add_layer(CurrentLoading(K=K_2, r=r_2, mu_r=1, alpha=alpha2))
        model.add_layer(MagneticLayer(r=r_3, mu_r=mu_r3))
        model.build()
        model.solve()    
    
        M_numeric = model.total_torque()
        y_torque.append(M_numeric)
        # print("M_numeric = ", M_numeric)
        
        if False:
            p_plot = PlanePlot(model)
            p_plot.contour(dr=100, dt=100, style="jet")
    
    # plot torque over angle
    if True:
        y_torque = np.array(y_torque)
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 7))
        plt.plot(x_angle, y_torque)
    print("INFO: finished load angle sweep.")

#%% Analytic vs. Numeric
if False:
    for i in range(1, 1001):
        p_i = np.random.randint(1, 30)
        K_i1 = K_1 * np.random.rand()
        K_i2 = K_2 * np.random.rand()
        r_i1 = np.random.rand() * np.random.randint(1, 10) + 0.1
        r_i2 = r_i1 * np.random.randint(2, 10)
        r_i3 = r_i2 * np.random.randint(2, 10)
        mu_ir1 = 10000 + 1000 * np.random.randint(-5, 5)
        mu_ir3 = 10000 + 1000 * np.random.randint(-5, 5)
        alphai1 = (np.pi/2) * np.random.randint(0, 50)*0.01
        alphai2 = (np.pi/2) * np.random.randint(0, 50)*0.01
        
        model = Model(p_i)
        model.add_layer(CurrentLoading(K=K_i1, r=r_i1, mu_r=mu_ir1, alpha=alphai1))
        model.add_layer(CurrentLoading(K=K_i2, r=r_i2, mu_r=1, alpha=alphai2))
        model.add_layer(MagneticLayer(r=r_i3, mu_r=mu_ir3))
        model.build()
    
        x_numeric = model.solve()
        x_analytic= [
        analytic_solution_K1(p_i, K1=K_i1, r1=r_i1, r2=r_i2, r3=r_i3,  
                                      mu_r1=mu_ir1, mu_r3=mu_ir3),
        analytic_solution_K2(p_i, K2=K_i2, r1=r_i1, r2=r_i2, r3=r_i3,  
                                  mu_r1=mu_ir1, mu_r3=mu_ir3)    
            ]
        
        close_default = np.allclose(x_numeric, x_analytic)
        
        M_analytic = analytic_torque_on_K1(p_i, l, K1=K_i1, r1=r_i1, 
                                           a1_K2=x_analytic[1][0],
                                           alpha1=alphai1, alpha2=alphai2)
        M_numeric = model.total_torque()
        
        close_torque_np = np.allclose(M_numeric, M_analytic)
        close_torque_math = math.isclose(M_numeric, M_analytic, abs_tol=1e-6)
        
        if close_default and close_torque_math: #close_torque_np:
            print("=", end="")
        else: 
            print("\n ! \n")
            print(f"Solution is {close_default}")
            a1_n, b1_n, a2_n, b2_n, a3_n, b3_n, a4_n, b4_n = x_numeric[0]
            a1_a, b1_a, a2_a, b2_a, a3_a, b3_a, a4_a, b4_a = x_analytic[0]
            print("Solution: ")
            print(f"a1 - analytic: {a1_a}{(24 - len(str(a1_a))) * ' '}- numeric: {a1_n}")
            print(f"b1 - analytic: {b1_a}{(24 - len(str(b1_a))) * ' '}- numeric: {b1_n}")
            print(f"a2 - analytic: {a2_a}{(24 - len(str(a2_a))) * ' '}- numeric: {a2_n}")
            print(f"b2 - analytic: {b2_a}{(24 - len(str(b2_a))) * ' '}- numeric: {b2_n}")
            print(f"a3 - analytic: {a3_a}{(24 - len(str(a3_a))) * ' '}- numeric: {a3_n}")
            print(f"b3 - analytic: {b3_a}{(24 - len(str(b3_a))) * ' '}- numeric: {b3_n}")
            print(f"a4 - analytic: {a4_a}{(24 - len(str(a4_a))) * ' '}- numeric: {a4_n}")
            print(f"b4 - analytic: {b4_a}{(24 - len(str(b4_a))) * ' '}- numeric: {b4_n}")
            
            print("")
            
            a1_n, b1_n, a2_n, b2_n, a3_n, b3_n, a4_n, b4_n = x_numeric[1]
            a1_a, b1_a, a2_a, b2_a, a3_a, b3_a, a4_a, b4_a = x_analytic[1]
            print("Solution: ")
            print(f"a1 - analytic: {a1_a}{(24 - len(str(a1_a))) * ' '}- numeric: {a1_n}")
            print(f"b1 - analytic: {b1_a}{(24 - len(str(b1_a))) * ' '}- numeric: {b1_n}")
            print(f"a2 - analytic: {a2_a}{(24 - len(str(a2_a))) * ' '}- numeric: {a2_n}")
            print(f"b2 - analytic: {b2_a}{(24 - len(str(b2_a))) * ' '}- numeric: {b2_n}")
            print(f"a3 - analytic: {a3_a}{(24 - len(str(a3_a))) * ' '}- numeric: {a3_n}")
            print(f"b3 - analytic: {b3_a}{(24 - len(str(b3_a))) * ' '}- numeric: {b3_n}")
            print(f"a4 - analytic: {a4_a}{(24 - len(str(a4_a))) * ' '}- numeric: {a4_n}")
            print(f"b4 - analytic: {b4_a}{(24 - len(str(b4_a))) * ' '}- numeric: {b4_n}")
            
            print("\np =", p_i)
            print(f"Torque is {close_torque_math}")
            print("M_numeric = ", M_numeric)
            print("M_analytic =", M_analytic)
    
            break
        if i % 100 == 0:
            print("")
    
#%% Debug
if False:
        x_analytic= [
        analytic_solution_K1(p_i, K1=K_i1, r1=r_i1, r2=r_i2, r3=r_i3,  
                                      mu_r1=mu_ir1, mu_r3=mu_ir3),
        analytic_solution_K2(p_i, K2=K_i2, r1=r_i1, r2=r_i2, r3=r_i3,  
                                  mu_r1=mu_ir1, mu_r3=mu_ir3)    
            ]
        
        model.total_torque()