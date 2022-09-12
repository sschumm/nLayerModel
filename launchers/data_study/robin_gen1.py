# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot
from analytics.precalcs import kb, kd, kp, K



#%% -------------------- init parameters --------------------

# [-]
Q = 468
m = 3
p = 39
q = 2
a  = 156
mu_r = 1000
Nf = 63
Ns = 89


# [A]
I_f = 416.897
I_a = 6.1654 * 1e3


# [m]
hp = 0.0481
hyr = 0.0462
tau_p = 0.2485 
wp = 0.0982 
wp_side = 0.012

l_i = 1.2003 # lfe = 1.1173 # L = 1.2023 
rri = 2.9532 
rro = rri + hyr + hp # == 3.0475 
r_f = rro 
# rdelta = 3.0661
r_a = 3.0847
rsi = 3.0847
rso = 3.2543
 

# Rotor Winding Factor
# breadth factor
sigma_r = (tau_p - wp - wp_side) / tau_p
kr_b = kb(n=1, sigma=sigma_r)


# Stator Winding Factor
ks_d = kd(m=m, q=q) # zone factor 
ks_p = kp(w=5, tp=6) # pitching factor

# [A/m]
# A_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf) # 1 layer field winding
A_r = K(m=1, I=I_f, d=2*r_f, N=2*2*p*Nf) # 2 layer field winding
A_s = K(m=m, I=I_a, d=2*r_a, N=Ns) # given: A_s=1.6987*1e5 [A/m]

A_r_amplitude = np.sqrt(2) * A_r * kr_b
A_s_amplitude = np.sqrt(2) * A_s * ks_d * ks_p # np.sin((pi/2) * (5/6))


n_syn = 8.33   # [rpm]

# [rad/s]
w_syn = 0.8723 
gamma = 2.9044
alpha_f = pi * 0.5
alpha_a = pi * 0.0


#%% -------------------- model setup --------------------

model = Model(p=p, l=l_i)

# Air -| |- Iron -| |- Air -|k|- Air -|k|s- Iron -||- Env 

model.add_layer(AirLayer(r=rri))
model.add_layer(MagneticLayer(r=rro, 
                              mu_r=mu_r))
model.add_layer(CurrentLoading(K=A_r_amplitude, 
                               r=r_f, 
                               alpha=alpha_f, 
                               mu_r=1
                               ))
model.add_layer(CurrentLoading(K=A_s_amplitude,
                               r=r_a,
                               alpha=alpha_a,
                               mu_r=1.
                               ))
# model.add_layer(AirLayer(r=rsi))
model.add_layer(MagneticLayer(r=rso, 
                              mu_r=mu_r))

model.build()
model.solve()
model.total_torque()

#%% -------------------- evaluation --------------------

print(f"{A_r = } [A/m]")
print(f"{A_s = } [A/m]")
print(f"{A_r_amplitude = } [A/m]")
print(f"{A_s_amplitude = } [A/m]")

# [MW]
P_out_target = -7
P_out_model = w_syn * model.Mneg / 1e6

print(f"{P_out_target = } [MW]")
print(f"{P_out_model  = } [MW]")

# [MN]
Torque = model.Mpos/1e6
print(f"{Torque = } [MN]")

# [T]
B = model.get_B_data(np.linspace(rri, rso, 1000), np.linspace(0, 2*pi, 1000))
Bmax = np.max(np.sqrt(B.Br**2 + B.Bt**2))
print(f"{Bmax = } [T]")

#%% -------------------- create plots --------------------
plt = PlanePlot(model, fgsz=150)
# plt.contour(dr=1000, dt=1000, style="jet")
# plt.fluxplot(dr=1000, dt=1000, lvls=15)

rmp = RadialMultiPlot(model)
# rmp.multiplot(["Az", "Br", "Ht"])

#%% angle sweep
if False:
    print("INFO: computing load angle sweep...")
    x_angle = (np.pi/2) * np.linspace(0, 2, 100)
    y_torque = list()
    alpha_a = pi * 0.0

    for i in x_angle:
        alpha_f = i
    
        model = Model(p=p, l=l_i)
        
        # Air -| |- Iron -| |- Air -|k|- Air -|k|s- Iron -||- Env 
        
        model.add_layer(AirLayer(r=rri))
        model.add_layer(MagneticLayer(r=rro, 
                                      mu_r=mu_r))
        model.add_layer(CurrentLoading(K=A_r_amplitude, 
                                       r=r_f, 
                                       alpha=alpha_f, 
                                       mu_r=1
                                       ))
        model.add_layer(CurrentLoading(K=A_s_amplitude,
                                       r=r_a,
                                       alpha=alpha_a,
                                       mu_r=1.
                                       ))
        # model.add_layer(AirLayer(r=rsi))
        model.add_layer(MagneticLayer(r=rso, 
                                      mu_r=mu_r))
        
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