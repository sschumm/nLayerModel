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
Q = 348
m = 3
p = 29
q = 2
a  = 29
mu_r = 1000
Nf = 103
Nc = 27
Ns = 108


# [A]
I_f = 1.1527 * 1e3
I_a = 5.9685 * 1e3


# [m]
tau_p = 0.3322 
tau_Q = 0.0554
hyr = 0.0757
hys = 0.0635

hp = 0.0481
wp = 0.139 
wp_side = 0.012

hQ = 0.0827 
bQ = 0.0471



l_i = 0.7691 # lfe = 0.686 # L = 0.771 
rri = 2.9054 
rro = rri + hyr # == 2.9811
r_f = rro + hp
# rdelta = 3.0478
rso = 3.2126
rsi = rso - hys # given: 3.0664 - but:   rso-hys == 3.1491
r_a = rsi - hQ/2 # given: 3.0664 - but:  rso-hys-hQ/2 == 3.10775


# Rotor Winding Factor
sigma_r = bQ / tau_Q
kr_b = kb(n=1, sigma=sigma_r) # breadth factor


# Stator Winding Factor
ks_d = kd(m=m, q=q) # zone factor 
ks_p = kp(w=6, tp=6) # pitching factor

# [A/m]
A_r = K(m=1, I=I_f, d=2*r_f, N=2*2*p*Nf) 
A_s = K(m=m, I=I_a, d=2*r_a, N=Ns) # given: A_s=2.0074*1e5 [A/m]

A_r_amplitude = np.sqrt(2) * A_r * kr_b
A_s_amplitude = np.sqrt(2) * A_s * ks_d * ks_p 


n_syn = 8.33   # [rpm]

# [rad/s]
w_syn = 0.8723 
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

# print(f"{A_r = } [A/m]")
# print(f"{A_s = } [A/m]")
# print(f"{A_r_amplitude = } [A/m]")
# print(f"{A_s_amplitude = } [A/m]")

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
# print(f"{Bmax = } [T]")

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