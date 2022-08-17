# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot
from data.precalculations import kb, kd, K



#%% -------------------- init parameters --------------------

# [-]
Q = 288
m = 3
p = 12
q = 4
a  = 1
mu_r = 1e6
sigma = 0.5
Nf_pp = 1400
Na_ps = 22


# [A]
I_f = 680
I_a = 760


# [m]
lfe = 0.8 
rri = 1.6
rro = 1.672
r_f = 1.756
r_a = 1.92
rsi = 1.974
rso = 2.05


# [-]
k_d = kd(m=m, q=q) 
k_b = kb(n=1, sigma=sigma)

# [A/m]
A_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf_pp)
A_s = K(m=m, I=I_a, d=2*r_a, N=q*p*Na_ps) # given: A_s = 400 [kA/m]

A_r_amplitude = np.sqrt(2) * A_r
A_s_amplitude = np.sqrt(2) * A_s * k_d * k_b


n_rated = 10   # [rpm]

# [rad/s]
w_rated = 2*pi*n_rated/60 
alpha_f = pi * 0.5
alpha_a = pi * 0.0




#%% -------------------- model setup --------------------

model = Model(p=p, l=lfe)

# Air -- Iron -- Air -|- Air -|- Air -- Iron -- Air 

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
model.add_layer(AirLayer(r=rsi))
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
P_out_target = 10 
P_out_model = w_rated * model.M / 1e6

print(f"{P_out_target = } [MW]")
print(f"{P_out_model  = } [MW]")

# [MN]
Torque = model.M/1e6
print(f"{Torque = } [MN]")

# B = model.get_B_data(np.linspace(rri, rso, 1000), np.linspace(0, 2*pi, 1000))
# Bmax = np.max(B.Bt)

# print(f"{Bmax = } [T]")

#%% -------------------- create plots --------------------
plt = PlanePlot(model, fgsz=70)
# plt.contour(dr=1000, dt=1000, style="jet")



