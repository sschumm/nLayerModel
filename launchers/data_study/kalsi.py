# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot
from analytics.precalcs import kb, kd, K



#%% -------------------- init parameters --------------------

# [-]
Q = 288
m = 3
p = 12
q = 4
a = 1 # probably due to high voltage of 4.5 [kV]
mu_r = 1e4
Nf_pc = 800 # 800 or 600 - ??? in text and table different values are given
Na_pc = 11


# [A]
I_f = 1000
I_a = 1360


# [m]
lfe = 0.69 # < 2
rri = 2.072
rro = 2.170
r_f = 2.22 + (2.238 - 2.22)/2
r_a = 2.368+ (2.375 -2.368)/2
rsi = 2.426
rso = 2.5

# Rotor Winding Factor
sigma_r = 2/3
kr_b = kb(n=1, sigma=sigma_r)

# Stator Winding Factor
ks_d = kd(m=m, q=q)

# [A/m]
A_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf_pc)   # single-layer winding in rotor
A_s = K(m=m, I=I_a, d=2*r_a, N=q*p*Na_pc/a) # single-layer winding in stator

A_r_amplitude = np.sqrt(2) * A_r * kr_b
A_s_amplitude = np.sqrt(2) * A_s * ks_d 


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
P_out_target = 10.5
P_out_model = w_rated * model.Mneg / 1e6

print(f"{P_out_target = } [MW]")
print(f"{P_out_model  = } [MW]")

# [MN]
Torque = model.Mneg/1e6
print(f"{Torque = } [MN]")

# [T]
B = model.get_B_data(np.linspace(rri, rso, 1000), np.linspace(0, 2*pi, 1000))
Bmax = np.max(np.sqrt(B.Br**2 + B.Bt**2))
print(f"{Bmax = } [T]")

#%% -------------------- create plots --------------------
plt = PlanePlot(model, fgsz=70)
# plt.contour(dr=1000, dt=1000, style="jet")
# plt.fluxplot(dr=1000, dt=1000, lvls=15)




