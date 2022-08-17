# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot
from data.precalculations import kb, kd, K



#%% -------------------- init parameters --------------------

# [-]
Q = 468
m = 3
p = 39
q = 2
a  = 156
mu_r = 1e6
sigma = 0.5 # sigma=0 -> breadth=0, sigma=1 -> breadth=pi
Nf = 63
Ns = 89


# [A]

# --!--
I_f = 416.897
I_f = 2.8644 * 1e4
# ----

I_a = 6.1654 * 1e3


# [m]
rdelta = 3.0661
lfe = 1.1173  # l_i = 1.2003 # L = 1.2023 
rri = 2.9532 
rro = 3.0475 
r_f = rdelta
r_a = 3.0847
rsi = 3.0847
rso = 3.2543


# [-]
k_d = kd(m=m, q=q) 
k_b = kb(n=1, sigma=sigma)

# [A/m]

# --!--
A_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf)
A_r = K(m=1, I=I_f, d=2*r_f, N=Nf)
A_r = K(m=1, I=I_f, d=2*r_f, N=2*Nf)
# ----

A_s = K(m=m, I=I_a, d=2*r_a, N=Ns) # given: A_s=1.6987*1e5 [A/m]

A_r_amplitude = np.sqrt(2) * A_r
A_s_amplitude = np.sqrt(2) * A_s * k_d * k_b


n_syn = 8.33   # [rpm]

# [rad/s]
w_syn = 0.8723 
gamma = 2.9044
alpha_f = pi * 0.5 # gamma * 0.5 
alpha_a = pi * 0.0

# rdelta = 3.0661, 


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
P_out_target = 7
P_out_model = w_syn * model.M / 1e6

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

rmp = RadialMultiPlot(model)
# rmp.multiplot(["Az", "Br", "Ht"])

