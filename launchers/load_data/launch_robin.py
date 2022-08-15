# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)

from scipy.constants import pi
from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from data.precalculations import K


# -------------------- init parameters --------------------
Q = 468
m = 3
p = 39
q = 2

Nf = 63
Ns = 89
Nc = 89
a  = 156

lfe = 1.1173  # l_i = 1.2003 # L = 1.2023
rri = 2.9532
rro = 3.0475
# rdelta = 3.0661
rsi = 3.0847
rso = 3.2543

I_f = 416.897
I_s = 6.1654 * 1e3

n_syn = 8.33   # [rpm]
w_syn = 0.8723 # [1/s]
P_out_target = 7e6 # -6938312.88102393

A_r = K(m=1, I=I_f, d=2*rro, Ns=Nf)
A_s = K(m=m, I=I_s, d=2*rsi, Ns=Ns) # A_s = 1.6987 * 1e5

A_r_amplitude = np.sqrt(2) * A_r
A_s_amplitude = np.sqrt(2) * A_s


# -------------------- model setup --------------------

model = Model(p=p, l=lfe)

model.add_layer(AirLayer(r=rri))
model.add_layer(CurrentLoading(K=A_r_amplitude, 
                               r=rro, 
                               alpha=pi*0.5, 
                               mu_r=1e5
                               ))
model.add_layer(CurrentLoading(K=A_s_amplitude,
                               r=rsi,
                               alpha=pi*0.0,
                               mu_r=1.
                               ))
model.add_layer(MagneticLayer(r=rso, 
                              mu_r=1e5
                              ))

model.build()
model.solve()
model.total_torque()

# -------------------- post processing --------------------
P_out_model = w_syn * model.M
print("P_out_target =", P_out_target)
print("P_out_model  =", P_out_model)