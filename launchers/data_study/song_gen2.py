# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=5)
import tikzplotlib as tkz
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
a  = 1
mu_r = 1e6
Nf_pp = 3220 # =/= gen1
Na_ps = 22


# [A]
I_f = 340 # =/= gen1
I_a = 760


# [m]
lfe = 0.8 
rri = 1.55 # =/= gen1
rro = 1.632 # =/= gen1
r_f = 1.736 # =/= gen1
r_a = 1.92
rsi = 1.974
rso = 2.05

# Rotor Winding Factor
# breadth factor, params taken from fig. 2
coil_width = 0.12 # [m]
coil_side_distance = 0.18 # [m]
sigma_r = coil_width*2/(r_f*pi/p) # coil_width / (coil_width+coil_side_distance)
kr_b = kb(n=1, sigma=sigma_r)

# Stator Winding Factor
# zone factor == given winding factor: 0.9577
ks_d = kd(m=m, q=q)

# [A/m]
A_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf_pp)
A_s = K(m=m, I=I_a, d=2*r_a, N=q*p*Na_ps) # given: A_s = 400 [kA/m]

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
P_out_target = 10 
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


#%% eval with varying r_f

def f(K_r_amplitude, r_rw):
    model = Model(p=p, l=lfe)
    model.add_layer(AirLayer(r=rri))
    model.add_layer(MagneticLayer(r=rro, 
                                  mu_r=mu_r))
    model.add_layer(CurrentLoading(K=K_r_amplitude, 
                                    r=r_rw, 
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
    return model.Mpos
    
x = np.linspace(0.9521640091116172, 1.0478359908883828, 100) #(rro -> r_rw_außen)
y = []

for i in x:
    r_rw = r_f * i
    sigma_r = coil_width*2/(r_f*pi/p)
    kr_b = kb(n=1, sigma=sigma_r)
    K_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf_pp)
    K_r_amplitude = np.sqrt(2) * A_r * kr_b
    
    P = w_rated * f(K_r_amplitude, r_rw) * 1e-6
    y.append(P)
    
import matplotlib.pyplot as pyplot
fig = pyplot.figure(dpi=300)
ax = pyplot.subplot()
pyplot.grid()
pyplot.plot(x, y)

tkz.clean_figure()
tkz.save("msm_vary_r_rw.tex")



#%% varying alpha

def g(alpha_f):
    model = Model(p=p, l=lfe)
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
    return model.Mpos
    
x = np.linspace(pi/4, pi/2, 100)
y = []

for i in x:
    a = i
        
    P = w_rated * g(a) * 1e-6
    y.append(P)
    
import matplotlib.pyplot as pyplot
fig = pyplot.figure(dpi=300)
ax = pyplot.subplot()
pyplot.grid()
pyplot.plot(x, y)


#%% eval with varying w

def h(K_r_amplitude):
    model = Model(p=p, l=lfe)
    model.add_layer(AirLayer(r=rri))
    model.add_layer(MagneticLayer(r=rro, 
                                  mu_r=mu_r))
    model.add_layer(CurrentLoading(K=K_r_amplitude, 
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
    return model.Mpos
    
x = np.linspace(0.95, 1.03, 100)
y = []

for i in x:
    width = coil_width * i
    sigma_r = width*2/(r_f*pi/p)
    kr_b = kb(n=1, sigma=sigma_r)
    K_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf_pp)
    K_r_amplitude = np.sqrt(2) * A_r * kr_b
    
    P = w_rated * h(K_r_amplitude) * 1e-6
    y.append(P)
    
import matplotlib.pyplot as pyplot
fig = pyplot.figure(dpi=300)
ax = pyplot.subplot()
pyplot.grid()
pyplot.plot(x, y)


#%% eval with varying r_a

def k(r_sw):
    model = Model(p=p, l=lfe)
    model.add_layer(AirLayer(r=rri))
    model.add_layer(MagneticLayer(r=rro, 
                                  mu_r=mu_r))
    model.add_layer(CurrentLoading(K=A_r_amplitude, 
                                    r=r_f, 
                                    alpha=alpha_f, 
                                    mu_r=1
                                    ))
    model.add_layer(CurrentLoading(K=A_s_amplitude,
                                    r=r_sw,
                                    alpha=alpha_a,
                                    mu_r=1.
                                    ))
    model.add_layer(AirLayer(r=rsi))
    model.add_layer(MagneticLayer(r=rso, 
                                  mu_r=mu_r))

    model.build()
    model.solve()
    model.total_torque()
    return model.Mpos
    
x = np.linspace(0.95, 1.07, 100)
y = []

for i in x:
    r_sw = r_a * i
    
    P = w_rated * k(r_sw) * 1e-6
    y.append(P)
    
import matplotlib.pyplot as pyplot
fig = pyplot.figure(dpi=300)
ax = pyplot.subplot()
pyplot.grid()
pyplot.plot(x, y)
