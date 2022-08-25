# -*- coding: utf-8 -*-
#%% imports
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi

from modules.model import Model
from modules.layer import MagneticLayer, AirLayer, CurrentLoading
from modules.plot.radial import RadialMultiPlot
from modules.plot.plane import PlanePlot

from analytics.precalcs import kb, kd, kp, K, taup

from data.specs import Generator as gn
from data.specs import StatorWinding_Cu as sw
from data.specs import FieldWinding as fw

# define model
def iterate_model(l):    
    # ---------------- adjust params -------------------
    p = 20
    
    # selected dimensions
    l = l
    h_ys = 0.15 # [m]
    h_yr = 0.05 # [m]
    h_ws = 0.001 # [m]
    h_wr = 0.001 # [m]
    
    # derived dimensions
    d_so = np.sqrt(sw.d2so_L/l)
    r_so = d_so * 0.5
    r_si = r_so - h_ys
    r_ro = r_si - gn.delta_mag
    r_ri = r_ro - h_yr
    
    r_f = r_ro + 0.5 * h_wr
    r_a = r_si - 0.5 * h_ws
    
    tau_p = taup(d_si=r_si*2, p=p)
    # print(f"{tau_p/gn.delta_mag = }, should be > 3")
    
    # load angle
    alpha_f = pi * 0.5
    alpha_a = pi * 0.0
    
    # winding numbers
    Nf = 50
    Ns = 100
    
    # rotor winding factor
    kr_w = 1
    
    # stator winding factor
    ks_w = 1
    
    # currents
    I_f = 2e3
    I_a = 1e3
    
    # current loadings
    A_r = K(m=1, I=I_f, d=2*r_f, N=2*p*Nf)
    A_s = K(m=gn.m, I=I_a, d=2*r_a, N=Ns)
    
    A_r_amplitude = np.sqrt(2) * kr_w * A_r
    A_s_amplitude = np.sqrt(2) * ks_w * A_s
    # ---------------- build model ---------------------
    model = Model(p=p, l=l)
    model.add_layer(AirLayer(r=r_ri))
    model.add_layer(MagneticLayer(r=r_ro, 
                                  mu_r=gn.mu_r_yoke))
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
    model.add_layer(AirLayer(r=r_si))
    model.add_layer(MagneticLayer(r=r_so, 
                                  mu_r=gn.mu_r_yoke))
    model.build()
    model.solve()
    model.total_torque()
    return model, d_so


#%% vary generator length

lengths = list(np.linspace(1,5,100))
diams = list()
p_el_out = list()

for length in lengths:
    i_model, diam = iterate_model(length)
    P_el_out = gn.w_syn * i_model.Mneg / 1e6
    diams.append(diam)
    p_el_out.append(P_el_out)

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

plt.figure(figsize=(15, 15))
# c = plt.contourf(lengths, poles, p_el_out,
#                  levels=100)
# plt.colorbar(c)
ax = plt.axes(projection="3d")
ax.plot3D(lengths, diams, p_el_out)
ax.set_xlabel("Generator Length")
ax.set_ylabel("Generator Diameter")
ax.set_zlabel("Power Output")