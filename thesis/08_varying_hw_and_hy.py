# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from data.export import save_params

from design import n7_Model


# --- init generator with main parameter ---
p = 24
l_e = 0.3
r_so = sw.r_so(sw, l_e)


generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     k_fill_s=0.3,
                     k_fill_r=0.3,
                     B_yoke_max=1.7)


# --- init generator dimension via yoke and winding heights ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r =  0.01 * r_so, 0.02 * r_so #[factor * (gn.h_pole_frame*2) for factor in [1.2, 1.2]]

generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
    
# --- compute initial current loadings and with it the initial model ---
K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r

generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_p=0.866)
generator.show_results()
# generator.fast_flux()

#%%

h_wndg_s_0 = 0
h_wndg_s_step = 0.01
history_K_s = []
history_h_wndg_s = []
history_P = []

for i in range(30):
    h_wndg_s = h_wndg_s_0 + h_wndg_s_step * i
    
    generator.update_dimensions(h_wndg_s=h_wndg_s)
    K_s = generator.k_fill_s * generator.J_e_s * h_wndg_s
    generator.update_model_by_K(K_s = K_s, K_r = K_r)
    
    history_h_wndg_s.append(h_wndg_s)
    history_K_s.append(K_s)
    history_P.append(generator.P)
    
    
fig = plt.figure(dpi=300, figsize=(6,5))
ax1 = plt.subplot()
ax2 = ax1.twinx()
# ax.set_xlabel("B / T")
# ax.set_ylabel("J_e / Amm2")
# ax.set_xlim(0, 3)
# ax.set_ylim(0, 25)
# ax.set_xticks([i for i in range(0, 181, 20)])
# ax.set_yticks([i for i in range(0, 2501, 250)])
plt.grid()

ax1.plot(history_h_wndg_s, history_K_s) #, label=f"T = {lst_T[0]} / K")
ax2.plot(history_h_wndg_s, history_P) #, label=f"T = {lst_T[0]} / K")
# generator.fast_flux()