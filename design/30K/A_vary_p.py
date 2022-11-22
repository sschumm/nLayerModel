# -*- coding: utf-8 -*-
import sys
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_30K as sw
from data import StatorWinding_77K as tw
from data.export import save_params

from design import n7_Model


def get_k_fill_s(r_si, h_sw, p):
    nom = (pi*(r_si-h_sw)/(3*p)-gn.w_pole_frame-gn.r_bend_max)*(h_sw-2*gn.h_pole_frame)*3*p
    den = h_sw * (r_si - h_sw) * pi
    return nom / den


p = 20

d2l = 1

r_so = 3.75
l_e = d2l / (2*r_so)**2 # sw.l_e(sw, r_so) # 0.2845
print(f"{l_e =}")


generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     k_fill_s=0.18,
                     k_fill_r=0.2,
                     B_yoke_max=1.7)


# --- init generator dimension via yoke and winding heights ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.4, 1.5]]

generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
    
# --- compute initial current loadings and with it the initial model ---
K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.866)

# generator.show_results()
# generator.fast_flux()

# generator.keep_const_kr_b(kr_b = 0.6)
# generator.apply_coil_sizes_and_lift_factor(verbose = False)




#%%
generator.keep_const_kr_b(kr_b=0.8)
generator.apply_coil_sizes_and_lift_factor(verbose=False)
generator.k_fill_s = get_k_fill_s(generator.dims.r_si, generator.h_wndg_s, generator.p)
generator.coil_shapes(verbose=True)


lst_P = []
lst_K = []
lst_h_wndg_s = []
lst_h_wndg_r = []

lst_rs_bend = []
lst_rr_bend = []

h_wndg_s_0 = 2 * gn.h_pole_frame * 1.3
h_wndg_s_1 = h_wndg_s_0 * 2

h_wndg_r_0 = 2 * gn.h_pole_frame * 1.1
h_wndg_r_1 = h_wndg_r_0 * 2


n_iters_s = 9
n_iters_r = 20 

for idx_stator, iter_stator in enumerate(np.linspace(0, h_wndg_s_1-h_wndg_s_0, n_iters_s)):
       
    h_wndg_s = h_wndg_s_0 + iter_stator
    generator.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=True)
        
    hist_r_P = []
    hist_r_K = []
    hist_r_h = []
    hist_rs_bend = []
    hist_rr_bend = []
    for idx_rotor, iter_rotor in enumerate(np.linspace(0, h_wndg_r_1-h_wndg_r_0, n_iters_r)):
        
        h_wndg_r = h_wndg_r_0 + iter_rotor
        
        generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
        generator.k_fill_s = get_k_fill_s(generator.dims.r_si, generator.h_wndg_s, generator.p)
        generator.apply_coil_sizes_and_lift_factor(verbose=False)
        
        # generator.coil_shapes(verbose=True)
        
        hist_r_P.append(generator.P)
        hist_r_K.append(generator.K_s)
        # hist_r_h.append(generator.h_wndg_r)
        if idx_stator == 0:
            lst_h_wndg_r.append(h_wndg_r)
        generator.coil_shapes()
        hist_rs_bend.append(generator.coil.r_s_bend)
        hist_rr_bend.append(generator.coil.r_r_bend)
        
        sys.stdout.write(f"\rStator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r}    ")
        sys.stdout.flush()
        
    lst_P.append(hist_r_P)
    lst_K.append(hist_r_K)
    lst_h_wndg_s.append(h_wndg_s)
    lst_rs_bend.append(hist_rs_bend)
    lst_rr_bend.append(hist_rr_bend)
        
        
#%%    
for iters_s in range(n_iters_s):
    fig = plt.figure(dpi=300, figsize=(6,7))
    ax1 = plt.subplot()
    plt.grid()
    plt.title(f"{lst_h_wndg_s[iters_s] = }")
    ax1.scatter(lst_h_wndg_r, [v*1e-6 for v in lst_P[iters_s]])
    ax1.plot(lst_h_wndg_r, [v*1e-6 for v in lst_P[iters_s]])
    
    ax2 = ax1.twinx()
    # ax2.set_ylim(0.015, 0.08)
    ax2.scatter(lst_h_wndg_r, lst_rs_bend[iters_s], color = "red")
    ax2.scatter(lst_h_wndg_r, lst_rr_bend[iters_s], color = "green")
    ax2.plot(lst_h_wndg_r, [generator.gn.r_bend_max] * len(lst_h_wndg_r), color="black")
    plt.show()
    pass

#%%

# generator.fast_flux()
