# -*- coding: utf-8 -*-
#%% imports and functions
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, yoke_height, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model
from design import get_L_TPL2100 as L

plot_all = False
mdl, plt, res = [None]*3
r_so, r_si, r_sA, r_rF, r_ro, r_ri = [0]*6

def update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r):
    global r_so, r_si, r_sA, r_rF, r_ro, r_ri
    r_si = r_so - h_yoke_s
    r_sA = r_si - 0.5 * h_wndg_s
    r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
    r_ro = r_rF - 0.5 * h_wndg_r
    r_ri = r_ro - h_yoke_r
    return Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)
    
def find_K_with_L(dims: Main_Dims, params: Main_Params, verbose = False):
    
    h_wndg_s, h_wndg_r = 2*(dims.r_si-dims.r_sA), 2*(dims.r_rF-dims.r_ro)
    
    K_s_static = params.k_fill_s * h_wndg_s * sw.J_e
    K_r_static = params.k_fill_r * h_wndg_r * fw.J_e
    
    mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, 
                                         Ks=K_s_static, Kr=K_r_static)
    
    J_e_s = L(sw.T_HTS, res.B_s_c, 30) * gn.Ic_spec/(gn.A_tape*1.7)
    J_e_r = L(fw.T_HTS, res.B_r_c, 30) * gn.Ic_spec/(gn.A_tape*1.7)    
    K_s_lift = params.k_fill_s * J_e_s * h_wndg_s
    K_r_lift = params.k_fill_r * J_e_r * h_wndg_r
    
    K_s_hist = [K_s_static, K_s_lift]
    K_r_hist = [K_r_static, K_r_lift]
    
    for idx_lift in range(20): 
        
        K_s = K_s_hist[-2] * 0.3 + K_s_hist[-1] * 0.7
        K_r = K_r_hist[-2] * 0.3 + K_r_hist[-1] * 0.7
        
        mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, 
                                             Ks=K_s, Kr=K_r)
        
        J_e_s = L(sw.T_HTS, res.B_s_c, 30) * gn.Ic_spec/(gn.A_tape*1.7)
        J_e_r = L(fw.T_HTS, res.B_r_c, 30) * gn.Ic_spec/(gn.A_tape*1.7)    
        K_s_lift = params.k_fill_s * J_e_s * h_wndg_s
        K_r_lift = params.k_fill_r * J_e_r * h_wndg_r     
        K_s_hist.append(K_s_lift)
        K_r_hist.append(K_r_lift)
        
    # --- plot convergence ---
    if verbose:    
        pyplt.figure(dpi=500)
        pyplt.plot([i for i in range(len(K_s_hist))], K_s_hist)
        pyplt.figure(dpi=500)
        pyplt.plot([i for i in range(len(K_r_hist))], K_r_hist)
    return K_s, K_r
    
#%% initialize parameters
# --- define global parameters --- 
p = 20
l_e = 0.4 # [m]
r_so = np.sqrt(sw.d2so_L/l_e)/2 # [m]
k_fill_r, k_fill_s = 0.5, 0.5
B_yoke_max = 6 # [T]

# --- define inital dimensions ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = gn.delta_mag, gn.delta_mag
dims = update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
params = Main_Params(p, l_e, r_so, r_si, k_fill_r, k_fill_s, B_yoke_max)
dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

# --- derive initial current loadings ---
K_s = params.k_fill_s * sw.J_e * h_wndg_s
K_r = params.k_fill_r * fw.J_e * h_wndg_r

# --- build initial model ---
mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, Ks=K_s, Kr=K_r)
res.show("initial model")
if plot_all: plt.fluxplot(1000, 1000, lvls=10)

# --- initial model with lift factor applied ---
K_s, K_r = find_K_with_L(dims, params, False)
mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, Ks=K_s, Kr=K_r)
res.show("initial model with L applied")
if plot_all: plt.fluxplot(1000, 1000, lvls=10)

#%% adapt yokes
# --- adapt stator yoke ---
h_yoke_s = yoke_height(params.p, r_si, res.B_s, params.B_yoke_max)

# --- adapt rotor yoke ---
h_yoke_r = yoke_height(params.p, r_si, res.B_r, params.B_yoke_max)

# --- build adapted model ---
dims = update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, Ks=K_s, Kr=K_r)
res.show("initial model with L applied")
if plot_all: plt.fluxplot(1000, 1000, lvls=10)

#%% vary the current loading of the stator winding
# --- increase the winding height of the stator in a for loop ---
increasing_h_by_factor=0.2
break_when_P_at_factor=0.5

iter_s_P, iter_s_h = [res.P], [h_wndg_s]
for idx_h_wndg_s in range(30):
    h_wndg_s = iter_s_h[-1] + iter_s_h[0] * increasing_h_by_factor
    dims = update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
    K_s, K_r = find_K_with_L(dims, params, False)
    mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, Ks=K_s, Kr=K_r)
    # res.show_P()
    iter_s_P.append(res.P)
    iter_s_h.append(h_wndg_s)
    if iter_s_P[-1] < iter_s_P[0] * break_when_P_at_factor:
        break

# --- plot the results ---
if True:
    fig, ax1 = pyplt.subplots()
    ax2 = ax1.twinx()
    fig.tight_layout() 
    fig.dpi=500
    ax1.scatter([x for x in range(len(iter_s_P))], [y*1e-6 for y in iter_s_P], marker=".", c="b")
    ax1.set_ylabel('P [MW]', color="b")
    ax2.plot([x for x in range(len(iter_s_h))], iter_s_h, c="r")
    ax2.set_ylabel('h_wndg_s [m]', color="r")


# --- select an index ---
idx = 2






