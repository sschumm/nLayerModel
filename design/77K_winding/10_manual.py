# -*- coding: utf-8 -*-
#%% imports and functions
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model
from design import get_L_TPL2100 as L

mdl, plt, res = [None]*3
r_so, r_si, r_sA, r_rF, r_ro, r_ri = [0]*6

def update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r):
    global r_so, r_si, r_sA, r_rF, r_ro, r_ri
    r_si = r_so - h_yoke_s
    r_sA = r_si - 0.5 * h_wndg_s
    r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
    r_ro = r_rF - 0.5 * h_wndg_r
    r_ri = r_ro - h_yoke_r
    
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
    
    for iter_idx_lift in range(20): 
        
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
        pyplt.figure(dpi=1000)
        pyplt.plot([i for i in range(len(K_s_hist))], K_s_hist)
        pyplt.figure(dpi=1000)
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
update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
params = Main_Params(p, l_e, r_so, r_si, k_fill_r, k_fill_s, B_yoke_max)
dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

# --- derive initial current loadings ---
K_s = params.k_fill_s * sw.J_e * h_wndg_s
K_r = params.k_fill_r * fw.J_e * h_wndg_r

# --- build initial model ---
mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, Ks=K_s, Kr=K_r)

res.show()

K_s, K_r = find_K_with_L(dims, params, False)