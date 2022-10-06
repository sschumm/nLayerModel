# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model
from design import get_L_TPL2100 as L

K_s, K_r = 0, 0
Pel_goal = gn.Pel_out * 1
tolerance_lower, tolerance_upper = 0.95, 1.1
A_tape = (12*0.22) * 1e-6 # [m**2]
Ic_0 = 300 # [A]
r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000
all_dims, all_models = [], []
r_si, r_sA, r_rF, r_ro, r_ri = 0,0,0,0,0

p = 21
l_e = 0.4 # [m]
r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)

h_yoke_s = 0.02 * r_so
h_wndg_s = gn.delta_mag
h_wndg_r = gn.delta_mag
h_yoke_r = h_yoke_s

def update_dimensions():
    global r_si, r_sA, r_rF, r_ro, r_ri
    r_si = r_so - h_yoke_s
    r_sA = r_si - 0.5 * h_wndg_s
    r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
    r_ro = r_rF - 0.5 * h_wndg_r
    r_ri = r_ro - h_yoke_r

update_dimensions()


""" --- initial params --- """
init_params = Main_Params(p, l_e, r_so, r_si, k_fill_r=0.5, k_fill_s=0.5, B_yoke_max=6)
init_dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

J_e_s_init = sw.J_e * 0.9 # reduce given vals by factor
J_e_r_init = fw.J_e * 0.1 # reduce given vals by factor 

K_s_init = init_params.k_fill_s * h_wndg_s * J_e_s_init
K_r_init = init_params.k_fill_r * h_wndg_r * J_e_r_init

""" --- find model with static K_s and K_r --- """
if True:
    K_s = K_s_init
    K_r = K_r_init
    iter_static_P_history = []
    for iter_idx_static in range(50):
        
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
        iter_static_P_history.append(res.P)
        
        # --- adapt stator current loading ---
        if res.P < tolerance_upper * Pel_goal and res.P > tolerance_lower * Pel_goal:
            res.show(f"{iter_idx_static = }")
            break
        elif res.P < Pel_goal:
            K_s *= 1.1
            h_wndg_s = K_s / (init_params.k_fill_s * sw.J_e)
            update_dimensions()
        else:
            K_s *= 0.9
            h_wndg_s = K_s / (init_params.k_fill_s * sw.J_e)
            update_dimensions()
               
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
        iter_static_P_history.append(res.P)
        
        # --- adjust stator yoke ---        
        flux_p_pole = Phi_from(res.B_airgap, 
                                taup=taup(2*r_si, init_params.p), 
                                l_e=init_params.l_e)
        h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
        update_dimensions()
        
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
        iter_static_P_history.append(res.P)
        
        # --- adapt rotor current loading ---
        if res.P < tolerance_upper * Pel_goal and res.P > tolerance_lower * Pel_goal:
            res.show(f"{iter_idx_static = }")
            break
        elif res.P < Pel_goal:
            K_r *= 1.1
            h_wndg_r = K_r / (init_params.k_fill_r * fw.J_e)
            update_dimensions()
        else:
            K_r *= 0.9
            h_wndg_r = K_r / (init_params.k_fill_r * fw.J_e)
            update_dimensions()
        
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
        iter_static_P_history.append(res.P)
        
        # --- detect diverging loadings ---
        if iter_idx_static > 20 and iter_static_P_history[-1] < Pel_goal * 0.2:
            print(f"detected possible divergence of loadings in {iter_idx_static = }")
            pyplt.figure(dpi=1000)
            pyplt.plot([i for i in range(len(iter_static_P_history))], iter_static_P_history)
            break
        
    # --- adjust rotor yoke ---
    flux_p_pole = Phi_from(res.B_airgap, 
                           taup=taup(2*r_si, init_params.p), 
                           l_e=init_params.l_e)
    h_yoke_r = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
    update_dimensions()
    
    # --- build model with adjusted yokes ---
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    
    # res.show("model with adjusted yokes")
    # plt.fluxplot(dr,dt,lvls=10)
    ...

K_s_static = K_s
K_r_static = K_r

""" --- function to find K_s and K_r to fit critical current --- """
def find_K_with_L(K_s_static, K_r_static, K_s_lift, K_r_lift, verbose = False):
    iter_lift_K_s_history = [K_s_static] * 9 + [K_s_lift]
    iter_lift_K_r_history = [K_r_static] * 9 + [K_r_lift]
    
    for iter_idx_lift in range(100): 
        div = 10
        K_s = sum(iter_lift_K_s_history[-div:]) / div
        K_r = sum(iter_lift_K_r_history[-div:]) / div
        
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
          
        K_s_lift = init_params.k_fill_s * L(sw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7) * h_wndg_s
        K_r_lift = init_params.k_fill_r * L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7) * h_wndg_r        
        iter_lift_K_s_history.append(K_s_lift)
        iter_lift_K_r_history.append(K_r_lift)
    # --- plot convergence ---
    if verbose:    
        pyplt.figure(dpi=1000)
        pyplt.plot([i for i in range(len(iter_lift_K_s_history))], iter_lift_K_s_history)
        pyplt.plot([i for i in range(len(iter_lift_K_r_history))], iter_lift_K_r_history)
    return K_s, K_r

""" --- find model with B dependent engineering current density --- """
if True:
    
    # --- find densities and loadings for B_airgap from lift factor ---
    J_e_s = L(sw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)
    J_e_r = L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)    
    K_s_lift = init_params.k_fill_s * J_e_s * h_wndg_s
    K_r_lift = init_params.k_fill_r * J_e_r * h_wndg_r
    
    K_s, K_r = find_K_with_L(K_s_static, K_r_static, K_s_lift, K_r_lift)
        
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)  
    res.show("result after using the lift factor")
    # plt.quiver(dr=20, dt=200, scale=250, width=0.001)
    # plt.fluxplot(dr, dt, lvls=10)
    ...


""" --- increase h_wndg_r and apply lift factor --- """
if True:
    print("")
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    
    J_e_s = L(sw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)
    J_e_r = L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)   
    
    for iter_dyn in range(10):
                
        update_dimensions()
        
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
        
        J_e_s = L(sw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)
        J_e_r = L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)    
        K_s_lift = init_params.k_fill_s * J_e_s * h_wndg_s
        K_r_lift = init_params.k_fill_r * J_e_r * h_wndg_r
        
        K_s, K_r = find_K_with_L(K_s, K_r, K_s_lift, K_r_lift, verbose=False)
            
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r) 
        
        res.show_P(f"model with increasing {np.round(h_wndg_r, 3) = }") 
        
        # h_wndg_s *= 1.05
        h_wndg_r *= 1.05
        K_s = init_params.k_fill_s * J_e_s * h_wndg_s
        K_r = init_params.k_fill_r * J_e_r * h_wndg_r
        
    # res.show_B_airgap("")
    # K_r = init_params.k_fill_r * 
    ...



""" --- find model with B dependent engineering current density --- """
if False:
    
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    
    # --- find densities and loadings for B_airgap from lift factor ---
    J_e_s = L(sw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)
    J_e_r = L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)    
    K_s_lift = init_params.k_fill_s * J_e_s * h_wndg_s
    K_r_lift = init_params.k_fill_r * J_e_r * h_wndg_r
    
    K_s, K_r = find_K_with_L(K_s, K_r, K_s_lift, K_r_lift, verbose=False)
        
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)  
    
    res.show("after increasing h_wndg_r")
    # plt.quiver(dr=20, dt=200, scale=250, width=0.001)
    # plt.fluxplot(dr, dt, lvls=10)
    ...
    


""" --- sth --- """
if False:
    for iter_idx_dynamic in range(14):
        
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_r)
        J_e_r = L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7) 
        res.show_P(f"{iter_idx_dynamic = }")
        
        # --- adapt rotor current loading ---
        if res.P < tolerance_upper * Pel_goal and res.P > tolerance_lower * Pel_goal:
            res.show(f"{iter_idx_dynamic = }")
            break
        elif res.P < Pel_goal:
            h_wndg_r *= 1.15
            K_r = h_wndg_r * init_params.k_fill_r * J_e_r
            update_dimensions()
        else:
            h_wndg_r *= 0.95
            K_r = h_wndg_r * init_params.k_fill_r * J_e_r
            update_dimensions()

        
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    
    J_e_s = L(sw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)
    J_e_r = L(fw.T_HTS, res.B_airgap, 30) * Ic_0/(A_tape*1.7)    
    K_s_lift = init_params.k_fill_s * J_e_s * h_wndg_s
    K_r_lift = init_params.k_fill_r * J_e_r * h_wndg_r
    
    K_s, K_r = find_K_with_L(K_s, K_r, K_s_lift, K_r_lift, verbose=True)

    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    res.show()
    