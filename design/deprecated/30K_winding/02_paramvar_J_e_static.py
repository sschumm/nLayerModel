# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from prettytable import PrettyTable
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_30K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model

K_s, K_r = 0, 0
Pel_goal = gn.Pel_out * 1.2
tolerance_lower, tolerance_upper = 0.95, 1.1
r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000
r_si, r_sA, r_rF, r_ro, r_ri = 0,0,0,0,0
h_yoke_s, h_wndg_s, h_wndg_r, h_yoke_r = 0,0,0,0

B_yoke_max, k_fill_s, k_fill_r = 20, 0.5, 0.5

def update_dimensions():
    global r_si, r_sA, r_rF, r_ro, r_ri, h_yoke_s, h_wndg_s, h_wndg_r, h_yoke_r
    r_si = r_so - h_yoke_s
    r_sA = r_si - 0.5 * h_wndg_s
    r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
    r_ro = r_rF - 0.5 * h_wndg_r
    r_ri = r_ro - h_yoke_r

""" --- iterate to create table --- """
if 1:
    
    set_of_pole_pairs = [p for p in range(10, 20,1)]
    set_of_gen_length = [0.05, 0.1, 0.15, 0.2, 0.25]
    results_table = PrettyTable(["p | l_e ->"] + set_of_gen_length)
    
    for idx_p, p in enumerate(set_of_pole_pairs):
        print(f"iteration {idx_p} of {len(set_of_pole_pairs)-1}")
        table_row = [p]
        for idx_l_e, l_e in enumerate(set_of_gen_length):
            r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
            
            h_yoke_s = 0.02 * r_so
            h_wndg_s = gn.delta_mag
            h_wndg_r = gn.delta_mag
            h_yoke_r = h_yoke_s
              
            update_dimensions()
                    
            """ --- initial params --- """
            init_params = Main_Params(p, l_e, r_so, r_si, k_fill_r=k_fill_r, k_fill_s=k_fill_s, B_yoke_max=B_yoke_max)
            init_dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)
            
            J_e_s_init = sw.J_e * 1 # reduce given vals by factor
            J_e_r_init = fw.J_e * 1 # reduce given vals by factor 
            
            K_s_init = init_params.k_fill_s * h_wndg_s * J_e_s_init
            K_r_init = init_params.k_fill_r * h_wndg_r * J_e_r_init
            
            """ --- find model with static K_s and K_r --- """
            K_s = K_s_init
            K_r = K_r_init
            iter_static_P_history = []
            solution = False
            for iter_idx_static in range(50):
                
                mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                                      p   =init_params.p, 
                                                      l   =init_params.l_e,
                                                      Ks  =K_s,
                                                      Kr  =K_r)
                iter_static_P_history.append(res.P)
                
                # --- adapt rotor current loading ---
                if res.P < tolerance_upper * Pel_goal and res.P > tolerance_lower * Pel_goal:
                    # res.show(f"{iter_idx_static = }")
                    solution = True
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
                
                # --- adapt stator current loading ---
                if res.P < tolerance_upper * Pel_goal and res.P > tolerance_lower * Pel_goal:
                    # res.show(f"{iter_idx_static = }")                    
                    solution = True
                    break
                elif res.P < Pel_goal:
                    K_s *= 1.05
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
                
                # # --- adjust stator yoke ---        
                # flux_p_pole = Phi_from(res.B_s, 
                #                         taup=taup(2*r_si, init_params.p), 
                #                         l_e=init_params.l_e)
                # h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
                # update_dimensions()
                
                # mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                #                                       p   =init_params.p, 
                #                                       l   =init_params.l_e,
                #                                       Ks  =K_s,
                #                                       Kr  =K_r)
                # iter_static_P_history.append(res.P)
                
                # --- detect diverging loadings ---
                if iter_idx_static > 20 and iter_static_P_history[-1] < Pel_goal * 0.2:
                    # print(f"detected possible divergence of loadings in {iter_idx_static = }")
                    # pyplt.figure(dpi=1000)
                    # pyplt.plot([i for i in range(len(iter_static_P_history))], iter_static_P_history)
                    solution = False
                    break
            
            if solution:
                # --- adjust stator yoke ---        
                flux_p_pole = Phi_from(res.B_s, 
                                        taup=taup(2*r_si, init_params.p), 
                                        l_e=init_params.l_e)
                h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
                update_dimensions()
                
                mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                                      p   =init_params.p, 
                                                      l   =init_params.l_e,
                                                      Ks  =K_s,
                                                      Kr  =K_r)
                
                # --- adjust rotor yoke ---
                flux_p_pole = Phi_from(res.B_r, 
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
                table_row.append(res.get_some())
            else: 
                table_row.append("diverged...")

        results_table.add_row(table_row)
    
    print(results_table)

""" --- build one model --- """
if 0:
    p = 36
    l_e = 0.2
    
    if True:
        r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
        
        h_yoke_s = 0.02 * r_so
        h_wndg_s = gn.delta_mag
        h_wndg_r = gn.delta_mag
        h_yoke_r = h_yoke_s
          
        update_dimensions()
                
        """ --- initial params --- """
        init_params = Main_Params(p, l_e, r_so, r_si, k_fill_r=k_fill_r, k_fill_s=k_fill_s, B_yoke_max=B_yoke_max)
        init_dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)
        
        J_e_s_init = sw.J_e * 1 # reduce given vals by factor
        J_e_r_init = fw.J_e * 1 # reduce given vals by factor 
        
        K_s_init = init_params.k_fill_s * h_wndg_s * J_e_s_init
        K_r_init = init_params.k_fill_r * h_wndg_r * J_e_r_init
        
        """ --- find model with static K_s and K_r --- """
        K_s = K_s_init
        K_r = K_r_init
        iter_static_P_history = []
        solution = False
        for iter_idx_static in range(50):
            
            mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                                  p   =init_params.p, 
                                                  l   =init_params.l_e,
                                                  Ks  =K_s,
                                                  Kr  =K_r)
            iter_static_P_history.append(res.P)
            
            # --- adapt stator current loading ---
            if res.P < tolerance_upper * Pel_goal and res.P > tolerance_lower * Pel_goal:
                # res.show(f"{iter_idx_static = }")                    
                solution = True
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
            flux_p_pole = Phi_from(res.B_s, 
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
                # res.show(f"{iter_idx_static = }")
                solution = True
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
                solution = False
                break
        
        if solution:
            # --- adjust stator yoke ---        
            flux_p_pole = Phi_from(res.B_s, 
                                    taup=taup(2*r_si, init_params.p), 
                                    l_e=init_params.l_e)
            h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
            update_dimensions()
            
            mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                                  p   =init_params.p, 
                                                  l   =init_params.l_e,
                                                  Ks  =K_s,
                                                  Kr  =K_r)
            
            # --- adjust rotor yoke ---
            flux_p_pole = Phi_from(res.B_r, 
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
            res.show()
    
    if True: plt.fluxplot(dr, dt, lvls=10)
    if True: plt.quiver(dr=20, dt=200, scale=250, width=0.001)