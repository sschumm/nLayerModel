# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_30K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model

plot_nth_model = False

r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000
all_dims, all_models = [], []
r_si, r_sA, r_rF, r_ro, r_ri = 0,0,0,0,0

p = 14
l_e = 0.05 # [m]
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

K_s = init_params.k_fill_s * sw.J_e * h_wndg_s
K_r = init_params.k_fill_r * fw.J_e * h_wndg_r

""" --- initial model --- """
if plot_nth_model: 
    all_dims.append(init_dims)
    mdl, plt, res = create_n_Layer_model(dims=init_dims,
                                         p   =init_params.p, 
                                         l   =init_params.l_e, 
                                         Ks  =K_s,
                                         Kr  =K_r)
    all_models.append((mdl, plt, res))

#%% --- 1. adaption: varyK_s to get an initial idea --- 
if 0:
    h_wndg_r = gn.delta_mag
    K_r = init_params.k_fill_r * fw.J_e * h_wndg_r
    for i in [i/10 for i in range(3,20, 4)]:
        K_si = K_s * i
        h_wndg_s = K_si / (init_params.k_fill_s * sw.J_e)
        update_dimensions()
           
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e, 
                                              Ks  =K_si,
                                              Kr  =K_r)
        res.show_P(f"{i = }")
        # plt.fluxplot(dr, dt, lvls=10)

#%% --- 2. adaption: vary K_r to get an initial idea ---
if 0:
    h_wndg_s = gn.delta_mag
    K_s = init_params.k_fill_s * sw.J_e * h_wndg_s
    for i in [i/10 for i in range(3,20, 4)]:
        K_ri = K_r * i
        h_wndg_r = K_ri / (init_params.k_fill_r * fw.J_e)
        update_dimensions()
           
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e, 
                                              Ks  =K_s,
                                              Kr  =K_ri)
        res.show_P(f"{i = }")
        # plt.fluxplot(dr, dt, lvls=10) 
        
#%% """
if 1:
    h_wndg_s = gn.delta_mag 
    h_wndg_r = gn.delta_mag 
    h_yoke_s = 0.02 * r_so *2
    h_yoke_r = h_yoke_s 
    K_s = init_params.k_fill_s * sw.J_e * h_wndg_s
    K_r = init_params.k_fill_r * fw.J_e * h_wndg_r

    update_dimensions()
       
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e, 
                                          Ks  =K_s,
                                          Kr  =K_r)
    res.show_P()
    res.show_B_airgap()
    plt.fluxplot(dr, dt, lvls=10)

# #%%
# """ print model specs """
# for idx, mdls in enumerate(all_models):
#     # all_dims[idx].show(f"iteration {idx}")
#     mdls[2].show_P(f"iteration {idx}")
    
# #%%
# """ fluxplot models """
# if plot_nth_model:    
#     for mdls in all_models:
#         mdls[1].fluxplot(dr, dt, lvls=10)

# #%%
# """ quiver models """
# if plot_nth_model:
#     for mdls in all_models:
#         mdls[1].quiver(dr=20, dt=200, scale=250, width=0.001)
    
# #%%
# # all_models[-1][1].fluxplot(dr, dt, lvls=10)



