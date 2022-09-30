# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model
from design import get_L_TPL2100 as L

plot_nth_model = False
if not plot_nth_model:
    print("not plotting models...")

Pel_goal = gn.Pel_out * 0.2
A_tape = (12*0.22) * 1e-6 # [m**2]
Ic_0 = 300 # [A]
r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000
all_dims, all_models = [], []
r_si, r_sA, r_rF, r_ro, r_ri = 0,0,0,0,0

p = 12
l_e = 0.8 # [m]
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
init_params = Main_Params(p, l_e, r_so, r_si, k_fill_r=0.5, k_fill_s=0.5, B_yoke_max=12)
init_dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

Bc_s_init = 1.2
Bc_r_init = 1.7

data_Bc_s, data_Bc_r = [Bc_s_init, Bc_s_init], [Bc_r_init, Bc_r_init]


""" --- 1. adaption: vary K_r to get an initial idea --- """
if False:
    for i in np.arange(1, 0.3, -0.1):
        K_ri = K_r * i
        h_wndg_r = K_ri / (init_params.k_fill_r * J_e_r)
        update_dimensions()
           
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_ri)
        
        # --- adjust rotor yoke ---
        flux_p_pole = Phi_from(res.B_r, 
                               taup=taup(2*r_si, init_params.p), 
                               l_e=init_params.l_e)
        h_yoke_r = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
        
        # --- adjust stator yoke ---        
        flux_p_pole = Phi_from(res.B_s, 
                                taup=taup(2*r_si, init_params.p), 
                                l_e=init_params.l_e)
        h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
        
        update_dimensions()
        """ --- decrease K_r model --- """
        all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
        mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_ri)
        all_models.append((mdl, plt, res))
        
        """ print model specs """
        for idx, mdls in enumerate(all_models):
            # all_dims[idx].show(f"iteration {idx}")
            # mdls[2].show_B_r(f"iteration {idx} - B_r =")
            # mdls[2].show_B_s(f"iteration {idx} - B_s =")
            mdls[2].show(f"iteration {idx}")

#%%
""" --- 2. --- """
# if True:
""" --- init iteration values --- """
J_e_s = L(sw.T_HTS, 0.5*data_Bc_s[-2] + 0.5*data_Bc_s[-1], 30) * Ic_0/(A_tape*1.7)
J_e_r = L(fw.T_HTS, 0.5*data_Bc_r[-2] + 0.5*data_Bc_r[-1], 30) * Ic_0/(A_tape*1.7)
print(f"{J_e_s/1e6 = }, {J_e_r/1e6 = }")

K_s = init_params.k_fill_s * J_e_s * h_wndg_s
K_r = init_params.k_fill_r * J_e_r * h_wndg_r

#%%
""" --- adapt K_r and h_wndg_r --- """
counter = 0
while True:
    counter += 1
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
   

    if res.P < 1.1 * Pel_goal and res.P > 0.9 * Pel_goal:
        print(f"{counter = }")
        res.show()
        break
    elif res.P < Pel_goal:
        K_s *= 1.1
        h_wndg_s = K_s / (init_params.k_fill_s * J_e_s)
        update_dimensions()
    else:
        K_s *= 0.9
        h_wndg_s = K_s / (init_params.k_fill_s * J_e_s)
        update_dimensions()
        

    if res.P < 1.1 * Pel_goal and res.P > 0.9 * Pel_goal:
        print(f"{counter = }")
        res.show()
        break
    elif res.P < Pel_goal:
        K_r *= 1.1
        h_wndg_r = K_r / (init_params.k_fill_r * J_e_r)
        update_dimensions()
    else:
        K_r *= 0.9
        h_wndg_r = K_r / (init_params.k_fill_r * J_e_r)
        update_dimensions()
    
    mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)

#%%
update_dimensions()
""" --- model with adapted K_r and h_wndg_r --- """   
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

# --- adjust stator yoke ---        
flux_p_pole = Phi_from(res.B_s, 
                        taup=taup(2*r_si, init_params.p), 
                        l_e=init_params.l_e)
h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)

update_dimensions()
""" --- model with adapted yokes --- """
mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                      p   =init_params.p, 
                                      l   =init_params.l_e,
                                      Ks  =K_s,
                                      Kr  =K_r)

data_Bc_r.append(res.B_r)
data_Bc_s.append(res.B_s)
res.show()



#%%
""" print model specs """
# for idx, mdls in enumerate(all_models):
#     # all_dims[idx].show(f"iteration {idx}")
#     mdls[2].show(f"iteration {idx}")
#     # mdls[2].show_B_r(f"iteration {idx} - B_r =")
#     # mdls[2].show_B_s(f"iteration {idx} - B_s =")
    
#%%
""" fluxplot models """
if plot_nth_model:    
    for mdls in all_models:
        mdls[1].fluxplot(dr, dt, lvls=10)

#%%
""" quiver models """
if plot_nth_model:
    for mdls in all_models:
        mdls[1].quiver(dr=20, dt=200, scale=250, width=0.001)
    
#%%
# all_models[-1][1].fluxplot(dr, dt, lvls=10)



