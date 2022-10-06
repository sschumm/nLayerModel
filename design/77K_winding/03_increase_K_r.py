# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model

plot_nth_model = True
if not plot_nth_model:
    print("not plotting models...")

r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000
all_dims, all_models = [], []
r_si, r_sA, r_rF, r_ro, r_ri = 0,0,0,0,0

p = 20
l_e = 0.5 # [m]
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
init_params = Main_Params(p, l_e, r_so, r_si, B_yoke_max=6)
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

#%%
""" --- 1. adaption: vary K_r to get an initial idea --- """
if False:
    for i in np.arange(1,0.3, -0.1):
        K_ri = K_r * i
        h_wndg_r = K_ri / (init_params.k_fill_r * fw.J_e)
        update_dimensions()
           
        mdl, plt, res = create_n_Layer_model(dims=Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri), 
                                              p   =init_params.p, 
                                              l   =init_params.l_e, 
                                              Kr  =K_ri)
        # --- adjust stator yoke ---
        flux_data = mdl.get_B_data(r=np.array([r_sA]), 
                                    t=np.linspace(0,taup(2*r_si, init_params.p), 400))
        flux_at_aw = np.max(np.abs(flux_data.Br))
        flux_p_pole = Phi_from(flux_at_aw, 
                                taup=taup(2*r_si, init_params.p), 
                                l_e=init_params.l_e)
    
        h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
        
        # --- adjust rotor yoke ---
        flux_data = mdl.get_B_data(r=np.array([r_rF]), 
                                   t=np.linspace(0,taup(2*r_si, init_params.p), 400))
        flux_at_fw = np.max(np.abs(flux_data.Br))
        flux_p_pole = Phi_from(flux_at_fw, 
                               taup=taup(2*r_si, init_params.p), 
                               l_e=init_params.l_e)
        
        h_yoke_r = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)
    
    
        update_dimensions()
        """ --- decrease K_r model --- """
        all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
        mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_ri)
        all_models.append((mdl, plt, res))

#%%
""" --- select K_r and with that h_yoke_r, h_wndg_r, h_yoke_s --- """
K_r = 1270.4e3
h_wndg_r = K_r / (init_params.k_fill_r * fw.J_e)

h_yoke_r = 0.0774
h_yoke_s = 0.0618

update_dimensions()
""" --- 2. model --- """
if plot_nth_model: 
    all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
    mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    all_models.append((mdl, plt, res))

#%%
""" --- 2. adaption: increase K_s to find tipping point --- """
if False:
    for i in np.arange(1, 5, 0.2):
        K_si = K_s * i
        h_wndg_s = K_si / (init_params.k_fill_s * sw.J_e)
        
        update_dimensions()
        """ --- increase K_s model --- """
        all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
        mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_si,
                                              Kr  =K_r)
        all_models.append((mdl, plt, res))
    pyplt.figure(dpi=1000)
    pyplt.ticklabel_format(style='plain')
    pyplt.xlabel("K_s [kA]")
    pyplt.ylabel("P [MW]")
    pyplt.scatter([m[2].K_s * 1e-3 for m in all_models],
                  [m[2].P * 1e-6 for m in all_models]) 
    
#%%    
""" --- select K_s from slightly before tipping point --- """
K_s = 600e3
h_wndg_s = K_s / (init_params.k_fill_s * sw.J_e)

update_dimensions()
""" --- 3. model --- """
if plot_nth_model: 
    all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
    mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    all_models.append((mdl, plt, res))

#%%
""" 3. adaption: increase K_r to achieve power goal """
if False:
    for i in np.arange(1.4,1.5,0.02):
        K_ri = K_r * i
        h_wndg_r = K_ri / (init_params.k_fill_r * fw.J_e)
        
        update_dimensions()
        """ --- increase K_r model --- """
        all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
        mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                              p   =init_params.p, 
                                              l   =init_params.l_e,
                                              Ks  =K_s,
                                              Kr  =K_ri)
        all_models.append((mdl, plt, res))
    pyplt.figure(dpi=1000)
    pyplt.ticklabel_format(style='plain')
    pyplt.xlabel("K_r [kA]")
    pyplt.ylabel("P [MW]")
    pyplt.scatter([m[2].K_r * 1e-3 for m in all_models],
                  [m[2].P * 1e-6 for m in all_models]) 

#%%
""" --- select K_s from slightly before tipping point --- """
K_r = 1900e3
h_wndg_r = K_r / (init_params.k_fill_r * fw.J_e)

update_dimensions()
""" --- 4. model --- """
if plot_nth_model: 
    all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
    mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    all_models.append((mdl, plt, res))

#%%
""" 4. adaption: adjust rotor yoke """
if plot_nth_model: 
    flux_data = mdl.get_B_data(r=np.array([r_rF]), 
                               t=np.linspace(0,taup(2*r_si, init_params.p), 400))
    flux_at_fw = np.max(np.abs(flux_data.Br))
    flux_p_pole = Phi_from(flux_at_fw, 
                           taup=taup(2*r_si, init_params.p), 
                           l_e=init_params.l_e)

    h_yoke_r = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)

    update_dimensions()
    """ --- 4. model --- """
    all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
    mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                          p   =init_params.p, 
                                          l   =init_params.l_e,
                                          Ks  =K_s,
                                          Kr  =K_r)
    all_models.append((mdl, plt, res))

#%%
""" print model specs """
for idx, mdls in enumerate(all_models):
    # all_dims[idx].show(f"iteration {idx}")
    mdls[2].show(f"iteration {idx}")
    
#%%
""" fluxplot models """
for mdls in all_models:
    mdls[1].fluxplot(dr, dt, lvls=10)

#%%
""" quiver models """
for mdls in all_models:
    mdls[1].quiver(dr=20, dt=200, scale=250, width=0.001)
    
#%%
# all_models[-1][1].fluxplot(dr, dt, lvls=10)



