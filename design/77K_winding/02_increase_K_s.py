# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from scipy.constants import pi
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot#, RadialMultiPlot
from analytics import taup, K_from, Phi_from#, B_from, 
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, Main_Results, create_n_Layer_model


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
all_dims.append(init_dims)
mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                     p   =init_params.p, 
                                     l   =init_params.l_e, 
                                     Ks  =K_s,
                                     Kr  =K_r)
all_models.append((mdl, plt, res))


""" --- 1. adaption: rotor current loading using analytics --- """
mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s)
flux_data = mdl.get_B_data(r=np.array([all_dims[-1].r_rF]), 
                              t=np.linspace(0,init_params.pole_pitch, 400))
flux_at_fw = np.max(np.abs(flux_data.Br))
K_r = K_from(flux_at_fw, C_e=sw.C_e)

h_wndg_r = K_r / (init_params.k_fill_r * fw.J_e)

update_dimensions()

""" --- 1. model --- """
all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s,
                                      Kr  =K_r)
all_models.append((mdl, plt, res))


""" 2. adaption: adjust rotor yoke --- """
flux_data = mdl.get_B_data(r=np.array([all_dims[-1].r_rF]), 
                           t=np.linspace(0,init_params.pole_pitch, 400))
flux_at_fw = np.max(np.abs(flux_data.Br))
flux_p_pole = Phi_from(flux_at_fw, 
                       taup=taup(2*all_dims[-1].r_si, init_params.p), 
                       l_e=init_params.l_e)

h_yoke_r = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)

update_dimensions()
""" --- 2. model --- """
all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s,
                                      Kr  =K_r)
all_models.append((mdl, plt, res))


""" 3. adaption: adjust stator yoke --- """
flux_data = mdl.get_B_data(r=np.array([all_dims[-1].r_sA]), 
                           t=np.linspace(0,taup(2*all_dims[-1].r_si, init_params.p), 400))
flux_at_aw = np.max(np.abs(flux_data.Br))
flux_p_pole = Phi_from(flux_at_aw, 
                       taup=taup(2*all_dims[-1].r_si, init_params.p), 
                       l_e=init_params.l_e)

h_yoke_s = (0.5 * flux_p_pole) / (init_params.l_e * init_params.B_yoke_max)

update_dimensions()
""" --- 3. model --- """
all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s,
                                      Kr  =K_r)
all_models.append((mdl, plt, res))


""" 4. adaption: increase stator current loading --- """
for i in range(10):
    K_s *= 1.05
    h_wndg_s = K_s / (init_params.k_fill_s * sw.J_e)
    
    update_dimensions()  
    """ --- ith. model --- """
    all_dims.append(Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri))
    mdl, plt, res = create_n_Layer_model(dims=all_dims[-1], 
                                          p   =init_params.p, 
                                          l   =init_params.l_e, 
                                          Ks  =K_s,
                                          Kr  =K_r)
    all_models.append((mdl, plt, res)) 


# plt.fluxplot(dr, dt, lvls=10)
# plt.quiver(dr=20, dt=200, scale=180, width=0.001)
# res.show()


for idx, mdls in enumerate(all_models):
    # all_dims[idx].show(f"iteration {idx}")
    mdls[2].show(f"iteration {idx}")
    
#%%
# for mdls in all_models:
#     mdls[1].fluxplot(dr, dt, lvls=10)

#%%
all_models[-1][1].fluxplot(dr, dt, lvls=10)
