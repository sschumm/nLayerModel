# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplt
from scipy.constants import pi
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot#, RadialMultiPlot
from analytics import taup, K_from#, B_from, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, Main_Results, create_n_Layer_model


r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000


p = 20
l_e = 0.5 # [m]
r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)

h_yoke_s = 0.02 * r_so
h_wndg_s = gn.delta_mag
h_wndg_r = gn.delta_mag
h_yoke_r = h_yoke_s

r_si = r_so - h_yoke_s
r_sA = r_si - 0.5 * h_wndg_s
r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
r_ro = r_rF - 0.5 * h_wndg_r
r_ri = r_ro - h_yoke_r

""" --- initial parameter --- """
init_params = Main_Params(p, l_e, r_so, r_si)
init_dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

K_s = init_params.k_fill_s * sw.J_e * h_wndg_s
K_r = init_params.k_fill_r * fw.J_e * h_wndg_r

mdl, plt, res = create_n_Layer_model(dims=init_dims, 
                                     p   =init_params.p, 
                                     l   =init_params.l_e, 
                                     Ks  =K_s,
                                     Kr  =K_r)
plt.fluxplot(dr, dt, lvls=10)

""" --- 1. adaption: rotor current loading --- """
mdl, plt, res = create_n_Layer_model(dims=init_dims, 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s)
flux_data_Ks = mdl.get_B_data(r=np.array([init_dims.r_rF]), 
                              t=np.linspace(0,init_params.pole_pitch, 400))
# pyplt.figure(figsize=(10,10))
# pyplt.plot(flux_data_Ks.T, flux_data_Ks.Br)
flux_at_fw = np.max(np.abs(flux_data_Ks.Br))
K_r = K_from(flux_at_fw, C_e=sw.C_e)

h_wndg_r = K_r / (init_params.k_fill_r * fw.J_e)

r_si = r_so - h_yoke_s
r_sA = r_si - 0.5 * h_wndg_s
r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
r_ro = r_rF - 0.5 * h_wndg_r
r_ri = r_ro - h_yoke_r

dims_iter1 = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

mdl, plt, res = create_n_Layer_model(dims=dims_iter1, 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s,
                                      Kr  =K_r)
plt.fluxplot(dr, dt, lvls=10)

""" --- 2. adaption: --- """
mdl, plt, res = create_n_Layer_model(dims=dims_iter1, 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Kr  =K_r)
flux_data_Kr = mdl.get_B_data(r=np.array([init_dims.r_sA]), 
                              t=np.linspace(0,init_params.pole_pitch, 400))
# pyplt.figure(figsize=(10,10))
# pyplt.plot(flux_data_Kr.T, flux_data_Kr.Br)
flux_at_sw = np.max(np.abs(flux_data_Kr.Br))
K_s = K_from(flux_at_sw, C_e=sw.C_e)

h_wndg_s = K_s / (init_params.k_fill_s * sw.J_e)

r_si = r_so - h_yoke_s
r_sA = r_si - 0.5 * h_wndg_s
r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
r_ro = r_rF - 0.5 * h_wndg_r
r_ri = r_ro - h_yoke_r

dims_iter2 = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

mdl, plt, res = create_n_Layer_model(dims=dims_iter2, 
                                      p   =init_params.p, 
                                      l   =init_params.l_e, 
                                      Ks  =K_s,
                                      Kr  =K_r)
plt.fluxplot(dr, dt, lvls=10)




# plt.fluxplot(dr, dt, lvls=10)
# plt.quiver(dr=20, dt=200, scale=180, width=0.001)
res.show()




