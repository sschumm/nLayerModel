# -*- coding: utf-8 -*-
#%% imports and inits
import numpy as np
import matplotlib.pyplot as pyplt
from analytics import taup, Phi_from
from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from design import Main_Params, Main_Dims, create_n_Layer_model
from design import get_L_TPL2100 as L

r_so, r_si, r_sA, r_rF, r_ro, r_ri = 0,0,0,0,0,0

def update_dimensions(h_yoke_s, h_wndg_s, h_wndg_r, h_yoke_r):
    global r_so, r_si, r_sA, r_rF, r_ro, r_ri
    r_si = r_so - h_yoke_s
    r_sA = r_si - 0.5 * h_wndg_s
    r_rF = r_sA - 0.5 * h_wndg_s - gn.delta_mag - 0.5 * h_wndg_r
    r_ro = r_rF - 0.5 * h_wndg_r
    r_ri = r_ro - h_yoke_r
    
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
update_dimensions(h_yoke_s, h_wndg_s, h_wndg_r, h_yoke_r)
params = Main_Params(p, l_e, r_so, r_si, k_fill_r, k_fill_s, B_yoke_max)
dims = Main_Dims(r_so, r_si, r_sA, r_rF, r_ro, r_ri)

# --- derive initial current loadings ---
K_s = params.k_fill_s * sw.J_e * h_wndg_s
K_r = params.k_fill_r * fw.J_e * h_wndg_r

# --- build initial model ---
mdl, plt, res = create_n_Layer_model(dims, params.p, params.l_e, Ks=K_s, Kr=K_r)


#%% comparing flux computations
flux_data = mdl.get_B_data(r=np.array([r_sA]), t=np.linspace(0,np.pi/params.p, 400))

flux_using_integration = np.trapz(np.abs(flux_data.Br.flatten())*params.l_e, np.linspace(0,taup(2*r_sA, params.p), 400))
flux_using_estimation  = Phi_from(np.max(np.abs(flux_data.Br)), taup(2*r_sA, params.p), params.l_e)

print(f"using integration: B = {flux_using_integration} [T]")
print(f"using estimation:  B = {flux_using_estimation} [T]")

# pyplt.figure(dpi=1000)
# pyplt.plot(np.linspace(0,taup(2*r_sA, params.p), 400), np.abs(flux_data.Br.flatten()))

# plt.fluxplot(1000, 1000, lvls=10)
# plt.quiver(dr=20, dt=300, scale=250, width=0.002)