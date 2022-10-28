# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot, RadialMultiPlot
from analytics import B_from, Phi_from, taup, K_from
from data import Generator as gn
from data import FieldWinding as fw
# --------- using the standard copper specs -----------------
from data import StatorWinding_Cu as sw
r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000

class main_dimensions():
    
    def __init__(self, r_so, r_si, r_sA, r_rF, r_ro, r_ri):
        self.r_so = r_so
        self.r_si = r_si
        self.r_sA = r_sA
        self.r_rF = r_rF
        self.r_ro = r_ro
        self.r_ri = r_ri
     
    def show(self, header="which iteration?"):
        print("---", header, "---")
        print(f"r_so = {np.round(self.r_so, 3)} [m]")
        print(f"r_si = {np.round(self.r_si, 3)} [m]")
        print(f"r_sA = {np.round(self.r_sA, 3)} [m]")
        print(f"r_rF = {np.round(self.r_rF, 3)} [m]")
        print(f"r_ro = {np.round(self.r_ro, 3)} [m]")
        print(f"r_ri = {np.round(self.r_ri, 3)} [m]")
        
        
class results():
    
    def __init__(self, M, P, K_r, K_s, h_wndg_r, h_wndg_s, h_yoke_r):
        self.M = np.round(M, 2)
        self.P = np.round(P, 2)
        self.K_r = np.round(K_r, 2)
        self.K_s = np.round(K_s, 2)
        self.h_wndg_r = np.round(h_wndg_r, 2)
        self.h_wndg_s = np.round(h_wndg_s, 2)                 
        self.h_yoke_r = np.round(h_yoke_r, 2)


def height_rotor_yoke(Ce, K, taup, l, ByokeMax):
    B_d_a = (np.sqrt(2) * Ce) / (pi**2 * K)
    Phi_h = (2 * taup * l * B_d_a) / pi
    Phi_y = 0.5 * Phi_h
    h_yoke_r = Phi_y / (l * ByokeMax)
    return h_yoke_r


def create_n_Layer_model(dims, p, l, Ks=False, Kr=False, **kwargs):
    
    is_stator_aircore = kwargs.get("is_stator_aircore", True)

    r_so = dims.r_so
    r_si = dims.r_si
    r_sA = dims.r_sA
    r_rF = dims.r_rF
    r_ro = dims.r_ro
    r_ri = dims.r_ri

    mdl = Model(p=p, l=l)
    mdl.add_layer(AirLayer(r=r_ri))
    mdl.add_layer(MagneticLayer(r=r_ro, 
                                mu_r=gn.mu_r_yoke))
    if Kr:    
        mdl.add_layer(CurrentLoading(K=Kr*np.sqrt(2), 
                                    r=r_rF, 
                                    alpha=0.5, 
                                    mu_r=1.
                                    ))
    else:
        Kr = None
    if Ks:
        mdl.add_layer(CurrentLoading(K=Ks*np.sqrt(2),
                                    r=r_sA,
                                    alpha=0.0,
                                    mu_r=1.
                                    ))
    else:
        Ks = None
    if is_stator_aircore:
        mdl.add_layer(AirLayer(r=r_si))
    else:
        mdl.add_layer(AirLayer(r=r_sA))
    mdl.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl.build()
    mdl.solve()
    mdl.total_torque()
    p_plt = PlanePlot(mdl)
    res = results(mdl.Mpos, mdl.Mpos*gn.w_syn, Kr, Ks, 
                  h_wndg_r=r_rF-r_ro, h_wndg_s=r_si-r_sA,
                  h_yoke_r=r_ro-r_ri)
    return mdl, p_plt, res
    

#%% ---------------- d2so_L ----------------------------
if False:
    r_so = 2 # [m]
    l_e  = sw.d2so_L/(2 * r_so)**2
else:
    l_e = 0.3 # [m]
    r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
    
#%% -------------- 1. assumptions -----------------
p = 26
k_fill_r, k_fill_s = 0.7, 0.5
r_si = 0.98 * r_so # [m]
B_yoke_max = 2.5 # [Tesla]                                <----- B_yoke_max

# --------------1. derivations -----------------
pole_pitch = taup(d_si=2 * r_si, p=p)


# ----------------------------------------------------------- iterations ---
models_counter = 0
iters_stator = 0
solutions, solving_dims = [], []
""" ------------ init stator iteration ------------- """
K_s0 = 4e5
K_s1 = K_s0
while True:
# --- stator iteration ---
    iters_stator += 1
    K_s = K_s0 * 0.3 + K_s1 * 0.7
    h_wndg_s = K_s / (k_fill_s * sw.J_e)
    # ---
    
    """ ------------ init rotor iteration ------------- """
    K_r0 = 5e5
    K_r1 = K_r0
    data_rP = []
    while True:
    # --- rotor iteration ---
        K_r = K_r0 * 0.3 + K_r1 * 0.7
        h_wndg_r = K_r / (k_fill_r * fw.J_e)
        # ---
        
        # ------ assume h_yoke_r with design formulas and compute a model ------
        h_yoke_r = height_rotor_yoke(sw.C_e, K_r, pole_pitch, l_e, ByokeMax=B_yoke_max)
        dims = main_dimensions(r_so, r_si, 
                               r_sA= r_si - h_wndg_s, 
                               r_rF= r_si - h_wndg_s - gn.delta_mag - h_wndg_r*0.5, 
                               r_ro= r_si - h_wndg_s - gn.delta_mag - h_wndg_r, 
                               r_ri= r_si - h_wndg_s - gn.delta_mag - h_wndg_r - h_yoke_r)
        mdl, plt, res = create_n_Layer_model(dims, p=p, l=l_e, Ks=K_s, Kr=K_r, 
                                             is_stator_aircore=False)
        # ------ get B data from the 1st model and recompute h_yoke_r -----
        mdl_B_data = mdl.get_B_data(r=np.array([dims.r_rF]), t=np.linspace(0,pole_pitch, 100))
        new_B_d_a = np.max(np.abs(mdl_B_data.Br))
        h_yoke_r = (r_si * new_B_d_a) / (p * B_yoke_max)
        dims = main_dimensions(r_so, r_si, 
                               r_sA= r_si - h_wndg_s, 
                               r_rF= r_si - h_wndg_s - gn.delta_mag - h_wndg_r*0.5, 
                               r_ro= r_si - h_wndg_s - gn.delta_mag - h_wndg_r, 
                               r_ri= r_si - h_wndg_s - gn.delta_mag - h_wndg_r - h_yoke_r)    
        mdl, plt, res = create_n_Layer_model(dims, p=p, l=l_e, Ks=K_s, Kr=K_r,
                                             is_stator_aircore=False)
        # ---------------------------------------------------------------------
        
        
        models_counter += 2
        """ --- break from rotor iteration --- """
        data_rP.append(res.P)
        if data_rP[-1] != max(data_rP):
            break
        if res.P <= 7 *1e6:
            # --- up the rotor current loading ---
            K_r0 = K_r
            K_r1 = 1.05 * K_r
            # ---
        else:
            solutions.append({"Pout": res.P, "K_s": K_s, "K_r": K_r})
            solving_dims.append(dims)
            break       
    # --- end rotor iteration ---
    
    
    """ --- break from stator iteration --- """
    # --- up the stator current loading ---
    K_s0 = K_s
    K_s1 = 1.05 * K_s
    # ---
    if len(solutions) >= 1:
        if solutions[-1]["K_r"] == K_r0:
            print(f"break due to minimum {K_r0 = } reached")
            break
        if len(solutions) > 10:
            print(f"break due to {len(solutions) = }")
            break
    if h_wndg_s > 0.2 * r_si:
        print(f"break due to {h_wndg_s = }")
        break
    if models_counter > 5000:
        print(f"break due to {models_counter = }")
        break
    
# --- end stator iteration ---

solutions

#%%
"""
for idx in range(0,2):
    mdl, plt, res = create_n_Layer_model(dims=solving_dims[idx],
                                         p=p, l=l_e,
                                         Ks=solutions[idx]["K_s"], 
                                         Kr=solutions[idx]["K_r"],
                                         is_stator_aircore=False)
    plt.fluxplot(dr, dt, lvls=10)
    # plt.quiver(dr=20, dt=200, scale=180, width=0.001)
"""
