# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot#, RadialMultiPlot
from analytics import taup#, B_from, Phi_from, K_from
from data import Generator as gn
from data import FieldWinding as fw
# --------- using the standard copper specs -----------------
from data import StatorWinding_Cu as sw
r_so, l_e, dims = None, None, None
mdl, plt, res = None, None, None
dr, dt = 1000, 1000

class outer_params():
    
    def __init__(self, p, l_e, r_so, r_si, k_fill_r=0.5, k_fill_s=0.5, B_yoke_max=2.5):
        self.p = p
        self.l_e = l_e
        self.r_so = r_so
        self.r_si = r_si
        
        self.k_fill_r = k_fill_r
        self.k_fill_s = k_fill_s
        self.B_yoke_max = B_yoke_max
        
    def show(self, header="which iteration?"):
        print("---", header, "---")
        print(f"p   = {self.p}")
        print(f"l_e = {np.round(self.l_e, 3)} [m]")
        
        

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
    
    def __init__(self, M, P, K_r, K_s, h_wdng_r, h_wdng_s, h_yoke_r):
        self.M = np.round(M, 2)
        self.P = np.round(P, 2)
        self.K_r = np.round(K_r, 2)
        self.K_s = np.round(K_s, 2)
        self.h_wdng_r = np.round(h_wdng_r, 2)
        self.h_wdng_s = np.round(h_wdng_s, 2)                 
        self.h_yoke_r = np.round(h_yoke_r, 2)


def height_rotor_yoke(Ce, K, tau_p, l, ByokeMax):
    B_d_a = (np.sqrt(2) * Ce) / (pi**2 * K)
    Phi_h = (2 * tau_p * l * B_d_a) / pi
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
                  h_wdng_r=r_rF-r_ro, h_wdng_s=r_si-r_sA,
                  h_yoke_r=r_ro-r_ri)
    return mdl, p_plt, res
    

def find_current_loadings(params):
    p = params.p
    l_e = params.l_e
    r_so = params.r_so
    r_si = params.r_si
    k_fill_r = params.k_fill_r
    k_fill_s = params.k_fill_s
    B_yoke_max = params.B_yoke_max
    
    # derivations from outside params
    pole_pitch = taup(d_si=2 * r_si, p=p)
    
    
    models_counter = 0
    iters_stator = 0
    solutions, solving_dims = [], []
    """ ------------ init stator iteration ------------- """
    K_s = 4e5
    while True:
    # --- stator iteration ---
        iters_stator += 1
        h_wndg_s = K_s / (k_fill_s * sw.J_e)
        # ---
        
        """ ------------ init rotor iteration ------------- """
        K_r = 5e5
        data_rP = []
        while True:
        # --- rotor iteration ---
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
                break # due to tipping point
            if res.P <= 7 *1e6:
                # --- increase the rotor current loading ---
                K_r *= 1.02
                # ---
            else:
                solutions.append({"Pout": np.round(res.P/1e6, 2), "K_s": np.round(K_s,2), "K_r": np.round(K_r, 2), "p": p, "l_e": np.round(l_e,2)})
                solving_dims.append(dims)
                break       
        # --- end rotor iteration ---
        
        
        """ --- break from stator iteration --- """
        # --- increase the stator current loading ---
        K_s *= 1.02
        # ---
        if len(solutions) >= 2:
            if solutions[-1]["K_r"] == solutions[-2]["K_r"]:
                print(f"INFO: break due to minimum {K_r = } reached")
                break
            if len(solutions) > 10:
                print(f"INFO: break due to {len(solutions) = }")
                break
        if h_wndg_s > 0.2 * r_si:
            print(f"INFO: break due to {h_wndg_s = }")
            break
        if models_counter > 5000:
            print(f"INFO: break due to {models_counter = }")
            break
        if len(solutions) == 1:
            print(f"INFO: found solution for {p = }")
            break
        
    # --- end stator iteration ---
    print(f"INFO: solved {models_counter} models")
    return solutions, solving_dims


# #%% ---- default simulation


# """ init params """
# # ---------------- d2so_L ----------------------------
# # r_so = 2 # [m]
# # l_e  = sw.d2so_L/(2 * r_so)**2
# l_e = 0.3 # [m]
# r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
    
# # -------------- p -- k -- r_si -- B_yoke_max -----------------
# p = 26
# r_si = 0.98 * r_so # [m]
# # k_fill_r, k_fill_s = 0.7, 0.5
# # B_yoke_max = 2.5 # [Tesla]

# """ find current loadings """
# pars = outer_params(p, l_e, r_so, r_si)
# solutions, solving_dims = find_current_loadings(pars)

# """ create plots """
# for idx in range(0,2):
#     mdl, plt, res = create_n_Layer_model(dims=solving_dims[idx],
#                                          p=p, l=l_e,
#                                          Ks=solutions[idx]["K_s"], 
#                                          Kr=solutions[idx]["K_r"],
#                                          is_stator_aircore=False)
#     # plt.fluxplot(dr, dt, lvls=10)
#     plt.quiver(dr=20, dt=200, scale=180, width=0.001)


#%% ---- simulation with pole count set

all_solutions, all_solving_dims = [], []
pole_counts = range(18, 26, 2)
lengths = np.arange(0.3, 0.7, 0.1)
for p in pole_counts:
    for l_e in lengths:
        # l_e = 0.3 # [m]
        r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
        r_si = 0.98 * r_so # [m]
        
        solutions, solving_dims = find_current_loadings(outer_params(p, l_e, r_so, r_si))
        all_solutions.append(solutions)
        all_solving_dims.append(solving_dims)
        
        # print(solutions)
    
    
    
#%%
for idx, sol in enumerate(all_solutions):
    if len(sol) != 0:
        mdl, plt, res = create_n_Layer_model(dims=all_solving_dims[idx][0], 
                                             p=all_solutions[idx][0]["p"], 
                                             l=all_solutions[idx][0]["l_e"], 
                                             Ks=all_solutions[idx][0]["K_s"], 
                                             Kr=all_solutions[idx][0]["K_r"],
                                             is_stator_aircore=False)
        plt.quiver(dr=20, dt=200, scale=180, width=0.001)
        
    
    
    
    
    
    
    
    
    