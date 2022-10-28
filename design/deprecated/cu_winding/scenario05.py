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


def create_n_Layer_model(dims, Ks=True, Kr=True, **kwargs):
    
    p = kwargs.get("p")
    l = kwargs.get("l")
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
    l_e = 0.8 # [m]
    r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
    
#%% -------------- 1. assumptions -----------------
p = 12
k_fill_r, k_fill_s = 0.5, 0.5
r_si = 0.98 * r_so # [m]
B_yoke_max = 2.5 # [Tesla]                                <----- B_yoke_max

# --------------1. derivations -----------------
pole_pitch = taup(d_si=2 * r_si, p=p)


# ------------ iterate -------------
iter_s=0
K_s0 = 2e5
K_s1 = K_s0
data_P = []
data_K_s = []


while True:
    K_s = K_s0 * 0.3 + K_s1 * 0.7
    h_wndg_s = K_s / (k_fill_s * sw.J_e)
    
    K_r = 1e6
    h_wndg_r = K_r / (k_fill_r * fw.J_e)
    h_yoke_r = height_rotor_yoke(sw.C_e, K_r, pole_pitch, l_e, ByokeMax=B_yoke_max)
    
    dims = main_dimensions(r_so, r_si, 
                           r_sA= r_si - h_wndg_s, 
                           r_rF= r_si - h_wndg_s - gn.delta_mag - h_wndg_r*0.5, 
                           r_ro= r_si - h_wndg_s - gn.delta_mag - h_wndg_r, 
                           r_ri= r_si - h_wndg_s - gn.delta_mag - h_wndg_r - h_yoke_r)
    
    mdl, plt, res = create_n_Layer_model(dims, Ks=K_s, Kr=K_r, p=p, l=l_e,
                                         is_stator_aircore=False)
    
    mdl_B_data = mdl.get_B_data(r=np.array([dims.r_rF]), t=np.linspace(0,pole_pitch, 100))
    new_B_d_a = np.max(np.abs(mdl_B_data.Br))
    h_yoke_r = (r_si * new_B_d_a) / (p * B_yoke_max)
    
    dims = main_dimensions(r_so, r_si, 
                           r_sA= r_si - h_wndg_s, 
                           r_rF= r_si - h_wndg_s - gn.delta_mag - h_wndg_r*0.5, 
                           r_ro= r_si - h_wndg_s - gn.delta_mag - h_wndg_r, 
                           r_ri= r_si - h_wndg_s - gn.delta_mag - h_wndg_r - h_yoke_r)
    
    mdl, plt, res = create_n_Layer_model(dims, Ks=K_s, Kr=K_r, p=p, l=l_e,
                                         is_stator_aircore=False)
    
    data_P.append(res.P)
    data_K_s.append(res.K_s)
    print(f"{res.P = }, {res.K_s = }")
    
    if res.P <= 7 *1e6:
        iter_s+=1
        K_s0 = K_s1
        K_s1 *= 1.1
    else:
        break
    if data_P[-1] != max(data_P):
        break
    

#%%
plt.fluxplot(dr, dt, lvls=10)
plt.quiver(dr=20, dt=200, scale=180, width=0.001)




























