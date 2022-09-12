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
r_so, l_e = None, None
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
              


def create_n_Layer_model(dims, Ks=True, Kr=True, **kwargs):
    
    p = kwargs.get("p")
    l = kwargs.get("l")
    K_s_amplitude = kwargs.get("K_s_a")
    K_r_amplitude = kwargs.get("K_r_a")

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
        mdl.add_layer(CurrentLoading(K=K_r_amplitude, 
                                    r=r_rF, 
                                    alpha=0.5, 
                                    mu_r=1.
                                    ))
    if Ks:
        mdl.add_layer(CurrentLoading(K=K_s_amplitude,
                                    r=r_sA,
                                    alpha=0.0,
                                    mu_r=1.
                                    ))
    mdl.add_layer(AirLayer(r=r_si))
    mdl.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl.build()
    mdl.solve()
    mdl.total_torque()
    p_plt = PlanePlot(mdl)
    return mdl, p_plt
    

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

stop_iter_stator, iterate_stator_loading = 0, True
K_s0 = 2e5
K_s1 = K_s0

while iterate_stator_loading:
    stop_iter_rotor, iterate_rotor_loading  = 0, True
    K_r0 = 2e5
    K_r1 = K_r0
    
    K_s = np.round(0.3 * K_s0  + 0.7 * K_s1, 2)
    # B_d1_a_pre = B_from(K_s, sw.C_e, kw=1)
    # flux_main = Phi_from(B_d1_a_pre, pole_pitch, l_e=l_e)
    # flux_yoke = flux_main * 0.5
    
    h_yoke_s = r_so - r_si
    h_wndg_s = K_s / (k_fill_s * sw.J_e)

    h_yoke_r = h_yoke_s # flux_yoke / (B_yoke_max * l_e)
    
    r_sA = r_si - h_wndg_s
    # r_rF = r_sA - gn.delta_mag - h_wndg_r * 0.5 # due to air core winding in rotor
    # r_ro = r_rF - h_wndg_r * 0.5
    # r_ri = r_ro - h_yoke_r


    while iterate_rotor_loading:
        K_r = np.round(0.3 * K_r0 + 0.7 * K_r1, 2)
            
        h_wndg_r = K_r / (k_fill_r * fw.J_e)
        r_rF = r_sA - gn.delta_mag - h_wndg_r * 0.5
        r_ro = r_rF - h_wndg_r * 0.5
        r_ri = r_ro - h_yoke_r

        dim_iter = main_dimensions(r_so=r_so, r_si=r_sA, 
                                   r_sA=r_sA, r_rF=r_rF, 
                                   r_ro=r_ro, r_ri=r_ri)
        
        mdl_iter, plt_iter = create_n_Layer_model(dims=dim_iter, 
                                                  Ks = True, Kr = True,
                                                  K_s_a = np.sqrt(2) * K_s,
                                                  K_r_a = np.sqrt(2) * K_r,
                                                  p=p, l=l_e)
        
        P_iter = gn.w_syn * mdl_iter.Mpos/1e6
        
        if stop_iter_rotor >= 20:
            iterate_rotor_loading = False
            print("   Stopping Rotor iteration...")
        elif np.isclose(P_iter, gn.Pel_out/1e6, atol=0.5):
            iterate_rotor_loading = False
            print(f"   Found solution: {P_iter = }.")
            dim_iter.show(f"dimensions with {K_s=} and {K_r=}")
        else:
            print(f"   iter {stop_iter_rotor}{(3 - len(str(stop_iter_rotor))) * ' '} with {K_r = }{(12 - len(str(K_r))) * ' '} and {P_iter = }.")
            stop_iter_rotor += 1
            K_r0 = K_r1
            if P_iter < gn.Pel_out/1e6:
                K_r1 *= 1.1
            else:
                K_r1 *= 0.9
    
    if stop_iter_stator >= 10:
        iterate_stator_loading = False
        print("Stopping Stator iteration...")
    else:
        print(f"iter {stop_iter_stator}{(3 - len(str(stop_iter_rotor))) * ' '} with {K_s = }{(12 - len(str(K_s))) * ' '} and {P_iter = }.")
        stop_iter_stator += 1
        K_s0 = K_s1
        if P_iter < gn.Pel_out/1e6:
            K_s1 *= 1.1
        else:
            K_s1 *= 0.9

#%%
plt_iter.fluxplot(dr, dt, lvls=15)
plt_iter.quiver(dr=20, dt=200, scale=100, width = 0.001)





# lengths = list(np.linspace(1,5,100))
# diams = list()
# p_el_out = list()

# for length in lengths:
#     i_model, diam = iterate_model(length)
#     P_el_out = gn.w_syn * i_model.Mneg / 1e6
#     diams.append(diam)
#     p_el_out.append(P_el_out)

# import matplotlib.pyplot as plt
# from mpl_toolkits import mplot3d

# plt.figure(figsize=(15, 15))
# # c = plt.contourf(lengths, poles, p_el_out,
# #                  levels=100)
# # plt.colorbar(c)
# ax = plt.axes(projection="3d")
# ax.plot3D(lengths, diams, p_el_out)
# ax.set_xlabel("Generator Length")
# ax.set_ylabel("Generator Diameter")
# ax.set_zlabel("Power Output")

