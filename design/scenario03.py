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
K_s = 1.7*1e5 # [A/m]
B_yoke_max = 2.5 # [Tesla]                                <----- B_yoke_max

# --------------1. derivations -----------------
pole_pitch = taup(d_si=2 * r_si, p=p)
B_d1_a_pre = B_from(K_s, sw.C_e, kw=1)
flux_main = Phi_from(B_d1_a_pre, pole_pitch, l_e=l_e)
flux_yoke = flux_main * 0.5

h_yoke_s = r_so - r_si
h_yoke_r = flux_yoke / (B_yoke_max * l_e)
h_wndg_s = K_s / (k_fill_s * sw.J_e)
h_wndg_r = h_wndg_s

r_sA = r_si - h_wndg_s
r_rF = r_sA - gn.delta_mag - h_wndg_r * 0.5 # due to air core winding in rotor
r_ro = r_rF - h_wndg_r * 0.5
r_ri = r_ro - h_yoke_r

#%% ------------ 1. computations -----------------
# dim_iter_1 = main_dimensions(r_so=r_so, r_si=r_sA, 
#                              r_sA=r_sA, r_rF=r_rF, 
#                              r_ro=r_ro, r_ri=r_ri)
# dim_iter_1.show(header="iteration 1:")

# mdl_iter_1, plt_iter_1 = create_n_Layer_model(dim_iter_1, 
#                                               Kr=False, 
#                                               p=p, l=l_e, 
#                                               K_s_a=np.sqrt(2)*K_s)
# plt_iter_1.fluxplot(dr, dt, lvls=15)
# B_iter_1_rF = mdl_iter_1.get_B_data(r=np.array([r_rF]), 
#                                     t=np.linspace(0, (2*pi)*(pole_pitch/(2*r_rF*pi)), 40)
#                                     )
# B_iter_1_d = mdl_iter_1.get_B_data(r=np.array([r_rF + gn.delta_mag*0.5]), 
#                                    t=np.linspace(0, (2*pi)*(pole_pitch/(2*r_rF*pi)), 40)
#                                    )
# B_d1_a_post = np.max(np.sqrt(B_iter_1_d.Br**2 + B_iter_1_d.Bt**2)) * np.sqrt(2)

# results of iteration 1:
#     - big rotor yoke
#     - flux density in airgap: radial (0.0, 0.51) [T], tangential (0.0, 0.25) [T]
#       -> big airgap leads to flux being not mainly radial
#     - Flux density in Airgap much lower 

#%% --------- fit rotor current loading ------------
cond, stop = True, 0
K_r0 = K_s
K_r1 = K_s
mdl_iter_2, plt_iter_2, dim_iter_2 = None, None, None
while cond:
    K_r = np.round(0.3 * K_r0 + 0.7 * K_r1, 2)
    
    h_wndg_r = K_r / (k_fill_r * fw.J_e)
    r_rF = r_sA - gn.delta_mag - h_wndg_r * 0.5
    r_ro = r_rF - h_wndg_r * 0.5
    r_ri = r_ro - h_yoke_r
    
    dim_iter_2 = main_dimensions(r_so=r_so, r_si=r_sA, 
                                 r_sA=r_sA, r_rF=r_rF, 
                                 r_ro=r_ro, r_ri=r_ri)
    # dim_iter_2.show(header=f"iteration {stop}:")
    
    mdl_iter_2, plt_iter_2 = create_n_Layer_model(dim_iter_2, 
                                                  Ks=False,
                                                  Kr=True,
                                                  K_s_a=None,
                                                  K_r_a=np.sqrt(2)*K_r,
                                                  p=p, l=l_e)
    # print(f"M = {mdl_iter_2.Mpos/1e6} [MNm]")
    # print(f"P = {gn.w_syn * mdl_iter_2.Mneg/1e6} [MW]")
    
    # B_iter_2_rF = mdl_iter_2.get_B_data(r=np.array([r_rF]), 
    #                                     t=np.linspace(0, (2*pi)*(pole_pitch/(2*r_rF*pi)), 40)
    #                                     )
    B_iter_2_d = mdl_iter_2.get_B_data(r=np.array([r_rF + h_wndg_r*0.5 + gn.delta_mag*0.5]), 
                                       t=np.linspace(0, (2*pi)*(pole_pitch/(2*r_rF*pi)), 40)
                                       )
    B_d2_a_post = np.max(np.sqrt(B_iter_2_d.Br**2 + B_iter_2_d.Bt**2)) * np.sqrt(2)
    
    if np.isclose(B_d1_a_pre, B_d2_a_post, atol=0.01):
        cond = False
        print(f"Found Solution with {K_r = }{(12 - len(str(K_r))) * ' '} and   {B_d2_a_post = }")
    elif stop >= 100:
        cond = False
        print("Aborting due to stop condition")
    else:
        stop += 1
        print(f"iteration: {stop}{(3 - len(str(stop))) * ' '} with {K_r = }{(12 - len(str(K_r))) * ' '} and   {B_d2_a_post = }")
        if B_d2_a_post < B_d1_a_pre:
            K_r0 = K_r1
            K_r1 = K_r1 * 1.1
        else:
            K_r0 = K_r1
            K_r1 = K_r1 * 0.9
    
#%%
dim_iter_2.show(header=f"iteration {stop}:")
# plt_iter_2.fluxplot(dr, dt, lvls=15)
rad_iter_2 = RadialMultiPlot(mdl_iter_2)
# rad_iter_2.multiplot()

#%%

mdl_iter_3, plt_iter_3 = create_n_Layer_model(dim_iter_2, 
                                              Ks=True,
                                              Kr=True,
                                              K_s_a=np.sqrt(2)*K_s,
                                              K_r_a=np.sqrt(2)*K_r,
                                              p=p, l=l_e)

plt_iter_3.fluxplot(dr, dt, lvls=15)
plt_iter_3.quiver(dr=20, dt=200, scale=100, width=0.001)
print(f"M = {mdl_iter_3.Mpos/1e6} [MNm]")
print(f"P = {gn.w_syn * mdl_iter_3.Mneg/1e6} [MW]")


#%% --------- fit stator current loading ------------
B_d4_a_pre = B_from(K_r, sw.C_e, kw=1)


mdl_iter_4, plt_iter_4, dim_iter_4 = None, None, None

    
mdl_iter_4, plt_iter_4 = create_n_Layer_model(dim_iter_2, 
                                              Ks=True,
                                              Kr=False,
                                              K_s_a=np.sqrt(2)*K_s,
                                              K_r_a=None,
                                              p=p, l=l_e)

B_iter_4_d = mdl_iter_4.get_B_data(r=np.array([r_sA - gn.delta_mag*0.5]), 
                                   t=np.linspace(0, (2*pi)*(pole_pitch/(2*r_sA*pi)), 80)
                                   )
B_d4_a_post = np.max(np.sqrt(B_iter_4_d.Br**2 + B_iter_4_d.Bt**2)) * np.sqrt(2)

# plt_iter_4.fluxplot(dr, dt, lvls=15)
# plt_iter_4.quiver(dr=20, dt=200, scale=100, width=0.001)



