# -*- coding: utf-8 -*-
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as plt

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from data.export import save_params

from design import n7_Model
from modules import Model, AirLayer, MagneticLayer, CurrentLoading
from modules.plot import PlanePlot


#%% --- fluxplot to describe yoke height --- 

if 0: 
    n_syn = gn.n_syn
    
    p, l_e = 5, 0.3
    
    r1, r2, r3, r4, r5, r6 = 1.4, 2, 2.3, 2.8, 3, 3.6 # [m]
    mu_r = 1e5 # [-]
    
    alpha_r = 0.0*pi # [rad] =  0 [°]
    alpha_s = 0.38*pi # [rad] = 90 [°]
    
    K_r_hat = 8e5 # [A/m]
    K_s_hat = 1e6 # [A/m]
    
    
    model = Model(p, l_e)
    
    # model.add_layer(AirLayer(r=r3))
    # model.add_layer(MagneticLayer(r=r2, mu_r=mu_r))
    # model.add_layer(CurrentLoading(K=K_r_hat, r=r3, alpha=alpha_r))
    model.add_layer(CurrentLoading(K=K_s_hat, r=r4, alpha=alpha_s))
    model.add_layer(AirLayer(r=r5))
    model.add_layer(MagneticLayer(r=r6, mu_r=mu_r))
    
    model.build()
    model.solve()
    
    
    # --- fluxplot ---
    r_min, r_max, margin = r4, r6, 0.0
    
    plt_plane = PlanePlot(model, fgsz=50)
    plt_plane.fluxplot(dr=70, dt=80, lvls=9, lw=3,
                        r_min=r_min-margin, r_max=3.6+margin,
                        t_min=0, t_max=np.pi/4,
                        custom_dims=True,
                        y0=0, x0=0,
                        show_cbar=True, 
                        show_axis=True,
                        show_borders=True,
                        transparent=True,
                        padding=0,
                        pdf=False, pdf_dpi=300, 
                        svg=False, svg_dpi=300)


#%% --- without lift factor ---

if 0:
    # --- init generator with main parameter ---
    p = 24
    l_e = 0.3
    r_so = sw.r_so(sw, l_e)
    
    generator = n7_Model(p, l_e, r_so, 
                         gn=gn, fw=fw, sw=sw,
                         k_fill_s=0.3,
                         k_fill_r=0.3,
                         B_yoke_max=1.7)
    
    
    # --- init generator dimension via yoke and winding heights ---
    h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
    h_wndg_s, h_wndg_r =  0.01 * r_so, 0.02 * r_so #[factor * (gn.h_pole_frame*2) for factor in [1.2, 1.2]]
    
    generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
        
    # --- compute initial current loadings and with it the initial model ---
    K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
    K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
    
    generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_p=0.866)
    generator.show_results()
    # generator.fast_flux()
    
    h_wndg_s_0 = 0
    h_wndg_s_step = 0.01
    history_K_s = []
    history_h_wndg_s = []
    history_P = []
    
    for i in range(30):
        h_wndg_s = h_wndg_s_0 + h_wndg_s_step * i
        
        generator.update_dimensions(h_wndg_s=h_wndg_s)
        K_s = generator.k_fill_s * generator.J_e_s * h_wndg_s
        generator.update_model_by_K(K_s = K_s, K_r = K_r)
        
        history_h_wndg_s.append(h_wndg_s)
        history_K_s.append(K_s)
        history_P.append(generator.P)
        
        
    fig = plt.figure(dpi=300, figsize=(6,5))
    ax1 = plt.subplot()
    ax2 = ax1.twinx()
    # ax.set_xlabel("B / T")
    # ax.set_ylabel("J_e / Amm2")
    # ax.set_xlim(0, 3)
    # ax.set_ylim(0, 25)
    # ax.set_xticks([i for i in range(0, 181, 20)])
    # ax.set_yticks([i for i in range(0, 2501, 250)])
    plt.grid()
    
    ax1.plot(history_h_wndg_s, history_K_s) #, label=f"T = {lst_T[0]} / K")
    ax2.plot(history_h_wndg_s, history_P) #, label=f"T = {lst_T[0]} / K")
    # generator.fast_flux()
    
#%% --- apply lift factor ---

if 1:
    p = 30
    r_so = 3.75
    l_e = sw.l_e(sw, r_so) # 0.267
    
    generator = n7_Model(p, l_e, r_so, 
                         gn=gn, fw=fw, sw=sw,
                         k_fill_s=0.1,
                         k_fill_r=0.3,
                         B_yoke_max=1.7)

    # --- init generator dimension via yoke and winding heights ---
    h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
    h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.5, 1.5]]
    generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
        
    # --- compute initial current loadings and with it the initial model ---
    K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
    K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
    generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.866)
    generator.keep_const_kr_b(kr_b = 0.65)
    
    generator.apply_lift_factor(pretty=False)
    # generator.show_results()
    
    # # --- adapt stator yoke ---
    # h_yoke_s = yoke_height(generator.p, generator.dims.r_si, 
    #                        generator.B_s, generator.B_yoke_max)    
    # # --- adapt rotor yoke ---
    # h_yoke_r = yoke_height(generator.p, generator.dims.r_ro, 
    #                        generator.B_r, generator.B_yoke_max)

    # generator.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r)
    # generator.apply_lift_factor(pretty=True)
    # # generator.show_results()
    
    
    # # --- random height adaption ---
    # generator.update_dimensions(h_wndg_s = 0.035)
    # generator.apply_lift_factor(pretty=True)
#%%
    generator.update_dimensions(h_wndg_r = 0.035)
    generator.apply_lift_factor(pretty=True)
    
    
    