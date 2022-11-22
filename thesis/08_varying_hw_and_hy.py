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

#%% --- apply lift factor ---

if 0:
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
    generator.update_dimensions(h_wndg_r = 0.035)
    generator.apply_lift_factor(pretty=True)
    
#%% --- K_s and P depending on h_wndg_s without lift factor ---

if 0:
    # --- init generator with main parameter ---
    p = 30
    l_e = 0.267
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
    ax1.set_xlabel("h_wndg_s / m")
    ax1.set_ylabel("K_s / kA/m")
    ax2.set_ylabel("P / MW")
    # ax.set_xlim(0, 3)
    # ax.set_ylim(0, 25)
    # ax.set_xticks([i for i in range(0, 181, 20)])
    # ax.set_yticks([i for i in range(0, 2501, 250)])
    plt.grid()
    
    ax1.plot(history_h_wndg_s, [K*1e-3 for K in history_K_s], color = "blue", label=f"K")
    ax2.plot(history_h_wndg_s, [P*1e-6 for P in history_P], color = "red", label=f"P")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="lower right")
    
    # generator.fast_flux()
    pass
    
#%% --- K_s and P depending on h_wndg_s with lift factor ---

if 0:
    # --- init generator with main parameter ---
    p = 30
    l_e = 0.267
    r_so = sw.r_so(sw, l_e)
    
    generator = n7_Model(p, l_e, r_so, 
                         gn=gn, fw=fw, sw=sw,
                         k_fill_s=0.2,
                         k_fill_r=0.2,
                         B_yoke_max=1.7)
    
    
    # --- init generator dimension via yoke and winding heights ---
    h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
    h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.2, 1.2]]
    
    generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
        
    # --- compute initial current loadings and with it the initial model ---
    K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
    K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
    
    generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_p=0.866)
    generator.apply_lift_factor(verbose=False)
    
    # --- adapt stator yoke ---
    h_yoke_s = yoke_height(generator.p, generator.dims.r_si, 
                            generator.B_s, generator.B_yoke_max)    
    # --- adapt rotor yoke ---
    h_yoke_r = yoke_height(generator.p, generator.dims.r_ro, 
                            generator.B_r, generator.B_yoke_max)

    generator.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r)#, keep_const_kr_b=True)
    generator.apply_lift_factor(verbose=False)
    
    # generator.fast_flux()

    h_wndg_s_0 = 0
    h_wndg_s_step = 0.01
    history_K_s = []
    history_h_wndg_s = []
    history_P = []
    
    for i in range(50):
        h_wndg_s = h_wndg_s_0 + h_wndg_s_step * i
        
        generator.update_dimensions(h_wndg_s=h_wndg_s)
        generator.apply_lift_factor()
        # K_s = generator.k_fill_s * generator.J_e_s * h_wndg_s
        # generator.update_model_by_K(K_s = K_s, K_r = K_r)
        
        history_h_wndg_s.append(h_wndg_s)
        history_K_s.append(generator.K_s)
        history_P.append(generator.P)
        
        
    fig = plt.figure(dpi=300, figsize=(6,5))
    ax1 = plt.subplot()
    ax2 = ax1.twinx()
    ax1.set_xlabel("h_wndg_s / m")
    ax1.set_ylabel("K_s / kA/m")
    ax2.set_ylabel("P / MW")
    # ax.set_xlim(0, 3)
    # ax.set_ylim(0, 25)
    # ax.set_xticks([i for i in range(0, 181, 20)])
    # ax.set_yticks([i for i in range(0, 2501, 250)])
    plt.grid()
    
    ax1.plot(history_h_wndg_s, [K*1e-3 for K in history_K_s], color = "blue", label="K")
    ax2.plot(history_h_wndg_s, [P*1e-6 for P in history_P], color = "red", label="P")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="lower right")
    
    # generator.fast_flux()
    pass
    
#%% ---  ---

if 0:
    # --- init generator with main parameter ---
    p = 30
    l_e = 0.267
    r_so = sw.r_so(sw, l_e)
    
    generator = n7_Model(p, l_e, r_so, 
                         gn=gn, fw=fw, sw=sw,
                         k_fill_s=0.3,
                         k_fill_r=0.3,
                         B_yoke_max=1.7)
    
    
    # --- init generator dimension via yoke and winding heights ---
    h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
    h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.2, 1.2]]
    
    generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
        
    # --- compute initial current loadings and with it the initial model ---
    K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
    K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
    
    generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_p=0.866)
    generator.apply_lift_factor(verbose=False)
    
    # --- adapt stator yoke ---
    h_yoke_s = yoke_height(generator.p, generator.dims.r_si, 
                            generator.B_s, generator.B_yoke_max)    
    # --- adapt rotor yoke ---
    h_yoke_r = yoke_height(generator.p, generator.dims.r_ro, 
                            generator.B_r, generator.B_yoke_max)

    generator.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r)#, keep_const_kr_b=True)
    generator.apply_lift_factor(verbose=False)
    
    generator.fast_flux()
    
    h_wndg_r_0 = 0
    h_wndg_r_step = 0.005
    history_K_r = []
    history_h_wndg_r = []
    history_P = []
    
    for i in range(40):
        h_wndg_r = h_wndg_r_0 + h_wndg_r_step * i
        
        generator.update_dimensions(h_wndg_r=h_wndg_r)
        generator.apply_lift_factor()
        # K_s = generator.k_fill_s * generator.J_e_s * h_wndg_s
        # generator.update_model_by_K(K_s = K_s, K_r = K_r)
        
        history_h_wndg_r.append(h_wndg_r)
        history_K_r.append(generator.K_r)
        history_P.append(generator.P)
        
        
    fig = plt.figure(dpi=300, figsize=(6,5))
    ax1 = plt.subplot()
    ax2 = ax1.twinx()
    ax1.set_xlabel("h_wndg_r / m")
    ax1.set_ylabel("K_r / kA/m")
    ax2.set_ylabel("P / MW")
    plt.grid()
    
    ax1.plot(history_h_wndg_r, [K*1e-3 for K in history_K_r], color = "blue", label="K")
    ax2.plot(history_h_wndg_r, [P*1e-6 for P in history_P], color = "red", label="P")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="lower right")
    
    generator.fast_flux()
    
    
#%%

if 1:
    p = 30
    l_e = 0.267
    r_so = sw.r_so(sw, l_e)
    
    geno = n7_Model(p, l_e, r_so, 
                    gn=gn, fw=fw, sw=sw,
                    k_fill_s=0.3,
                    k_fill_r=0.4,
                    B_yoke_max=1.7)
    
    geno.init_dimensions(h_yoke_s = 0.02 * r_so, h_yoke_r = 0.02 * r_so, 
                              h_wndg_s = 0.02, h_wndg_r = 0.02)
    
    geno.update_model_by_K(K_s = geno.k_fill_s * geno.h_wndg_s * sw.J_e,
                                K_r = geno.k_fill_r * geno.h_wndg_r * fw.J_e,
                                ks_p = 0.866)
    
    geno.apply_coil_sizes_and_lift_factor(verbose = False)
    geno.keep_const_kr_b(kr_b = 0.65)
    
    # geno.show_results()
    
    def adapt_yokes():
        # --- adapt stator yoke ---
        h_yoke_s = yoke_height(geno.p, geno.dims.r_si, geno.B_s, geno.B_yoke_max)    
        # --- adapt rotor yoke ---
        h_yoke_r = yoke_height(geno.p, geno.dims.r_ro, geno.B_r, geno.B_yoke_max)

        geno.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r, keep_const_kr_b=True)
        geno.apply_coil_sizes_and_lift_factor(verbose = False)

    adapt_yokes()
    geno.apply_coil_sizes_and_lift_factor(verbose = False)
    
    # geno.show_results()

        
    if 0:
        increasing_h_by_factor=0.005
        break_when_P_at_factor=0.5
    
        lst_P, lst_h_wndg_r, lst_Brc = [geno.P], [geno.h_wndg_r], [geno.B_r_c]
        for idx_h_wndg_r in range(30):
            h_wndg_r = lst_h_wndg_r[-1] + lst_h_wndg_r[0] * increasing_h_by_factor
            geno.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
            geno.apply_coil_sizes_and_lift_factor(verbose=False)
            adapt_yokes()
            geno.apply_coil_sizes_and_lift_factor(verbose=False)
            
            lst_P.append(geno.P)
            lst_h_wndg_r.append(h_wndg_r)
            lst_Brc.append(geno.B_r_c)
                    
            if lst_P[-1] < lst_P[0] * break_when_P_at_factor:
                break
        
        fig, ax1 = plt.subplots()
        plt.grid()
        ax2 = ax1.twinx()
        fig.tight_layout()         
        fig.dpi=500
        ax1.scatter([x for x in range(len(lst_P))], [y*1e-6 for y in lst_P], marker=".", c="b")
        ax1.plot([x for x in range(len(lst_P))], [y*1e-6 for y in lst_P], c="b")
        ax1.set_ylabel('P [MW]', color="b")
        ax2.scatter([x for x in range(len(lst_h_wndg_r))], [h*1e3 for h in lst_h_wndg_r], marker=".", c="r")
        ax2.plot([x for x in range(len(lst_h_wndg_r))], [h*1e3 for h in lst_h_wndg_r], c="r")
        ax2.set_ylabel('h_wndg_r [m]', color="r")
        # ax2.plot([x for x in range(len(lst_Brc))], lst_Brc, c="r")
        # ax2.set_ylabel('B_r_c [T]', color="r")
        
        tkz.clean_figure()
        tkz.save("aus_P_and_hrw_over_iter.tex")
        
    geno.update_dimensions(h_wndg_r = 0.022, keep_const_kr_b=True)
    geno.apply_coil_sizes_and_lift_factor(verbose = False)
    geno.show_results("this")
    
    if 0:
        increasing_h_by_factor=0.08
        break_when_P_at_factor=0.5
    
        lst_P, lst_h_wndg_s, lst_Bsc = [geno.P], [geno.h_wndg_s], [geno.B_s_c]
        for idx_h_wndg_s in range(30):
            h_wndg_s = lst_h_wndg_s[-1] + lst_h_wndg_s[0] * increasing_h_by_factor
            geno.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=True)
            geno.apply_coil_sizes_and_lift_factor(verbose=False)
            adapt_yokes()
            geno.apply_coil_sizes_and_lift_factor(verbose=False)
            
            # if geno.coil.r_s_bend > geno.gn.r_bend_max: print("Valid") 
            # else:  print("Invalid")
            lst_P.append(geno.P)
            lst_h_wndg_s.append(h_wndg_s)
            lst_Bsc.append(geno.B_s_c)
                    
            if lst_P[-1] < lst_P[0] * break_when_P_at_factor:
                break
        
        fig, ax1 = plt.subplots()
        plt.grid()
        ax2 = ax1.twinx()
        fig.tight_layout()  
        fig.dpi=500
        ax1.scatter([x for x in range(len(lst_P))], [y*1e-6 for y in lst_P], marker=".", c="b")
        ax1.plot([x for x in range(len(lst_P))], [y*1e-6 for y in lst_P], c="b")
        ax1.set_ylabel('P [MW]', color="b")
        ax2.scatter([x for x in range(len(lst_h_wndg_s))], [h*1e3 for h in lst_h_wndg_s], marker=".", c="r")
        ax2.plot([x for x in range(len(lst_h_wndg_s))], [h*1e3 for h in lst_h_wndg_s], c="r")
        ax2.set_ylabel('h_wndg_s [m]', color="r")
        # ax2.plot([x for x in range(len(lst_Bsc))], lst_Bsc, c="r")
        # ax2.set_ylabel('B_r_c [T]', color="r")
        
        tkz.clean_figure()
        tkz.save("aus_P_and_hsw_over_iter.tex")
    
    geno.update_dimensions(h_wndg_s = 0.036)
    geno.apply_coil_sizes_and_lift_factor(verbose=False)
    adapt_yokes()
    geno.apply_coil_sizes_and_lift_factor()
    geno.show_results()

#%%
geno.update_dimensions(h_wndg_s = 0.036)
geno.apply_coil_sizes_and_lift_factor(verbose=False)
adapt_yokes()
geno.apply_coil_sizes_and_lift_factor()
geno.show_results()