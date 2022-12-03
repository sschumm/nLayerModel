# -*- coding: utf-8 -*-
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as plt

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_30K as sw
from data.export import save_params

from design import n7_Model
from modules import Model, AirLayer, MagneticLayer, CurrentLoading
from modules.plot import PlanePlot

    
    
#%%

if 1:
    p = 8
    l_e = 0.5
    r_so = 1 #sw.r_so(sw, l_e)
    
    geno = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=2)
    h_pf = gn.h_pole_frame
    
    geno.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
    
    # geno.alpha_r = 0.4*pi
    
    geno.update_model_by_K(K_s = geno.k_fill_s * geno.h_wndg_s * sw.J_e,
                                K_r = geno.k_fill_r * geno.h_wndg_r * fw.J_e,
                                ks_p = 0.866)
    geno.apply_coil_thickness_ratio()
    geno.apply_lift_factor()    
    geno.adapt_yokes()    
    geno.apply_coil_thickness_ratio()
    geno.apply_lift_factor()    
    geno.show_results()
    
    geno.update_dimensions(h_wndg_s = 0.025)
    geno.apply_coil_thickness_ratio()
    geno.apply_lift_factor()
    geno.adapt_yokes()    
    geno.apply_coil_thickness_ratio()
    geno.apply_lift_factor()
    


#%% rotor
  
    if 1:
        increasing_h_by_factor=0.06
        break_when_P_at_factor=0.5
    
        lst_P, lst_h_wndg_r, lst_Brc = [geno.P], [geno.h_wndg_r], [geno.B_r_c]
        for idx_h_wndg_r in range(40):
            h_wndg_r = lst_h_wndg_r[-1] + lst_h_wndg_r[0] * increasing_h_by_factor
            geno.update_dimensions(h_wndg_r = h_wndg_r)
            geno.apply_coil_thickness_ratio()
            geno.apply_lift_factor()
            geno.adapt_yokes()
            geno.apply_coil_thickness_ratio()
            geno.apply_lift_factor()
            geno.adapt_yokes()
            geno.apply_coil_thickness_ratio()
            geno.apply_lift_factor()
            geno.adapt_yokes()
            geno.apply_coil_thickness_ratio()
            geno.apply_lift_factor()
            
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
        # tkz.save("aus_P_and_hrw_over_iter1_04pi.tex")
        pass

#%% stator

geno.update_dimensions(h_wndg_r=0.038, h_wndg_s=h_pf * 2.1)
# geno.update_dimensions(h_wndg_r=0.060, h_wndg_s=0.038)
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()    
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()    
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.show_results()


#%% stator

if 1:
    increasing_h_by_factor= 0.06
    break_when_P_at_factor=0.5

    lst_P, lst_h_wndg_s, lst_Bsc = [geno.P], [geno.h_wndg_s], [geno.B_s_c]
    for idx_h_wndg_s in range(40):
        h_wndg_s = lst_h_wndg_s[-1] + lst_h_wndg_s[0] * increasing_h_by_factor
        geno.update_dimensions(h_wndg_s = h_wndg_s)
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        
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
    # tkz.save("aus_P_and_hsw_over_iter2_04pi.tex")
    pass

#%% rotor

geno.update_dimensions(h_wndg_r=0.038, h_wndg_s=0.03)
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.show_results()


#%% rotor
  
if 1:
    increasing_h_by_factor=0.022
    break_when_P_at_factor=0.5

    lst_P, lst_h_wndg_r, lst_Brc = [geno.P], [geno.h_wndg_r], [geno.B_r_c]
    for idx_h_wndg_r in range(40):
        h_wndg_r = lst_h_wndg_r[-1] + lst_h_wndg_r[0] * increasing_h_by_factor
        geno.update_dimensions(h_wndg_r = h_wndg_r)
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor() 
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor() 
        geno.adapt_yokes()
        geno.apply_coil_thickness_ratio()
        geno.apply_lift_factor() 
        
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
    # tkz.save("aus_P_and_hrw_over_iter2_04pi.tex")
    pass
    
    
#%%    
geno.update_dimensions(h_wndg_r = 0.055, h_wndg_s = 0.03)
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()  
geno.show_results()
# geno.fast_flux()

#%%

print(f"{geno.l_e = }, {geno.k_fill_s = }, {geno.k_fill_r = }")

N_s = int(np.round(geno.guess_Ns()))
N_c_div_a = int(np.round(geno.guess_N_c_div_a(N_s = N_s)))

print(f"\n{N_s = } \n{N_c_div_a = }")

I_s = geno.P / (np.sqrt(3) * geno.gn.U_LL_N)
ampere_turns_stator = N_c_div_a * I_s

J_s_amplitude = np.sqrt(2) * ampere_turns_stator / geno.coil.A_sc
geno.N_s = N_s
geno.J_s_amplitude = J_s_amplitude

print(f"\n{I_s = } [A] \n{ampere_turns_stator = } [A] \n{J_s_amplitude = } [A/m2] \n{J_s_amplitude*1e-6 = }[A_mm2]\n")

#%%
I_s = geno.P/ (np.sqrt(3) * geno.gn.U_LL_N * 0.94)
N_c_div_a = geno.J_e_s * geno.coil.A_sc / I_s
N_s = int(2 * geno.p * geno.q * N_c_div_a)

print(f"{I_s = }, {N_c_div_a = }, {N_s = }")

#%%
geno.N_s = N_s
geno.J_s_amplitude = np.sqrt(2) * geno.J_e_s

#%% ------ export to comsol -----

save_params("from_python", geno) 

#%%

if 1:

    margin = 0.00
    # geno.fast_flux()
    geno.plt.fluxplot(dr=100, dt=80, lvls=10, lw=0.2,
                      r_min=geno.dims.r_ri-margin, r_max=geno.dims.r_so+margin,
                      t_min=0, t_max=np.pi/4,
                       vmin=0, vmax=7.5,
                      custom_dims=True,
                      y0=0, x0=0,
                      border_width=0.05,
                      show_cbar=False, 
                      show_axis=False,
                      show_borders=True,
                      transparent=True,
                      padding=0,
                      pdf=False, pdf_dpi=300, 
                      svg=True, svg_dpi=300)

#%% check A_tape

N_c = geno.coil.A_sc / gn.A_tape
a = N_c / N_c_div_a

print(f"{N_c = }, {a = }, {N_c_div_a = }, {N_c/a = }")

I_c_with_J = geno.J_e_s * gn.A_tape * 1.7
I_max = I_s / a

print(f"{I_c_with_J = }, {I_max = }")

#%% adapt a

a_new = 8
N_c_div_a_new = int(np.round(N_c / a_new, 2))
N_s_new = int(2 * geno.p * geno.q * N_c / a_new)
print(f"{N_c_div_a_new = }, {N_s_new = }, {geno.N_s = }")

cosPhi = 0.82
I_s_new = geno.P / (np.sqrt(3) * 690 * cosPhi)
I_c_with_J = geno.J_e_s * gn.A_tape * 1.7
I_max = I_s_new / a_new
print(f"{I_c_with_J = }, {I_max = }, {I_s_new = }")

ampere_turns_stator = N_c_div_a_new * I_s_new
J_s_amplitude = np.sqrt(2) * ampere_turns_stator / geno.coil.A_sc
print(f"\n{I_s = } [A] \n{ampere_turns_stator = } [A] \n{J_s_amplitude = } [A/m2] \n{J_s_amplitude*1e-6 = }[A_mm2] \n{geno.J_e_s*1e-6 = }[A_mm2] \n ")

#%% apply adaption

ampere_turns_stator_new = N_c_div_a_new * I_s_new
J_s_amplitude_new = np.sqrt(2) * ampere_turns_stator_new / geno.coil.A_sc
geno.N_s = N_s_new
geno.J_s_amplitude = J_s_amplitude_new