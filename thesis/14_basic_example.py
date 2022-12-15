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

    
    
#%%

if 1:
    p = 30
    l_e = 0.94
    r_so = 2 #sw.r_so(sw, l_e)
    
    geno = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=1.7)
    h_pf = 0.02 * r_so
    geno.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf/4, h_wndg_r=h_pf/4)
    
    geno.update_model_by_K(K_s = geno.k_fill_s * geno.h_wndg_s * sw.J_e,
                           K_r = geno.k_fill_r * geno.h_wndg_r * fw.J_e)
    geno.ks_d = 0.9577
    geno.kr_b = 0.8916
    
    geno.k_fill_r = 0.6
    geno.k_fill_s = 0.6
    
    geno.apply_lift_factor()   
    geno.adapt_yokes()
    geno.apply_lift_factor()
    
    geno.show_results()
    
    
#%%
geno.update_dimensions(h_wndg_s = 0.027)
    
    
#%%

geno.update_dimensions(h_wndg_r = 0)
geno.apply_lift_factor()   
geno.adapt_yokes()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.show_results()

#%%        
if 1:
    h_rw0 = 0.01
    increasing_h_by_factor=0.2
    break_when_P_at_factor=0.5

    lst_P, lst_h_wndg_r, lst_Brc, lst_Jer, lst_hys, lst_hyr, lst_Jes, lst_Bsc = [], [], [], [], [], [], [], []
    for idx_h_wndg_r in range(40):
        h_wndg_r = h_rw0 * increasing_h_by_factor * idx_h_wndg_r
        geno.update_dimensions(h_wndg_r = h_wndg_r)
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_lift_factor()
        
        lst_P.append(geno.P)
        lst_h_wndg_r.append(h_wndg_r)
        lst_Brc.append(geno.B_r_c)
        lst_Bsc.append(geno.B_s_c)
        lst_Jer.append(geno.J_e_r)
        lst_Jes.append(geno.J_e_s)
        lst_hys.append(geno.h_yoke_s)
        lst_hyr.append(geno.h_yoke_r)
                
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
    ax2.set_ylabel('hwndgr [m]', color="r")
    tkz.clean_figure()
    tkz.save("aus_basics_77K_Pel_hrw.tex")
    plt.show()

    fig, ax1 = plt.subplots()
    plt.grid()
    # ax2 = ax1.twinx()
    fig.tight_layout()         
    fig.dpi=500
    ax1.scatter([x for x in range(len(lst_hys))], [y*1e3 for y in lst_hys], marker=".", c="b", label="hys")
    ax1.plot([x for x in range(len(lst_hys))], [y*1e3 for y in lst_hys], c="b")
    ax1.set_ylabel('J [Amm2]', color="b")
    ax1.scatter([x for x in range(len(lst_hyr))], [y*1e3 for y in lst_hyr], marker=".", c="r", label="hyr")
    ax1.plot([x for x in range(len(lst_hyr))], [y*1e3 for y in lst_hyr], c="r")
    ax1.legend()
    tkz.clean_figure()
    tkz.save("aus_basics_77K_hyr_hys.tex")
    plt.show()

    fig, ax1 = plt.subplots()
    plt.grid()
    ax2 = ax1.twinx()
    fig.tight_layout()         
    fig.dpi=500
    ax1.scatter([x for x in range(len(lst_Jer))], [y*1e-6 for y in lst_Jer], marker=".", c="b")
    ax1.plot([x for x in range(len(lst_Jer))], [y*1e-6 for y in lst_Jer], c="b")
    ax1.set_ylabel('J [Amm2]', color="b")
    ax2.scatter([x for x in range(len(lst_Brc))], [y for y in lst_Brc], c="r", marker=".")
    ax2.plot([x for x in range(len(lst_Brc))], [y for y in lst_Brc], c="r")
    ax2.set_ylabel('J [Amm2]', color="r")
    tkz.clean_figure()
    tkz.save("aus_basics_77K_Jer_Brc.tex")
    plt.show()
    
    
    fig, ax1 = plt.subplots()
    plt.grid()
    ax2 = ax1.twinx()
    fig.tight_layout()         
    fig.dpi=500
    ax1.scatter([x for x in range(len(lst_Jes))], [y*1e-6 for y in lst_Jes], marker=".", c="b")
    ax1.plot([x for x in range(len(lst_Jes))], [y*1e-6 for y in lst_Jes], c="b")
    ax1.set_ylabel('J [Amm2]', color="b")
    ax2.scatter([x for x in range(len(lst_Bsc))], [y for y in lst_Bsc], c="r", marker=".")
    ax2.plot([x for x in range(len(lst_Bsc))], [y for y in lst_Bsc], c="r")
    ax2.set_ylabel('B [T]', color="r")
    tkz.clean_figure()
    tkz.save("aus_basics_77K_Jes_Bsc.tex")
    plt.show()
   



    
    # tkz.clean_figure()
    # tkz.save("aus_P_and_hrw_over_iter.tex")
        
#%%    

geno.update_dimensions(h_wndg_r=0.02, h_wndg_s=0)
geno.apply_lift_factor()   
geno.adapt_yokes()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.show_results()

#%%
if 1:
    h_sw0 = 0.01
    increasing_h_by_factor=0.2
    break_when_P_at_factor=0.5

    lst_P, lst_h_wndg_s, lst_Bsc, lst_Jes, lst_hys, lst_hyr = [], [], [], [], [], []
    for idx_h_wndg_s in range(40):
        h_wndg_s = h_sw0 * increasing_h_by_factor * idx_h_wndg_s
        geno.update_dimensions(h_wndg_s = h_wndg_s)
        geno.apply_lift_factor()
        geno.adapt_yokes()
        geno.apply_lift_factor()
        
        lst_P.append(geno.P)
        lst_h_wndg_s.append(h_wndg_s)
        lst_Bsc.append(geno.B_s_c)
        lst_Jes.append(geno.J_e_s)
        lst_hys.append(geno.h_yoke_s)
        lst_hyr.append(geno.h_yoke_r)
                
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
    ax2.set_ylabel('hwndgr [m]', color="r")
    plt.show()
    
    fig, ax1 = plt.subplots()
    plt.grid()
    # ax2 = ax1.twinx()
    fig.tight_layout()         
    fig.dpi=500
    ax1.scatter([x for x in range(len(lst_hys))], [y for y in lst_hys], marker=".", c="b", label="hys")
    ax1.plot([x for x in range(len(lst_hys))], [y for y in lst_hys], c="b")
    ax1.set_ylabel('J [Amm2]', color="b")
    ax1.scatter([x for x in range(len(lst_hyr))], [y for y in lst_hyr], marker=".", c="r", label="hyr")
    ax1.plot([x for x in range(len(lst_hyr))], [y for y in lst_hyr], c="r")
    # ax2.scatter([x for x in range(len(lst_Brc))], lst_Brc, c="r", marker=".")
    # ax2.plot([x for x in range(len(lst_Brc))], lst_Brc, c="r")
    # ax2.set_ylabel('Brc [T]', color="r")
    ax1.legend()
    plt.show()

    # fig, ax1 = plt.subplots()
    # plt.grid()
    # ax2 = ax1.twinx()
    # fig.tight_layout()         
    # fig.dpi=500
    # ax1.scatter([x for x in range(len(lst_P))], [y*1e-6 for y in lst_Jes], marker=".", c="b")
    # ax1.plot([x for x in range(len(lst_P))], [y*1e-6 for y in lst_Jes], c="b")
    # ax1.set_ylabel('J [Amm2]', color="b")
    # ax2.scatter([x for x in range(len(lst_Bsc))], lst_Bsc, c="r", marker=".")
    # ax2.plot([x for x in range(len(lst_Bsc))], lst_Bsc, c="r")
    # ax2.set_ylabel('Bsc [T]', color="r")
    # plt.show()
    
    # tkz.clean_figure()
    # tkz.save("aus_P_and_hsw_over_iter.tex")
    
    
#%%    
geno.update_dimensions(h_wndg_s = 0.027, h_wndg_r=0.03)
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_lift_factor()
geno.show_results()

# geno.fast_flux(vmin=0, vmax=4)

#%% ----- select generator and finish setup for export ------

if 1:
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

#%% ------ export to comsol -----

save_params("from_python", geno) 

#%%

if 1:

    margin = 0.00
    # geno.fast_flux()
    geno.plt.fluxplot(dr=100, dt=80, lvls=10, lw=1,
                      r_min=geno.dims.r_ri-margin, r_max=geno.dims.r_so+margin,
                      t_min=0, t_max=np.pi/15,
                      vmin=0, vmax=4,
                      custom_dims=True,
                      y0=0, x0=0,
                      border_width=0.2,
                      show_cbar=False, 
                      show_axis=False,
                      show_borders=True,
                      transparent=True,
                      padding=0,
                      pdf=False, pdf_dpi=300, 
                      svg=True, svg_dpi=300)
    
    
#%% adapt N_s

N_c = geno.coil.A_sc / gn.A_tape
a = N_c / N_c_div_a

print(f"{N_c = }, {a = }, {N_c_div_a = }, {N_c/a = }")

I_c_with_J = geno.J_e_s * gn.A_tape * 1.7
I_c_with_N = I_s / a

print(f"{I_c_with_J = }, {I_c_with_N = }")
#%%

a_new = 60
N_c_div_a_new = int(np.round(N_c / a_new, 2))
N_s_new = int(2 * geno.p * geno.q * N_c / a_new)
print(f"{N_c_div_a_new = }, {N_s_new = }")

I_c_with_J = geno.J_e_s * gn.A_tape * 1.7
I_c_with_N = I_s / a_new
print(f"{I_c_with_J = }, {I_c_with_N = }")

ampere_turns_stator_new = N_c_div_a_new * I_s * 1.1

J_s_amplitude_new = np.sqrt(2) * ampere_turns_stator_new / geno.coil.A_sc
geno.N_s = N_s_new
geno.J_s_amplitude = J_s_amplitude_new

print(f"\n{I_s = } [A] \n{ampere_turns_stator_new = } [A] \n{J_s_amplitude_new = } [A/m2] \n{J_s_amplitude_new*1e-6 = }[A_mm2]")

#%%
from scipy.constants import mu_0

def en_dens_trans_curr(H_p, H_m, i_m):
        
    e_ih = 2*mu_0*H_p**2 * (
        ((H_p * (3 + i_m**2)) / (3 * H_m)) -
        ((2 * H_p**2 * (1-i_m**3)) / (3 * H_m**2)) + 
        ((6*H_p**3 * i_m**2 * (1-i_m)**2) / (3 * H_m**2 * (H_m - H_p * i_m))) +
        ((6 * H_p**3 * i_m**2 * (1-i_m)**2) / (3*H_m**2 * (H_m - H_p * i_m))) -
        ((4 * H_p**4 * i_m**2 * (1-i_m)**3) / (3 * H_m**2 * (H_m - H_p * i_m)**2))
        )
    
    return e_ih

H_m = 2 / mu_0
H_p = geno.J_e_s * 1.7 * 1e-3 * 6  
i_m = 63.86 / 181.85 # I_m / I_c

e_ih = en_dens_trans_curr(H_p, H_m, i_m)
print(f"{e_ih = }")

geno.HTS_usage()
V_s = geno.HTS_volume_stator
print(f"{V_s = }")

P_ih = e_ih * geno.p * geno.gn.n_syn/60 * V_s
print(f"{P_ih = }")




