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


d_so = (sw.d2so_L / gn.thickness_ratio)**(1/3) # ~1.587
r_so = d_so / 2                                # ~0.794
l_e = d_so * gn.thickness_ratio                # ~0.397
# therefore approx:
d_so = 1.6
r_so = 0.8
l_e = 0.4
    
    
#%%

if 1:
    coilsizes = True
    p = 10
    
    geno = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=6)
    h_pf = gn.h_pole_frame
    geno.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=0.01, h_wndg_r=0.01)
    # geno.ks_d= 0.866
    
    geno.update_model_by_K(K_s = geno.k_fill_s * geno.h_wndg_s * sw.J_e,
                                K_r = geno.k_fill_r * geno.h_wndg_r * fw.J_e,
                                ks_d = 0.866)

    geno.k_fill_r = 0.45
    geno.k_fill_s = 0.4
    
    # if coilsizes:
    #     geno.apply_coil_sizes()
    # geno.apply_lift_factor()
    # geno.adapt_yokes()

    # # geno.apply_coil_thickness_ratio()
    # if coilsizes:
    #     geno.apply_coil_sizes()
    # geno.apply_lift_factor()
    
    # geno.fast_flux()
    # geno.show_results()
    
#%%
    geno.update_dimensions(h_wndg_s = 0.03, h_wndg_r=0.03)
    # geno.apply_coil_thickness_ratio()
    if coilsizes:
        geno.apply_coil_sizes()
    geno.apply_lift_factor()
    geno.adapt_yokes()
    
    if coilsizes:
        geno.apply_coil_sizes()
    geno.apply_lift_factor()

#%%        
    if 1:
        increasing_h_by_factor=0.03
        break_when_P_at_factor=0.5
    
        lst_P, lst_h_wndg_r, lst_Brc, lst_r_bend = [geno.P], [0.03], [geno.B_r_c], [geno.coil.r_r_bend]
        for idx_h_wndg_r in range(40):
            h_wndg_r = lst_h_wndg_r[-1] + lst_h_wndg_r[0] * increasing_h_by_factor
            geno.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=False)
            if coilsizes:
                geno.apply_coil_sizes(verbose=True)
            geno.apply_lift_factor()
            geno.adapt_yokes()
        
            # geno.apply_coil_thickness_ratio()
            if coilsizes:
                geno.apply_coil_sizes()
            geno.apply_lift_factor()
            
            lst_P.append(geno.P)
            lst_h_wndg_r.append(h_wndg_r)
            lst_Brc.append(geno.B_r_c)
            lst_r_bend.append(geno.coil.r_r_bend)
                    
            # if lst_P[-1] < lst_P[0] * break_when_P_at_factor:
            #     break
        
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
        ax2.scatter([x for x in range(len(lst_h_wndg_r))], [h*1e3 for h in lst_r_bend], marker=".", c="g")
        ax2.plot([x for x in range(len(lst_h_wndg_r))], [h*1e3 for h in lst_r_bend], c="g")
        ax2.set_ylabel('h_wndg_r [mm]', color="r")
#%%        
        fig, ax1 = plt.subplots()
        plt.grid()
        ax2 = ax1.twinx()
        fig.tight_layout()         
        fig.dpi=500
        ax1.scatter([x for x in range(len(lst_P))], [y*1e-6 for y in lst_P], marker=".", c="b")
        ax1.plot([x for x in range(len(lst_P))], [y*1e-6 for y in lst_P], c="b")
        ax1.set_ylabel('P [MW]', color="b")
        ax2.scatter([x for x in range(len(lst_Brc))], lst_Brc, marker=".", c="r")
        ax2.plot([x for x in range(len(lst_Brc))], lst_Brc, c="r")
        ax2.set_ylabel('B_r_c [T]', color="r")
        
        # tkz.clean_figure()
        # tkz.save("aus_P_and_hrw_over_iter.tex")
        
    # geno.update_dimensions(h_wndg_r = 0.030, keep_const_kr_b=False)
    # geno.apply_coil_thickness_ratio()
    # geno.apply_lift_factor()
    # geno.show_results("this")
#%%    

geno.update_dimensions(h_wndg_r=0.037, h_wndg_s=0.025)
if coilsizes:
    geno.apply_coil_sizes()
geno.apply_lift_factor()
geno.adapt_yokes()

# geno.apply_coil_thickness_ratio()
if coilsizes:
    geno.apply_coil_sizes()
geno.apply_lift_factor()


if 1:
    increasing_h_by_factor=0.08
    break_when_P_at_factor=0.5

    lst_P, lst_h_wndg_s, lst_Bsc, lst_s_bend = [geno.P], [geno.h_wndg_s], [geno.B_s_c], [geno.coil.r_s_bend]
    for idx_h_wndg_s in range(40):
        h_wndg_s = lst_h_wndg_s[-1] + lst_h_wndg_s[0] * increasing_h_by_factor
        geno.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=False)
        if coilsizes:
            geno.apply_coil_sizes(verbose=True)
        geno.apply_lift_factor()
        geno.adapt_yokes()
    
        # geno.apply_coil_thickness_ratio()
        if coilsizes:
            geno.apply_coil_sizes()
        geno.apply_lift_factor()
        
        # if geno.coil.r_s_bend > 0.01: # geno.gn.r_bend_max: 
        #     print("Valid") 
        # else:  print("Invalid")
        lst_P.append(geno.P)
        lst_h_wndg_s.append(h_wndg_s)
        lst_Bsc.append(geno.B_s_c)
        lst_s_bend.append(geno.coil.r_s_bend)
                
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
    ax2.scatter([x for x in range(len(lst_h_wndg_s))], [h*1e3 for h in lst_s_bend], marker=".", c="g")
    ax2.plot([x for x in range(len(lst_h_wndg_s))], [h*1e3 for h in lst_s_bend], c="g")
    # ax2.plot([x for x in range(len(lst_Bsc))], lst_Bsc, c="r")
    # ax2.set_ylabel('B_r_c [T]', color="r")
    
    # tkz.clean_figure()
    # tkz.save("aus_P_and_hsw_over_iter.tex")
    
    
#%%    
geno.update_dimensions(h_wndg_s = 0.044)
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.adapt_yokes()
geno.apply_coil_thickness_ratio()
geno.apply_lift_factor()
geno.show_results()