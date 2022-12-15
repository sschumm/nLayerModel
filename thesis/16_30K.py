# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import sys
import copy
import time
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


""" Untersuchung warum 1m3 zu klein ist """


d_so = (1 / gn.thickness_ratio)**(1/3) # ~1.5874
r_so = d_so / 2                                #  ~0.7938
l_e = d_so * gn.thickness_ratio                #  ~0.397
# therefore approx:
d_so = 1.6
r_so = 0.8
l_e = 0.4

#%% Prior Data Analysis
if 1:
    # ============================== Initialize Iteration ==============================
    p = 6
    job_init_time = time.time()
    
    gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=3)
    h_pf = gn.h_pole_frame
    gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
    gen.ks_d= 0.866
    
    # ============================== Stator Winding Iteration ==============================    
    lst_stator_iter = []
    h_sw0 = h_pf*2.1 #0.3 
    h_rw0 = h_pf*2.1 #0.3
    
    lst = np.array([0.024, 0.026, 0.028, 0.03, 0.032, 0.034])
    
    n_iters_s = 7
    n_times_s= 1.5
    n_iters_r = 20
    n_times_r=3 
    
    detail = 0.0001
    verbose = False
    
    for idx_stator, iter_stator in enumerate(lst): #np.linspace(0, n_times_s, n_iters_s)):
        gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_r=h_rw0,
                            h_wndg_s=iter_stator) #h_sw0 + h_sw0*iter_stator)
        gen.maximize_k_fill()
        
        # ============================== Rotor Winding Iteration ==============================
        lst_rotor_iter = []
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, n_times_r, n_iters_r)):      
            gen.update_dimensions(h_wndg_r = h_rw0 + h_rw0*iter_rotor)
            gen.apply_coil_thickness_ratio()
            gen.apply_lift_factor()
            gen.adapt_yokes()
            gen.apply_coil_thickness_ratio()
            # print(f"{gen.k_fill_r = }, {gen.k_fill_s = }")
            gen.apply_lift_factor()
            lst_rotor_iter.append(copy.deepcopy(gen))
        # ============================== End Rotor Winding Iteration ==============================
        lst_stator_iter.append(lst_rotor_iter)
    # ============================== End Stator Winding Iteration ==============================
    
    job_finish_time = time.time()
    job_runtime = job_finish_time - job_init_time
    total_runs = gen.runs
    sys.stdout.write(f"\rJob finished with {total_runs} solved models in {np.round(job_runtime,3)} seconds." + " "*80)
    sys.stdout.flush()
    
    ##%% Plot iteration of h_sw and h_rw with lift factor
    # fig = plt.figure(dpi=300, figsize=(6,7))
    # ax1 = plt.subplot()
    # plt.grid()
    
    # for lst_rotor_iter in lst_stator_iter:
    #     ax1.scatter([g.h_wndg_r for g in lst_rotor_iter], 
    #              [g.P for g in lst_rotor_iter], 
    #              label=f"h_sw= {np.round(lst_rotor_iter[0].h_wndg_s, 3)}")
    #     ax1.plot([g.h_wndg_r for g in lst_rotor_iter], 
    #              [g.P for g in lst_rotor_iter])
    
    # ax1.set_xlabel("h_rw")
    # ax1.set_ylabel("P")
    # ax1.legend(loc="upper left")
    # plt.title(f"Pole Pair Count p = {gen.p}")
    # # plt.show()

##%%
export_tkz = False
export_png = False

##%%
export_tkz = True
export_png = True
preambel = f"1208_1505_30K_ctr_Vcomparison_p{p}_"

##%%

##%% Plot iteration of h_sw and h_rw with lift factor
fig = plt.figure(dpi=300, figsize=(6,7))
ax1 = plt.subplot()
plt.grid()

for lst_rotor_iter in lst_stator_iter:
    ax1.scatter([g.h_wndg_r  * 1e3 for g in lst_rotor_iter], 
             [g.P * 1e-6 for g in lst_rotor_iter], 
             label=f"hsw= {np.round(lst_rotor_iter[0].h_wndg_s, 3) *1e3}")
    ax1.plot([g.h_wndg_r * 1e3 for g in lst_rotor_iter], 
             [g.P * 1e-6  for g in lst_rotor_iter])

ax1.set_xlabel("hrw")
ax1.set_ylabel("P")
ax1.legend(loc="upper left")
plt.title(f"Pole Pair Count p = {gen.p}")
# plt.show()

if export_png:
    plt.savefig(fname = preambel + "P_over_hrw_hsw.png")

if export_tkz:
    tkz.clean_figure()
    tkz.save(preambel + "P_over_hrw_hsw.tex")

#%%

""" init with d2l = 2m3 """


d_so = (2 / gn.thickness_ratio)**(1/3) # 2 ,~1.5874
r_so = d_so / 2                                # 1, ~0.7938
l_e = d_so * gn.thickness_ratio                # 0.5, ~0.397
# therefore approx:
d_so = 2 # 1.6
r_so = 1 # 0.8
l_e = 0.5 #0.4

#%% Prior Data Analysis
if 0:
    # ============================== Initialize Iteration ==============================
    p = 13
    job_init_time = time.time()
    
    gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=3)
    h_pf = gn.h_pole_frame
    gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
    gen.ks_d= 0.866
    
    # ============================== Stator Winding Iteration ==============================    
    lst_stator_iter = []
    h_sw0 = h_pf*2.1 #0.3 
    h_rw0 = h_pf*2.1 #0.3
    
    n_iters_s = 9
    n_times_s= 2
    n_iters_r = 20
    n_times_r=3 
    
    detail = 0.0001
    verbose = False
    
    for idx_stator, iter_stator in enumerate(np.linspace(0, n_times_s, n_iters_s)):
        gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_r=h_rw0,
                            h_wndg_s=h_sw0 + h_sw0*iter_stator)
        gen.maximize_k_fill()
        
        # ============================== Rotor Winding Iteration ==============================
        lst_rotor_iter = []
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, n_times_r, n_iters_r)):      
            gen.update_dimensions(h_wndg_r = h_rw0 + h_rw0*iter_rotor)
            gen.apply_coil_thickness_ratio()
            gen.apply_lift_factor()
            gen.adapt_yokes()
            gen.apply_coil_thickness_ratio()
            # print(f"{gen.k_fill_r = }, {gen.k_fill_s = }")
            gen.apply_lift_factor()
            lst_rotor_iter.append(copy.deepcopy(gen))
        # ============================== End Rotor Winding Iteration ==============================
        lst_stator_iter.append(lst_rotor_iter)
    # ============================== End Stator Winding Iteration ==============================
    
    job_finish_time = time.time()
    job_runtime = job_finish_time - job_init_time
    total_runs = gen.runs
    sys.stdout.write(f"\rJob finished with {total_runs} solved models in {np.round(job_runtime,3)} seconds." + " "*80)
    sys.stdout.flush()
    
    ##%% Plot iteration of h_sw and h_rw with lift factor
    fig = plt.figure(dpi=300, figsize=(6,7))
    ax1 = plt.subplot()
    plt.grid()
    
    for lst_rotor_iter in lst_stator_iter:
        ax1.scatter([g.h_wndg_r for g in lst_rotor_iter], 
                 [g.P for g in lst_rotor_iter], 
                 label=f"h_sw= {np.round(lst_rotor_iter[0].h_wndg_s, 3)}")
        ax1.plot([g.h_wndg_r for g in lst_rotor_iter], 
                 [g.P for g in lst_rotor_iter])
    
    ax1.set_xlabel("h_rw")
    ax1.set_ylabel("P")
    ax1.legend(loc="upper left")
    plt.title(f"Pole Pair Count p = {gen.p}")
    plt.show()