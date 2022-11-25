# -*- coding: utf-8 -*-
import sys
import copy
import time
import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from data.export import save_params

from design import n7_Model


d_so = (sw.d2so_L / gn.thickness_ratio)**(1/3) # ~3.915
r_so = d_so / 2                                # ~1.9575
l_e = d_so * gn.thickness_ratio                # ~0.9787
# therefore approx:
d_so = 4
r_so = 2
l_e = 0.94


#%% Prior Data Analysis
if 0:
    # ============================== Initialize Iteration ==============================
    p = 30
    job_init_time = time.time()
    
    gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw)
    h_pf = gn.h_pole_frame
    gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
    gen.ks_d= 0.866
    
    # ============================== Stator Winding Iteration ==============================    
    lst_stator_iter = []
    h_sw0 = h_pf*2.1 #0.3 
    h_rw0 = h_pf*2.1 #0.3
    
    n_iters_s = 9
    n_times_s= 3
    n_iters_r = 12
    n_times_r=3 
    
    detail = 0.0001
    verbose = False
    
    for idx_stator, iter_stator in enumerate(np.linspace(0, n_times_s, n_iters_s)):
        gen.update_dimensions(h_wndg_s= h_sw0 + h_sw0*iter_stator)
        gen.maximize_k_fill_s()
        
        # ============================== Rotor Winding Iteration ==============================
        lst_rotor_iter = []
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, n_times_r, n_iters_r)):      
            gen.update_dimensions(h_wndg_r = h_rw0 + h_rw0*iter_rotor)
            gen.maximize_k_fill_r()
            
            gen.apply_coil_sizes_and_lift_factor()
            gen.adapt_yokes()
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

#%% Pole Pair Iteration
# ============================== Initialize Iteration ==============================

job_init_time = time.time()
total_runs = 0

p = 6
P_target = gn.Pel_out

# ============================== Pole Pair Iteration ==============================
lst_p_iter = []
for p in range(8, 61, 5):

    gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw)
    h_pf = gn.h_pole_frame
    gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
    gen.ks_d= 0.866

    # ============================== Stator Winding Iteration ==============================    
    lst_stator_iter = []
    h_sw0 = h_pf*2.1 #0.3 
    h_rw0 = h_pf*2.1 #0.3    
    n_iters_s, n_times_s = 9, 3
    n_iters_r, n_times_r = 12, 3 
    detail, verbose = 0.0001, False
    
    for idx_stator, iter_stator in enumerate(np.linspace(0, n_times_s, n_iters_s)):
        # gen.update_dimensions(h_wndg_s= h_sw0 + h_sw0*iter_stator)
        gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_r=h_pf * 2.1,
                            h_wndg_s=h_sw0 + h_sw0*iter_stator)
        gen.maximize_k_fill_s()
        
        # ============================== Rotor Winding Iteration ==============================
        lst_rotor_iter = [copy.deepcopy(gen)]
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, n_times_r, n_iters_r)):
            gen.update_dimensions(h_wndg_r = h_rw0 + h_rw0*iter_rotor)
            gen.maximize_k_fill_r()
            gen.apply_coil_sizes_and_lift_factor()
            gen.adapt_yokes()
            
            if gen.P >= P_target:            
                # ============================== Find Pel_out Iteration ==============================
                x = []
                y = []
                
                streak_pos, streak_neg = 1, 1
                n_Pel_out = 40
                overshoot = False
                old = copy.copy(gen.h_wndg_r)
                for idx_Pel_out in range(100):
                    if (gen.P >= (1-detail)*P_target and gen.P <= (1+detail)*P_target):
                        sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Converged.         ")
                        sys.stdout.flush()
    
                        gen.HTS_usage()
                        gen.compute_weight()
                        lst_stator_iter.append(copy.deepcopy(gen))
                        break
                    
                    sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Convergence Step {idx_Pel_out}         ")
                    sys.stdout.flush()
                    
                    n0 = (n_Pel_out - min(idx_Pel_out, n_Pel_out))/n_Pel_out
                    n1 = min(idx_Pel_out, n_Pel_out)/n_Pel_out
                    
                    if overshoot:
                        if gen.P < P_target:
                            h_rw = old - old * (n0*0.0002 + n1*0.00002)
                        else:
                            h_rw = old + old * (n0*0.00015 + n1*0.000015)
                        gen.update_dimensions(h_wndg_r = h_rw)
                        gen.maximize_k_fill_r()
                        gen.apply_coil_sizes_and_lift_factor()
                        gen.adapt_yokes()   
                        overshoot = False
                        if verbose:
                            x.append(idx_Pel_out)
                            y.append(gen.P)
                        continue
                    
                    old = copy.copy(gen.h_wndg_r)
                    if gen.P > P_target:
                        streak_neg = 1                    
                        h_rw = gen.h_wndg_r - gen.h_wndg_r * ((n0*0.002 + n1*0.00002) * streak_pos)
                        streak_pos *= 2
                        
                        gen.update_dimensions(h_wndg_r = h_rw)
                        gen.maximize_k_fill_r()
                        gen.apply_coil_sizes_and_lift_factor()
                        gen.adapt_yokes()
                        
                        if gen.P < P_target:
                            streak_pos = 1
                            overshoot = True  
                            continue
                    else:
                        streak_pos = 1
                        h_rw = gen.h_wndg_r + gen.h_wndg_r * ((n0*0.0015 + n1*0.000015) * streak_neg)
                        streak_neg *=2
    
                        gen.update_dimensions(h_wndg_r = h_rw)
                        gen.maximize_k_fill_r()
                        gen.apply_coil_sizes_and_lift_factor()
                        gen.adapt_yokes()
                        
                        if gen.P > P_target:
                            streak_neg = 1 
                            overshoot = True   
                            continue
                    if verbose:
                        x.append(idx_Pel_out)
                        y.append(gen.P)
                    
                if verbose:
                    fig = plt.figure(dpi=300, figsize=(6,7))
                    ax1 = plt.subplot()
                    ax1.plot(x,y)
                    ax1.scatter(x,y)
                    plt.grid()
                    plt.show()
                break
                # ============================== End Pel_out Iteration ==============================
            else:
                lst_rotor_iter.append(copy.deepcopy(gen))
                if len(lst_rotor_iter)>1 and lst_rotor_iter[-1].P < lst_rotor_iter[-2].P:
                    sys.stdout.write(f"\rComputing p = {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: BREAK                                     ")
                    sys.stdout.flush()
                    break
                else:
                    sys.stdout.write(f"\rComputing p = {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r}                                       ")
                    sys.stdout.flush()
            pass
        # ============================== End Rotor Winding Iteration ==============================
        pass
    # ============================== End Stator Winding Iteration ==============================
    sys.stdout.write(f"\r\nPole Pairs p = {p} - {len(lst_stator_iter)} configurations." + 100*" " + "\n\n")
    sys.stdout.flush()
    if lst_stator_iter:
        # pick generator with least HTS length
        mini = min([geno.HTS_length for geno in lst_stator_iter])
        idx = [geno.HTS_length for geno in lst_stator_iter].index(mini)
        best = lst_stator_iter[idx]            
        lst_p_iter.append(best)
    
    total_runs += gen.runs
# ============================== Pole Pair Iteration ==============================
job_finish_time = time.time()
job_runtime = job_finish_time - job_init_time
sys.stdout.write(f"\rJob finished with {total_runs} solved models in {np.round(job_runtime,3)} seconds." + " "*80)
sys.stdout.flush()


#%% ------ plot HTS length over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_p_iter]
    lst_HTS_length = [gen.HTS_length for gen in lst_p_iter]
    lst_weight = [gen.weight for gen in lst_p_iter]
    
    fig = plt.figure(dpi=300, figsize=(10,7))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("HTS length / km", color = "blue")
    ax1.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "blue")
    ax1.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "blue") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    ax2.set_ylabel("Generator Weight / t", color = "red")