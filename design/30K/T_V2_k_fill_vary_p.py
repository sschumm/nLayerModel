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

""" This script contains an algorithm to iterate over pole pair counts p
    - 1st loop: p
    - 2nd loop: h_wndg_s
    - 3rd loop: h_wndg_r
    - 4th loop: converging to P_target
    
    This script respects the following restrictions
    - Volume = 2
    - ks_d = 0.866
    - B_yoke_max = 2 [T]
    - kr_b adaption
"""

d_so = (2/ gn.thickness_ratio)**(1/3) # ~1.587
r_so = d_so / 2                                # ~0.794
l_e = d_so * gn.thickness_ratio                # ~0.397
# therefore approx:
d_so = 2.0
r_so = 1.0
l_e = 0.5

#%% Prior Data Analysis
if 0:
    # ============================== Initialize Iteration ==============================
    p = 12
    job_init_time = time.time()
    coilsizes = True
    
    gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=3)
    h_pf = gn.h_pole_frame
    gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                        h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
    gen.k_fill_r = 0.5
    gen.k_fill_s = 0.5
    gen.ks_d= 0.866
    # gen.kr_b= 0.65
    
    # ============================== Stator Winding Iteration ==============================    
    lst_stator_iter = []
    h_sw0 = h_pf * 2.1 # 0.01 # 0.01 # h_pf*2.1 #0.3 
    h_rw0 = h_pf * 2.1 # h_pf*2.1 #0.3
    
    n_iters_s = 9
    n_times_s= 2
    n_iters_r = 12
    n_times_r=2
    
    for idx_stator, iter_stator in enumerate(np.linspace(0, n_times_s, n_iters_s)):
        gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_r=h_rw0,
                            h_wndg_s=h_sw0 + h_sw0*iter_stator)
        
        
        # ============================== Rotor Winding Iteration ==============================
        lst_rotor_iter = []
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, n_times_r, n_iters_r)):      
            gen.update_dimensions(h_wndg_r = h_rw0 + h_rw0*iter_rotor)

            
            if coilsizes:
                gen.apply_coil_sizes()
            gen.apply_lift_factor()
            gen.adapt_yokes()

            # gen.apply_coil_thickness_ratio()
            if coilsizes:
                gen.apply_coil_sizes()
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
        # ax2 = ax1.twinx()
        # ax2.scatter([g.h_wndg_r for g in lst_rotor_iter], 
        #          [g.kr_b for g in lst_rotor_iter], color="red")
        
    
    ax1.set_xlabel("h_rw")
    ax1.set_ylabel("P")
    ax1.legend(loc="upper left")
    plt.title(f"Pole Pair Count p = {gen.p}")
    plt.show()

#%% Pole Pair Iteration
# ============================== Initialize Iteration ==============================

if 1:
    job_init_time = time.time()
    total_runs = 0
    coilsizes = True
    
    P_target = 7.2e6 # gn.Pel_out
    
    # ============================== Pole Pair Iteration ==============================
    lst_optimize_l = []
    lst_optimize_m = []
    for p in range(6, 20, 1):
    
        gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw, B_yoke_max=3)
        h_pf = gn.h_pole_frame
        gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                            h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
        gen.ks_d= 0.866
        # gen.kr_b= 0.65
        gen.k_fill_r = 0.7
        gen.k_fill_s = 0.7
    
        # ============================== Stator Winding Iteration ==============================    
        lst_stator_iter = []
        h_sw0 = 0.01 # h_pf*2.1 #0.3 
        h_rw0 = 0.005 # h_pf*2.1 #0.3    
        n_iters_s, n_times_s = 100, 4
        n_iters_r, n_times_r = 30, 8
        detail, verbose = 0.0001, False
        
        for idx_stator, iter_stator in enumerate(np.linspace(0, n_times_s, n_iters_s)):
            # gen.update_dimensions(h_wndg_s= h_sw0 + h_sw0*iter_stator)
            gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_r=h_rw0,
                                h_wndg_s=h_sw0 + h_sw0*iter_stator)
            # gen.maximize_k_fill_s()
            
            # ============================== Rotor Winding Iteration ==============================
            lst_rotor_iter = [copy.deepcopy(gen)]
            for idx_rotor, iter_rotor in enumerate(np.linspace(0, n_times_r, n_iters_r)):
                gen.update_dimensions(h_wndg_r = h_rw0 + h_rw0*iter_rotor)
                # gen.apply_coil_thickness_ratio()
                if coilsizes:
                    gen.apply_coil_sizes()
                gen.apply_lift_factor()
                gen.adapt_yokes()

                # gen.apply_coil_thickness_ratio()
                if coilsizes:
                    gen.apply_coil_sizes()
                gen.apply_lift_factor()
                
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
                            sys.stdout.write(f"\rComputing p = {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Converged.         ")
                            sys.stdout.flush()
        
                            gen.HTS_usage()
                            gen.compute_weight()
                            if gen.h_wndg_r >= h_rw0:
                                lst_stator_iter.append(copy.deepcopy(gen))
                            break
                        
                        sys.stdout.write(f"\rComputing p = {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Convergence Step {idx_Pel_out}         ")
                        sys.stdout.flush()
                        
                        n0 = (n_Pel_out - min(idx_Pel_out, n_Pel_out))/n_Pel_out
                        n1 = min(idx_Pel_out, n_Pel_out)/n_Pel_out
                        
                        if overshoot:
                            if gen.P < P_target:
                                h_rw = old - old * (n0*0.0002 + n1*0.00002)
                            else:
                                h_rw = old + old * (n0*0.00015 + n1*0.000015)
                            gen.update_dimensions(h_wndg_r = h_rw)
                            # gen.apply_coil_thickness_ratio()
                            if coilsizes:
                                gen.apply_coil_sizes()
                            gen.apply_lift_factor()
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
                            # gen.apply_coil_thickness_ratio()
                            if coilsizes:
                                gen.apply_coil_sizes()
                            gen.apply_lift_factor()
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
                            # gen.apply_coil_thickness_ratio()
                            if coilsizes:
                                gen.apply_coil_sizes()
                            gen.apply_lift_factor()
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
        sys.stdout.write(f"\rPole Pairs p = {p} - {len(lst_stator_iter)} configurations." + 100*" " + "\n\n")
        sys.stdout.flush()
        if lst_stator_iter:
            # pick generator with least HTS length
            mini = min([geno.HTS_length for geno in lst_stator_iter])
            idx = [geno.HTS_length for geno in lst_stator_iter].index(mini)
            best = lst_stator_iter[idx]            
            lst_optimize_l.append(best)
            
            # pick generator with least weight
            mini = min([geno.weight for geno in lst_stator_iter])
            idx = [geno.weight for geno in lst_stator_iter].index(mini)
            best = lst_stator_iter[idx]            
            lst_optimize_m.append(best)
            
        
        total_runs += gen.runs
    # ============================== Pole Pair Iteration ==============================
    job_finish_time = time.time()
    job_runtime = job_finish_time - job_init_time
    sys.stdout.write(f"\rJob finished with {total_runs} solved models in {int(np.round(job_runtime,3))//60} minutes and {int(np.round(job_runtime,3))%60} seconds." + " "*80)
    sys.stdout.flush()



#%% --- no export ---

export_png = False
export_tkz = False

#%% --- export ? ---
export_png = True
export_tkz = True
preambel = "1208_1107_30K_kfill_"

#%% ------ plot optim HTS length over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_l]
    lst_HTS_length = [gen.HTS_length for gen in lst_optimize_l]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    # lst_ref_wt = [gen.reference_weight for gen in lst_optimize_l]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "blue", label="weight")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "blue")
    ax1.set_ylabel("Generator Weight / t")
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("HTS length / km")
    ax2.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "red", label="HTS_length")
    ax2.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "red") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green", label="weight old")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "weight_HTS_length.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "weight_HTS_length.tex")

#%% ------ plot weight of rotor and stator ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    lst_weight_rotor = [gen.weight_rotor for gen in lst_optimize_m]
    lst_weight_stator = [gen.weight_stator for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "blue", label="weight")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "blue")
    ax1.set_ylabel("Generator Weight / t")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor], color = "red", label="weight rotor")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor], color = "red")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator], color = "green", label="weight stator")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator], color = "green")
    
    ax1.legend(loc="upper right")
    
    # ax2 = ax1.twinx()
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    # ax2.set_ylabel("Generator Weight / t")
    
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green", label="weight old")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    # ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "weight_rotor_stator.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "weight_rotor_stator.tex")

#%% ------ plot weight of rotor yoke, hts and kryo ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_weight_rotor = [gen.weight_rotor for gen in lst_optimize_m]
    lst_weight_rotor_yoke = [gen.weight_rotor_yoke for gen in lst_optimize_m]
    lst_weight_rotor_wndg = [gen.weight_rotor_wndg for gen in lst_optimize_m]
    lst_weight_rotor_kryo = [gen.weight_rotor_kryo for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_yoke], color = "blue", label="weight rotor yoke")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_yoke], color = "blue")
    ax1.set_ylabel("Generator Weight / t")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_wndg], color = "red", label="weight rotor wndg")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_wndg], color = "red")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_kryo], color = "green", label="weight rotor kryo")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_kryo], color = "green")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor], color = "black", label="weight rotor ")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor], color = "black")
    
    ax1.legend(loc="upper right")
    
    # ax2 = ax1.twinx()
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    # ax2.set_ylabel("Generator Weight / t")
    
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green", label="weight old")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    # ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "weight_rotor_yoke_wndg_kryo.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "weight_rotor_yoke_wndg_kryo.tex")

#%% ------ plot weight of stator yoke, hts and kryo ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_weight_stator = [gen.weight_stator for gen in lst_optimize_m]
    lst_weight_stator_yoke = [gen.weight_stator_yoke for gen in lst_optimize_m]
    lst_weight_stator_wndg = [gen.weight_stator_wndg for gen in lst_optimize_m]
    lst_weight_stator_kryo = [gen.weight_stator_kryo for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_yoke], color = "blue", label="weight stator yoke")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_yoke], color = "blue")
    ax1.set_ylabel("Generator Weight / t")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_wndg], color = "red", label="weight stator wndg")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_wndg], color = "red")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_kryo], color = "green", label="weight stator kryo")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_kryo], color = "green")
    
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator], color = "black", label="weight stator ")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator], color = "black")
    
    ax1.legend(loc="upper right")
    
    # ax2 = ax1.twinx()
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    # ax2.set_ylabel("Generator Weight / t")
    
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green", label="weight old")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    # ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "weight_stator_yoke_wndg_kryo.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "weight_stator_yoke_wndg_kryo.tex")

#%% ------ plot optim HTS length, rotor length, stator length ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_l]
    lst_HTS_length = [gen.HTS_length for gen in lst_optimize_l]
    lst_HTS_length_rotor = [gen.HTS_length_rotor for gen in lst_optimize_l]
    lst_HTS_length_stator = [gen.HTS_length_stator for gen in lst_optimize_l]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "blue", label="total length")    
    ax1.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "blue")
    ax1.set_ylabel("HTS length / km")
    
    ax1.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length_rotor], color = "red", label="length rotor")    
    ax1.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length_rotor], color = "red")
    
    ax1.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length_stator], color = "green", label="length stator")    
    ax1.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length_stator], color = "green")
    
    ax1.legend(loc="upper left")
    
    # ax2 = ax1.twinx()
    # ax2.set_ylabel("HTS length / km")
    # ax2.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "red", label="length")
    # ax2.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "red") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    
    # ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green", label="weight old")    
    # ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_ref_wt], color = "green")
    
    # lines1, labels1 = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "HTS_length_rotor_stator.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "HTS_length_rotor_stator.tex")

#%% ------ plot width of coils and poles over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_w_rc = [gen.coil.w_rc*1e3 for gen in lst_optimize_m]
    lst_w_sc = [gen.coil.w_sc*1e3 for gen in lst_optimize_m]
    lst_w_rp = [gen.coil.w_rp*1e3 for gen in lst_optimize_m]
    lst_w_sp = [gen.coil.w_sp*1e3 for gen in lst_optimize_m]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("widt / mm", color = "blue")
    ax1.scatter(lst_pole_pairs, [l for l in lst_w_rc], color = "blue", label="w_rc")
    ax1.plot(lst_pole_pairs, [l for l in lst_w_rc], color = "blue")
    ax1.scatter(lst_pole_pairs, [l for l in lst_w_sc], color = "green", label="w_sc")
    ax1.plot(lst_pole_pairs, [l for l in lst_w_sc], color = "green")
    ax1.scatter(lst_pole_pairs, [l for l in lst_w_rp], color = "purple", label="w_rp")
    ax1.plot(lst_pole_pairs, [l for l in lst_w_rp], color = "purple")
    ax1.scatter(lst_pole_pairs, [l for l in lst_w_sp], color = "orange", label="w_sp")
    ax1.plot(lst_pole_pairs, [l for l in lst_w_sp], color = "orange")
    
    
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    ax2.set_ylabel("Generator Weight / t")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "vary_p_optim_m.png")
    
#%% ------ plot rotor critical Field, current loading, engineering current density over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_B_r_c = [gen.B_r_c for gen in lst_optimize_m]
    lst_B_s_c = [gen.B_s_c for gen in lst_optimize_m]
    lst_J_e_r = [gen.J_e_r for gen in lst_optimize_m]
    lst_J_e_s = [gen.J_e_s for gen in lst_optimize_m]
    # lst_K = [gen.K_r for gen in lst_optimize_m]
    
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("B_c / T")
    ax1.scatter(lst_pole_pairs, [B for B in lst_B_r_c], color = "blue", label="Brc")
    ax1.plot(lst_pole_pairs, [B for B in lst_B_r_c], color = "blue") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    ax1.scatter(lst_pole_pairs, [B for B in lst_B_s_c], color = "red", label="Bsc")    
    ax1.plot(lst_pole_pairs, [B for B in lst_B_s_c], color = "red")
    ax1.legend(loc="upper left")
    
    # ax2 = ax1.twinx()
    # ax2.scatter(lst_pole_pairs, [w*1e-6 for w in lst_J_e_r], color = "blue", marker= "x",label="J e r")    
    # ax2.plot(lst_pole_pairs, [w*1e-6 for w in lst_J_e_r], color = "blue")
    # ax2.scatter(lst_pole_pairs, [w*1e-6 for w in lst_J_e_s], color = "red", marker= "x",label="J e s")    
    # ax2.plot(lst_pole_pairs, [w*1e-6 for w in lst_J_e_s], color = "red")
    # ax2.set_ylabel("Engineering Current Density / A mm2")
    # ax2.set_ylabel("Engineering Current Density / A mm2")
    
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    if export_png: plt.savefig(fname = preambel + "Brc_Bsc_over_p.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "Brc_Bsc_over_p.tex")
    
#%% ------ plot rotor current loading, engineering current density over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_K_r= [gen.K_r for gen in lst_optimize_m]
    lst_J_e_r = [gen.J_e_r for gen in lst_optimize_m]   
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("K_r / kA m")
    ax1.scatter(lst_pole_pairs, [K*1e-3 for K in lst_K_r], color = "blue", label="Kr")
    ax1.plot(lst_pole_pairs, [K*1e-3 for K in lst_K_r], color = "blue") 
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-6 for w in lst_J_e_r], color = "red", label="J e r")    
    ax2.plot(lst_pole_pairs, [w*1e-6 for w in lst_J_e_r], color = "red")
    ax2.set_ylabel("Engineering Current Density / A mm2")
    
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="center left")
    if export_png: plt.savefig(fname = preambel + "Kr_Jer_over_p.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "Kr_Jer_over_p.tex")
 
#%% ------ plot stator current loading, engineering current density over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_K_s = [gen.K_s for gen in lst_optimize_m]
    lst_J_e_s = [gen.J_e_s for gen in lst_optimize_m]
        
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("K_s / kA m")
    ax1.scatter(lst_pole_pairs, [K*1e-3 for K in lst_K_s], color = "blue", label="Ks")    
    ax1.plot(lst_pole_pairs, [K*1e-3 for K in lst_K_s], color = "blue")
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-6 for w in lst_J_e_s], color = "red", label="J e s")    
    ax2.plot(lst_pole_pairs, [w*1e-6 for w in lst_J_e_s], color = "red")
    ax2.set_ylabel("Engineering Current Density / A mm2")
    
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
    if export_png: plt.savefig(fname = preambel + "Ks_Jes_over_p.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "Ks_Jes_over_p.tex")

#%% ------ plot optim weight and rotor and stator weight  over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_m_rotor = [gen.weight_rotor for gen in lst_optimize_m]
    lst_m_stator = [gen.weight_stator for gen in lst_optimize_m]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("m / t")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_m_rotor], color = "blue", label="m_r")
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_m_rotor], color = "blue") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_m_stator], color = "green", label="m_s")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_m_stator], color = "green")
    ax1.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    ax1.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    
    ax1.legend(loc="upper right")
    if export_png: plt.savefig(fname = preambel + "vary_p_weight_r_and_s.png")

#%% ------ plotfill factors  over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_kr = [gen.k_fill_r for gen in lst_optimize_m]
    lst_ks = [gen.k_fill_s for gen in lst_optimize_m]
        
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("k_fill")
    ax1.scatter(lst_pole_pairs, [K for K in lst_kr], color = "blue", label="k fill r")    
    ax1.plot(lst_pole_pairs, [K for K in lst_kr], color = "blue")
    
    ax1.scatter(lst_pole_pairs, [K for K in lst_ks], color = "red", label="k fill s")    
    ax1.plot(lst_pole_pairs, [K for K in lst_ks], color = "red")
    
    ax1.legend(loc="center left")

    if export_png: plt.savefig(fname = preambel + "k_fills_over_p.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "k_fills_over_p.tex")


#%% ------ plot optim weight and rotor and breadth factor over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_kr_b = [gen.kr_b for gen in lst_optimize_m]
    # lst_m_stator = [gen.weight_stator for gen in lst_optimize_m]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("kr_b / -")
    ax1.scatter(lst_pole_pairs, [B for B in lst_kr_b], color = "blue", label="kr_b")
    ax1.plot(lst_pole_pairs, [B for B in lst_kr_b], color = "blue") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    ax2.set_ylabel("Generator Weight / t")
    
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")  
    if export_png: plt.savefig(fname = preambel + "vary_p_breadth_factor.png")
    
    if export_tkz:
        tkz.clean_figure()
        tkz.save(preambel + "vary_p_breadth_factor.tex")
    
    
#%% ------ plot optim weight and winding height over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_h_rw = [gen.h_wndg_r for gen in lst_optimize_m]
    lst_h_sw = [gen.h_wndg_s for gen in lst_optimize_m]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("h / mm")
    ax1.scatter(lst_pole_pairs, [h*1e3 for h in lst_h_rw], color = "blue", label="h_rw")
    ax1.plot(lst_pole_pairs, [h*1e3 for h in lst_h_rw], color = "blue") #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    ax1.scatter(lst_pole_pairs, [h*1e3 for h in lst_h_sw], color = "green", label="h_sw")    
    ax1.plot(lst_pole_pairs, [h*1e3 for h in lst_h_sw], color = "green")
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red", label="weight new")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight], color = "red")
    ax2.set_ylabel("Generator Weight / t")
    
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    if export_png: plt.savefig(fname = preambel + "vary_p_winding_height.png")
    
#%% ------ plot optim weight, show different weights over pole pairs ------
if 1:
    lst_pole_pairs = [gen.p for gen in lst_optimize_m]
    lst_HTS_length = [gen.HTS_length for gen in lst_optimize_m]
    lst_weight = [gen.weight for gen in lst_optimize_m]
    lst_weight_rotor_yoke = [gen.weight_rotor_yoke for gen in lst_optimize_m]
    lst_weight_stator_yoke = [gen.weight_stator_yoke for gen in lst_optimize_m]
    lst_weight_rotor_wndg = [gen.weight_rotor_wndg for gen in lst_optimize_m]
    lst_weight_stator_wndg = [gen.weight_stator_wndg for gen in lst_optimize_m]
    
    fig = plt.figure(dpi=300, figsize=(8,5))
    ax1 = plt.subplot()
    plt.title("Pole Pair Iteration")
    plt.grid()
    
    ax1.set_xlabel("pole pairs")
    ax1.set_ylabel("HTS length / km", color = "blue")
    ax1.scatter(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length], color = "blue", label="HTS_length")
    ax1.plot(lst_pole_pairs, [l*1e-3 for l in lst_HTS_length]) #, label=f"{np.round(lst_h_wndg_s[idx],3)}")
    
    ax2 = ax1.twinx()
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_yoke], color = "red", label="rotor yoke")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_yoke], color = "red")
    ax2.set_ylabel("Weight / t")
    
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_wndg], color = "green", label="rotor winding")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_rotor_wndg], color = "green")
    
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_wndg], color = "black", label="stator winding")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_wndg], color = "black")
    
    ax2.scatter(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_yoke], color = "purple", label="stator yoke")    
    ax2.plot(lst_pole_pairs, [w*1e-3 for w in lst_weight_stator_yoke], color = "purple")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="upper center")
    
    if export_png: plt.savefig(fname = preambel + "vary_p_optim_m.png")