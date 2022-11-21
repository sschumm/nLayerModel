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


# --- ecoswing ---
p = 8
r_so = 3.75
l_e = sw.l_e(sw, r_so) # 0.267


generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     k_fill_s=0.4,
                     k_fill_r=0.4,
                     B_yoke_max=1.7)


# --- init generator dimension via yoke and winding heights ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.05, 1.05]]

generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
    
# --- compute initial current loadings and with it the initial model ---
K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.866)
generator.keep_const_kr_b(kr_b = 0.6)
generator.apply_coil_sizes_and_lift_factor(verbose = False)

def adapt_yokes():
    # --- adapt stator yoke ---
    h_yoke_s = yoke_height(generator.p, generator.dims.r_si, 
                            generator.B_s, generator.B_yoke_max)    
    # --- adapt rotor yoke ---
    h_yoke_r = yoke_height(generator.p, generator.dims.r_ro, 
                            generator.B_r, generator.B_yoke_max)

    generator.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r, keep_const_kr_b=True)
    generator.apply_coil_sizes_and_lift_factor(verbose = False)

adapt_yokes()
# generator.show_results("Results with adapted yokes")
# generator.fast_flux()  


#%% ------ visualize correlation of power with h_wndg_s and h_wndg_r ------ 
if 0:
    lst_P = []
    lst_h_wndg_s = []
    lst_h_wndg_r = []
    
    h_wndg_s_0 = 2 * gn.h_pole_frame * 1.3
    h_wndg_s_1 = h_wndg_s_0 * 2
    
    h_wndg_r_0 = 2 * gn.h_pole_frame * 1.1
    h_wndg_r_1 = h_wndg_r_0 * 2
    
    
    n_iters_s = 9
    n_iters_r = 20 
    
    for idx_stator, iter_stator in enumerate(np.linspace(0, h_wndg_s_1-h_wndg_s_0, n_iters_s)):
        
        h_wndg_s = h_wndg_s_0 + iter_stator
        generator.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=True)
        
        hist_r_P = []
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, h_wndg_r_1-h_wndg_r_0, n_iters_r)):
            
            h_wndg_r = h_wndg_r_0 + iter_rotor
            
            generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
            generator.apply_coil_sizes_and_lift_factor()
            adapt_yokes()
            
            hist_r_P.append(generator.P)
            
            sys.stdout.write(f"\rStator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r}    ")
            sys.stdout.flush()
            
            if idx_stator == 0:
                lst_h_wndg_r.append(h_wndg_r)
    
        lst_P.append(hist_r_P)
        lst_h_wndg_s.append(h_wndg_s)

#%% ------ plot ------
if 0:
    fig = plt.figure(dpi=300, figsize=(6,7))
    ax1 = plt.subplot()
    # ax2 = ax1.twinx()
    # ax.set_xlabel("B / T")
    # ax.set_ylabel("J_e / Amm2")
    # ax.set_xlim(0, 3)
    # ax.set_ylim(0, 25)
    # ax.set_xticks([i for i in range(0, 181, 20)])
    # ax.set_yticks([i for i in range(0, 2501, 250)])
    plt.title(f"p = {generator.p}")
    plt.grid()
    
    for idx, val in enumerate(lst_P):
        ax1.scatter(lst_h_wndg_r, [v*1e-6 for v in val], label=f"{np.round(lst_h_wndg_s[idx],2)}")
        ax1.plot(lst_h_wndg_r, [v*1e-6 for v in val], label=f"{np.round(lst_h_wndg_s[idx],3)}")
        
    ax1.legend(loc="lower right", ncol=3)
    
    
    # ax1.plot(lst_h_wndg_s, lst_kr_b, label="kr_b", color = "blue")
    # ax1.plot(lst_h_wndg_s, lst_k_fill_s, label="k_fill_s", color = "red")
    
    # ax2.plot(lst_h_wndg_s, lst_P, label="P", color = "green")
    # ax2.plot(lst_h_wndg_s, lst_K_s , label="K_s", color = "black")
    
    # lines1, labels1 = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines1 + lines2, labels1 + labels2, loc=0)
    plt.show()
    pass
    
#%% ------ break at 7 MW crossing ------
if 0:
    lst_P = []
    lst_h_wndg_s = []
    lst_h_wndg_r = []
    
    h_wndg_s_0 = 2 * gn.h_pole_frame * 1.3
    h_wndg_s_1 = h_wndg_s_0 * 2
    
    h_wndg_r_0 = 2 * gn.h_pole_frame * 1.1
    h_wndg_r_1 = h_wndg_r_0 * 2
    
    
    n_iters_s = 9
    n_iters_r = 20 
    
    detail = 0.0001
    verbose = True
    
    # ============================== Stator Winding Iteration ==============================
    for idx_stator, iter_stator in enumerate(np.linspace(0, h_wndg_s_1-h_wndg_s_0, n_iters_s)):
        
        h_wndg_s = h_wndg_s_0 + iter_stator
        generator.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=True)
        
        
        # ============================== Rotor Winding Iteration ==============================
        hist_r_P = []
        for idx_rotor, iter_rotor in enumerate(np.linspace(0, h_wndg_r_1-h_wndg_r_0, n_iters_r)):
            
            h_wndg_r = h_wndg_r_0 + iter_rotor
            
            generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
            generator.apply_coil_sizes_and_lift_factor()
            adapt_yokes()
            
            if generator.P >= gn.Pel_out:
                sys.stdout.write(f"\rStator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Searching config...   ")
                sys.stdout.flush()
                
                x = []
                y = []
                config = [h_wndg_r, h_wndg_r - (h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r - 2*(h_wndg_r_1-h_wndg_r_0)/n_iters_r]
                n = 15
                for i in range(50):
                    if (generator.P >= (1-detail)*gn.Pel_out and generator.P <= (1+detail)*gn.Pel_out):
                        lst_P.append(generator.P)
                        lst_h_wndg_s.append(generator.h_wndg_s)
                        lst_h_wndg_r.append(generator.h_wndg_r)
                        break
                    j = min(i, n)
                    fac = 5 * np.abs(gn.Pel_out - generator.P)/gn.Pel_out  
                    if generator.P > gn.Pel_out:
                        h_wndg_r= 0.25*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[1]) + 0.1*(((n-j)/n)*config[1]+(j/n)*h_wndg_r) + 0.1*(((n-j)/n)*config[2]+(j/n)*h_wndg_r) + 0.05*(((n-j)/n)*config[2]+(j/n)*((1-fac)*h_wndg_r + fac*config[2]))
                    else:
                        h_wndg_r= 0.35*h_wndg_r +0.55*((1-fac)*h_wndg_r + fac*config[0]) + 0.1*(((n-j)/n)*config[0]+(j/n)*h_wndg_r)
                        
                    generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
                    generator.apply_coil_sizes_and_lift_factor()
                    adapt_yokes()
                    if verbose:
                        x.append(i)
                        y.append(generator.P)
                    
                if verbose:
                    fig = plt.figure(dpi=300, figsize=(6,7))
                    ax1 = plt.subplot()
                    ax1.plot(x,y)
                    ax1.scatter(x,y)
                    plt.grid()
                    plt.show()
                break
            else:
                hist_r_P.append(generator.P)
                if len(hist_r_P)>1 and hist_r_P[-1] < hist_r_P[-2]:
                    sys.stdout.write(f"\rStator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: BREAK    ")
                    sys.stdout.flush()
                    break
                else:
                    sys.stdout.write(f"\rStator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r}       ")
                    sys.stdout.flush()
            

#%% ------ iterating over p ------
if 1:
    job_init_time = time.time()
    r_so = 3.75
    l_e = sw.l_e(sw, r_so) # 0.267
    # l_e = 1.142 
    # r_so = sw.r_so(sw, l_e) # 1.812
    total_runs = 0
    
    # ============================== Pole Pair Iteration ==============================
    lst_best_generators = []
    for p in range(8, 61, 1):

        generator = n7_Model(p, l_e, r_so, 
                             gn=gn, fw=fw, sw=sw,
                             k_fill_s=0.3,
                             k_fill_r=0.3,
                             B_yoke_max=1.7)

        # --- init generator dimension via yoke and winding heights ---
        h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
        h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.05, 1.05]]
        generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
            
        # --- compute initial current loadings and with it the initial model ---
        K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
        K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
        generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.866)
        generator.keep_const_kr_b(kr_b = 0.65)
        generator.apply_coil_sizes_and_lift_factor(verbose = False)
        adapt_yokes()


        # ============================== Stator Winding Iteration ==============================    
        lst_generators = []
        
        h_wndg_s_0 = 2 * gn.h_pole_frame * 1.1
        h_wndg_s_1 = h_wndg_s_0 * 2.5
        
        h_wndg_r_0 = 2 * gn.h_pole_frame * 1.1
        h_wndg_r_1 = h_wndg_r_0 * 2
        
        n_iters_s = 20
        n_iters_r = 20 
        
        detail = 0.0001
        verbose = False
        
        for idx_stator, iter_stator in enumerate(np.linspace(0, h_wndg_s_1-h_wndg_s_0, n_iters_s)):
            
            h_wndg_s = h_wndg_s_0 + iter_stator
            generator.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=True)
            
            
            # ============================== Rotor Winding Iteration ==============================
            hist_r_P = []
            for idx_rotor, iter_rotor in enumerate(np.linspace(0, h_wndg_r_1-h_wndg_r_0, n_iters_r)):
                
                h_wndg_r = h_wndg_r_0 + iter_rotor
                
                generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
                generator.apply_coil_sizes_and_lift_factor()
                adapt_yokes()
                
                if generator.P >= gn.Pel_out:
                    
                    # ============================== Find Pel_out Iteration ==============================
                    x = []
                    y = []
                    config = [h_wndg_r, h_wndg_r - (h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r - 2*(h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r*1.1]
                    n = 15
                    for i in range(50):
                        sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Convergence Step {i}         ")
                        sys.stdout.flush()
                        
                        if (generator.P >= (1-detail)*gn.Pel_out and generator.P <= (1+detail)*gn.Pel_out):
                            sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Converged.         ")
                            sys.stdout.flush()

                            generator.HTS_usage()
                            generator.compute_weight()
                            lst_generators.append(copy.deepcopy(generator))
                            break
                        j = min(i, n)
                        fac = 5 * np.abs(gn.Pel_out - generator.P)/gn.Pel_out  
                        if generator.P > gn.Pel_out:
                            h_wndg_r= 0.199*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[1]) + 0.1*(((n-j)/n)*config[1]+(j/n)*h_wndg_r) + 0.2*(((n-j)/n)*config[2]+(j/n)*((1-fac)*h_wndg_r + fac*config[2])) + 0.001*(np.random.randint(2)*config[2])
                        else:
                            h_wndg_r= 0.196*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[0]) + 0.1*(((n-j)/n)*config[0]+(j/n)*h_wndg_r) + 0.2*(((n-j)/n)*config[0]+(j/n)*((1-fac)*h_wndg_r + fac*config[3])) + 0.004*config[3]
                            
                        generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
                        generator.apply_coil_sizes_and_lift_factor()
                        adapt_yokes()
                        if verbose:
                            x.append(i)
                            y.append(generator.P)
                        
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
                    hist_r_P.append(generator.P)
                    if len(hist_r_P)>1 and hist_r_P[-1] < hist_r_P[-2]:
                        sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: BREAK                                     ")
                        sys.stdout.flush()
                        break
                    else:
                        sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r}                                       ")
                        sys.stdout.flush()
            # ============================== End Rotor Winding Iteration ==============================
            
        # ============================== End Stator Winding Iteration ==============================
        if lst_generators:
            # pick generator with least HTS length
            mini = min([geno.HTS_length for geno in lst_generators])
            idx = [geno.HTS_length for geno in lst_generators].index(mini)
            best = lst_generators[idx]    
            
            lst_best_generators.append(best)
        
        total_runs += generator.runs
    # ============================== Pole Pair Iteration ==============================
    job_finish_time = time.time()
    job_runtime = job_finish_time - job_init_time
    sys.stdout.write(f"\rJob finished with {total_runs} solved models in {np.round(job_runtime,3)} seconds." + " "*80)
    sys.stdout.flush()
    # print(f"\n{total_runs = }")

#%% ------ plot HTS length over pole pairs ------
lst_pole_pairs = [gen.p for gen in lst_best_generators]
lst_HTS_length = [gen.HTS_length for gen in lst_best_generators]
lst_weight = [gen.weight for gen in lst_best_generators]

if 1:
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
    
    # import tikzplotlib as tkz
    
    # tkz.clean_figure()
    # tkz.save("aus_LandW_over_p.tex")

    # plt.savefig(fname = "HTSlengthOverPolePairCount.png")
    # plt.show()

#%% ----- select generator and finish setup for export ------
gen = lst_best_generators[20]

if 1:
    
    print(f"{gen.l_e = }, {gen.k_fill_s = }, {gen.k_fill_r = }")
    
    N_s = int(np.round(gen.guess_Ns()))
    N_c_div_a = int(np.round(gen.guess_N_c_div_a(N_s = N_s)))
    
    print(f"\n{N_s = } \n{N_c_div_a = }")
    
    I_s = gen.P / (np.sqrt(3) * gen.gn.U_LL_N)
    ampere_turns_stator = N_c_div_a * I_s
    
    J_s_amplitude = np.sqrt(2) * ampere_turns_stator / gen.coil.A_sc
    gen.N_s = N_s
    gen.J_s_amplitude = J_s_amplitude
    
    print(f"\n{I_s = } [A] \n{ampere_turns_stator = } [A] \n{J_s_amplitude = } [A/m2] \n{J_s_amplitude*1e-6 = }[A_mm2]")

#%% ------ export to comsol -----

save_params("from_python", gen) 

