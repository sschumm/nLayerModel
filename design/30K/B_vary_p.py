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
from data import StatorWinding_30K as sw
from data import StatorWinding_77K as tw
from data.export import save_params

from design import n7_Model


# --- ecoswing ---
p = 30
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



#%% ------ iterating over p ------
if 0:
    P_goal = 3e6
    
    
    job_init_time = time.time()
    # r_so = 1 # 3.75
    # l_e = sw.l_e(sw, r_so) # 0.267
    l_e = 0.2
    r_so = sw.r_so(sw, l_e) # 1.812
    total_runs = 0
    
    # ============================== Pole Pair Iteration ==============================
    lst_best_generators = []
    for p in range(4, 61, 1):

        generator = n7_Model(p, l_e, r_so, 
                             gn=gn, fw=fw, sw=sw,
                             k_fill_s=0.3,
                             k_fill_r=0.3,
                             B_yoke_max=4)

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

#%%
P_goal = 3e6
l_e = 0.2
r_so = 1.5 # sw.r_so(sw, l_e) # 1.812
p=20

generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     B_yoke_max=4)

# --- init generator dimension via yoke and winding heights ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.3, 1.3]]
generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)

def get_k_fill_s(r_si, h_sw, p):
    nom = (pi*(r_si-h_sw)/(3*p)-gn.w_pole_frame-gn.r_bend_max)*(h_sw-2*gn.h_pole_frame)*3*p
    den = h_sw * (r_si - h_sw) * pi
    return nom / den

def get_k_fill_r(r_ro, h_rw, p):
    nom = 2*p*(((pi * r_ro)/(2 * p)) - gn.w_pole_frame - gn.r_bend_max)*(h_rw - 2*gn.h_pole_frame)
    den = h_rw * pi * r_ro
    return nom / den

generator.k_fill_s = get_k_fill_s(generator.dims.r_si, generator.h_wndg_s, generator.p)
generator.k_fill_r = get_k_fill_r(generator.dims.r_ro, generator.h_wndg_r, generator.p)
    
# --- compute initial current loadings and with it the initial model ---
K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r

generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.866)
# generator.keep_const_kr_b(kr_b = 0.65)
generator.apply_coil_sizes_and_lift_factor(verbose = False)
adapt_yokes()

# generator.show_results()

#%%
# ============================== Stator Winding Iteration ==============================    
lst_generators = []

h_wndg_s_0 = 2 * gn.h_pole_frame * 1.1
h_wndg_s_1 = h_wndg_s_0 * 2

h_wndg_r_0 = 2 * gn.h_pole_frame * 1.1
h_wndg_r_1 = h_wndg_r_0 * 2

n_iters_s = 9
n_iters_r = 20 

detail = 0.1
verbose = True

for idx_stator, iter_stator in enumerate(np.linspace(0, h_wndg_s_1-h_wndg_s_0, n_iters_s)):
    
    h_wndg_s = h_wndg_s_0 + iter_stator
    generator.update_dimensions(h_wndg_s = h_wndg_s, keep_const_kr_b=False)
    generator.k_fill_s = get_k_fill_s(generator.dims.r_si, generator.h_wndg_s, generator.p)
    
    
    # ============================== Rotor Winding Iteration ==============================
    hist_r_P = []
    for idx_rotor, iter_rotor in enumerate(np.linspace(0, h_wndg_r_1-h_wndg_r_0, n_iters_r)):
        
        h_wndg_r = h_wndg_r_0 + iter_rotor
        generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=False)
        
        generator.k_fill_r = get_k_fill_r(generator.dims.r_ro, generator.h_wndg_r, generator.p)
        generator.apply_coil_sizes_and_lift_factor()
        adapt_yokes()
        
        if generator.P >= P_goal:
            
            # ============================== Find Pel_out Iteration ==============================
            x = []
            y = []
            config = [h_wndg_r, h_wndg_r - (h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r - 2*(h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r*1.1]
            n = 15
            for i in range(50):
                sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Convergence Step {i}         ")
                sys.stdout.flush()
                
                if (generator.P >= (1-detail)*P_goal and generator.P <= (1+detail)*P_goal):
                    sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Converged.         ")
                    sys.stdout.flush()

                    generator.HTS_usage()
                    generator.compute_weight()
                    lst_generators.append(copy.deepcopy(generator))
                    break
                j = min(i, n)
                fac = 5 * np.abs(P_goal - generator.P)/P_goal  
                if generator.P > P_goal:
                    h_wndg_r= 0.199*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[1]) + 0.1*(((n-j)/n)*config[1]+(j/n)*h_wndg_r) + 0.2*(((n-j)/n)*config[2]+(j/n)*((1-fac)*h_wndg_r + fac*config[2])) + 0.001*(np.random.randint(2)*config[2])
                else:
                    h_wndg_r= 0.196*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[0]) + 0.1*(((n-j)/n)*config[0]+(j/n)*h_wndg_r) + 0.2*(((n-j)/n)*config[0]+(j/n)*((1-fac)*h_wndg_r + fac*config[3])) + 0.004*config[3]
                    
                generator.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=False)
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
    #     if lst_generators:
    #         # pick generator with least HTS length
    #         mini = min([geno.HTS_length for geno in lst_generators])
    #         idx = [geno.HTS_length for geno in lst_generators].index(mini)
    #         best = lst_generators[idx]    
            
    #         lst_best_generators.append(best)
        
    #     total_runs += generator.runs
    # # ============================== Pole Pair Iteration ==============================
    # job_finish_time = time.time()
    # job_runtime = job_finish_time - job_init_time
    # sys.stdout.write(f"\rJob finished with {total_runs} solved models in {np.round(job_runtime,3)} seconds." + " "*80)
    # sys.stdout.flush()