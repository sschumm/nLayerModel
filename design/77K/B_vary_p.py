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


# ============================== Initialize Iteration ==============================

job_init_time = time.time()
total_runs = 0

p = 20

# ============================== Pole Pair Iteration ==============================
# lst_best_generators = []
# for p in range(8, 61, 1):

gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw)
h_pf = gn.h_pole_frame
gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                    h_wndg_s=h_pf * 2.1, h_wndg_r=h_pf * 2.1)
gen.ks_d= 0.866

#     # --- init generator dimension via yoke and winding heights ---
#     h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
#     h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.05, 1.05]]
#     generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)
        
#     # --- compute initial current loadings and with it the initial model ---
#     K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
#     K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r
#     generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.866)
#     generator.keep_const_kr_b(kr_b = 0.65)
#     generator.apply_coil_sizes_and_lift_factor(verbose = False)
#     adapt_yokes()


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
        
        # if gen.P >= gn.Pel_out:
            
        #     # ============================== Find Pel_out Iteration ==============================
        #     x = []
        #     y = []
        #     config = [h_wndg_r, h_wndg_r - (h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r - 2*(h_wndg_r_1-h_wndg_r_0)/n_iters_r, h_wndg_r*1.1]
        #     n = 15
        #     for i in range(50):
        #         sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Convergence Step {i}         ")
        #         sys.stdout.flush()
                
        #         if (gen.P >= (1-detail)*gn.Pel_out and gen.P <= (1+detail)*gn.Pel_out):
        #             sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r} - Converged.         ")
        #             sys.stdout.flush()

        #             gen.HTS_usage()
        #             gen.compute_weight()
        #             lst_gens.append(copy.deepcopy(gen))
        #             break
        #         j = min(i, n)
        #         fac = 5 * np.abs(gn.Pel_out - gen.P)/gn.Pel_out  
        #         if gen.P > gn.Pel_out:
        #             h_wndg_r= 0.199*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[1]) + 0.1*(((n-j)/n)*config[1]+(j/n)*h_wndg_r) + 0.2*(((n-j)/n)*config[2]+(j/n)*((1-fac)*h_wndg_r + fac*config[2])) + 0.001*(np.random.randint(2)*config[2])
        #         else:
        #             h_wndg_r= 0.196*h_wndg_r + 0.5*((1-fac)*h_wndg_r + fac*config[0]) + 0.1*(((n-j)/n)*config[0]+(j/n)*h_wndg_r) + 0.2*(((n-j)/n)*config[0]+(j/n)*((1-fac)*h_wndg_r + fac*config[3])) + 0.004*config[3]
                    
        #         gen.update_dimensions(h_wndg_r = h_wndg_r, keep_const_kr_b=True)
        #         gen.apply_coil_sizes_and_lift_factor()
        #         adapt_yokes()
        #         if verbose:
        #             x.append(i)
        #             y.append(gen.P)
                
        #     if verbose:
        #         fig = plt.figure(dpi=300, figsize=(6,7))
        #         ax1 = plt.subplot()
        #         ax1.plot(x,y)
        #         ax1.scatter(x,y)
        #         plt.grid()
        #         plt.show()
        #     break
            # ============================== End Pel_out Iteration ==============================
        # else:
        #     hist_r_P.append(gen.P)
        #     if len(hist_r_P)>1 and hist_r_P[-1] < hist_r_P[-2]:
        #         sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: BREAK                                     ")
        #         sys.stdout.flush()
        #         break
        #     else:
        #         sys.stdout.write(f"\rComputing Pole Pair: {p} - Stator Iteration: {idx_stator+1} / {n_iters_s} - Rotor Iterations: {idx_rotor+1} / {n_iters_r}                                       ")
        #         sys.stdout.flush()
        lst_rotor_iter.append(copy.deepcopy(gen))
    # ============================== End Rotor Winding Iteration ==============================
    lst_stator_iter.append(lst_rotor_iter)
# ============================== End Stator Winding Iteration ==============================
    # if lst_generators:
    #     # pick generator with least HTS length
    #     mini = min([geno.HTS_length for geno in lst_generators])
    #     idx = [geno.HTS_length for geno in lst_generators].index(mini)
    #     best = lst_generators[idx]    
        
    #     lst_best_generators.append(best)
    
# total_runs += generator.runs
# ============================== Pole Pair Iteration ==============================
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