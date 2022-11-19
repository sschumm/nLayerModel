# -*- coding: utf-8 -*-
import sys
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
p = 20
r_so = 3.75
l_e = sw.l_e(sw, r_so) # 0.267


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
if 1:
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
            











