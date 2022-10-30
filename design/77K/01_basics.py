# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplot

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw
from data.export import save_params

from design import n7_Model

# --- init generator with main parameter ---
p = 25
l_e = 0.2
r_so = sw.r_so(sw, l_e)

generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     k_fill_s=0.3,
                     k_fill_r=0.2,
                     B_yoke_max=6.0)


# --- init generator dimension via yoke and winding heights ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = [factor * (gn.h_pole_frame*2) for factor in [1.5, 1.5]]

generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)


# --- compute initial current loadings and with it the initial model ---
K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r

generator.update_model_by_K(K_s = K_s, K_r = K_r, ks_d=0.9577)
generator.show_results("Results for the initial model")
# generator.fast_flux()


#%% --- apply the lift factor to be self-consistent ---
generator.apply_lift_factor(verbose = False)
generator.show_results("Results for initial dimensions with applied lift factor")
# generator.fast_flux()

#%% --- adapt yokes ---
# --- adapt stator yoke ---
h_yoke_s = yoke_height(generator.p, generator.dims.r_si, 
                       generator.B_s, generator.B_yoke_max)

# --- adapt rotor yoke ---
h_yoke_r = yoke_height(generator.p, generator.dims.r_ro, 
                       generator.B_r, generator.B_yoke_max)

generator.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r)
generator.apply_lift_factor(verbose = False)
generator.show_results("Results with adapted yokes")
generator.fast_flux()

#%% --- vary h_wndg_s ---
if 1:
    increasing_h_by_factor=0.01
    iter_P_h_wndg_s, iter_h_wndg_s = [generator.P], [generator.h_wndg_s]
    
    for idx_h_wndg_s in range(50):
        h_wndg_s = iter_h_wndg_s[-1] + iter_h_wndg_s[0] * increasing_h_by_factor
        if h_wndg_s <= 0 or h_wndg_s >= generator.dims.r_ro/2:
            break
        generator.update_dimensions(h_wndg_s=h_wndg_s)
        generator.apply_lift_factor()
        
        iter_h_wndg_s.append(generator.h_wndg_s)
        iter_P_h_wndg_s.append(generator.P)
        
    
    fig, axs = pyplot.subplots(1)
    fig.tight_layout() 
    fig.dpi=500
    fig.set_figheight(4)
    fig.set_figwidth(5)
    
    ax1 = axs
    ax1.grid()
    ax1.plot(iter_h_wndg_s, [i*1e-6 for i in iter_P_h_wndg_s])
    ax1.set_xlabel("h_wndg_s [m]")
    ax1.set_ylabel('P [MW]', color="b")

#%% --- select h_wndg_s ---
if 1:
    h_wndg_s = 0.029
    generator.update_dimensions(h_wndg_s=h_wndg_s)
    generator.apply_lift_factor(verbose=True)
    h_yoke_s = yoke_height(generator.p, generator.dims.r_si, 
                           generator.B_s, generator.B_yoke_max)
    generator.update_dimensions(h_yoke_s=h_yoke_s)
    generator.apply_lift_factor()
    generator.show_results("Results for adapted h_wndg_s")
    # generator.fast_flux()

#%% --- vary h_wndg_r ---
if 1:
    increasing_h_by_factor=0.05
    iter_P_h_wndg_r, iter_h_wndg_r = [generator.P], [generator.h_wndg_r*0.8]
    
    for idx_h_wndg_r in range(100):
        h_wndg_r = iter_h_wndg_r[-1] + iter_h_wndg_r[0] * increasing_h_by_factor
        if h_wndg_r <= 0 or h_wndg_r >= generator.dims.r_ro/2:
            break
        generator.update_dimensions(h_wndg_r=h_wndg_r)
        generator.apply_lift_factor()
        
        iter_h_wndg_r.append(generator.h_wndg_r)
        iter_P_h_wndg_r.append(generator.P)
        
    
    fig, axs = pyplot.subplots(1)
    fig.tight_layout() 
    fig.dpi=500
    fig.set_figheight(4)
    fig.set_figwidth(5)
    
    ax1 = axs
    ax1.grid()
    ax1.plot(iter_h_wndg_r, [i*1e-6 for i in iter_P_h_wndg_r])
    ax1.set_xlabel("h_wndg_r [m]")
    ax1.set_ylabel('P [MW]', color="b")

#%% --- select h_wndg_r ---
if 1:
    h_wndg_r = 0.03
    generator.update_dimensions(h_wndg_r=h_wndg_r)
    generator.apply_lift_factor(verbose=True)
    h_yoke_r = yoke_height(generator.p, generator.dims.r_ro, 
                           generator.B_r, generator.B_yoke_max)
    generator.update_dimensions(h_yoke_r=h_yoke_r)
    generator.apply_lift_factor()
    generator.show_results("Results for adapted h_wndg_r")
    # generator.fast_flux()

#%% --- check coil dimensions ---
generator.coil_shapes()
generator.coil.show()


#%%
save_params("from_python", generator)










