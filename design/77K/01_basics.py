# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplot

from scipy.constants import pi

from analytics import yoke_height

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw

from design import n7_Model

# --- init generator with main parameter ---
p = 18
l_e = 0.3
r_so = sw.r_so(sw, l_e)

generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     k_fill_s=0.3,
                     k_fill_r=0.3,
                     B_yoke_max=6.0)


# --- init generator dimension via yoke and winding heights ---
h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = 0.5 * gn.delta_mag, 0.5 * gn.delta_mag

generator.init_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)


# --- compute initial current loadings and with it the initial model ---
K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r

generator.update_model_by_K(K_s = K_s, K_r = K_r)
generator.show_results("Results for the initial model")
if 0:generator.fast_flux()


#%% --- apply the lift factor to be self-consistent ---
generator.apply_lift_factor(verbose = False)
generator.show_results("Results for initial dimensions with applied lift factor")
if 0: generator.fast_flux()

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
if 0: generator.fast_flux()

#%% --- vary h_wndg_r ---
if 1:
    increasing_h_by_factor=0.2
    iter_P_h_wndg_r, iter_h_wndg_r = [generator.P], [generator.h_wndg_r]
    
    for idx_h_wndg_r in range(100):
        h_wndg_r = iter_h_wndg_r[-1] + iter_h_wndg_r[0] * increasing_h_by_factor
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
    h_wndg_r = 0.1
    generator.update_dimensions(h_wndg_r=h_wndg_r)
    generator.apply_lift_factor(verbose=True)
    generator.show_results("Results for adapted h_wndg_r")
    if 0: generator.fast_flux()


#%% --- vary h_wndg_s ---
if 1:
    # --- interating ---
    increasing_h_by_factor=0.1
    iter_P_h_wndg_s, iter_h_wndg_s = [generator.P], [generator.h_wndg_s]
    
    for idx_h_wndg_s in range(50):
        h_wndg_s = iter_h_wndg_s[-1] + iter_h_wndg_s[0] * increasing_h_by_factor
        generator.update_dimensions(h_wndg_s=h_wndg_s)
        generator.apply_lift_factor()
        
        iter_h_wndg_s.append(h_wndg_s)
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
    h_wndg_s = 0.05 
    generator.update_dimensions(h_wndg_s=h_wndg_s)
    generator.apply_lift_factor(verbose=True)
    generator.show_results("Results for adapted h_wndg_s")
    if 0: generator.fast_flux()
        



#%% --- check coil dimensions ---
generator.coil_shapes()
generator.coil.show()














