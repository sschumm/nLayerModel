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
from data.export import save_params

from design import n7_Model


#%% Ausgehend von Schlankheitsgrad

d2l = 2 # m3
l_d = 0.25 # -

d_so = (d2l / l_d)**(1/3)
r_so = d_so / 2
l_e = l_d * d_so

print(f"{r_so = }")
print(f"{l_e = }")

#%% Schlankheitsgrad berechnen

d2l = 2 # m3

r_so = 1
d_so = 2* r_so
l_e = d2l / (d_so)**2
l_d = l_e / d_so

print(f"{l_d = }")

#%% Init generator Model
p = 20
gen = n7_Model(p, l_e, r_so, gn=gn, fw=fw, sw=sw)
h_pf = gn.h_pole_frame
gen.init_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, 
                    h_wndg_s=h_pf * 2.1, 
                    h_wndg_r=h_pf * 2.1
                    )

def get_k_fill_s(r_si, h_sw, p):
    nom = (pi*(r_si-h_sw)/(3*p)-gn.w_pole_frame-gn.r_bend_max)*(h_sw-2*gn.h_pole_frame)*3*p
    den = h_sw * (r_si - h_sw) * pi
    return nom / den

def get_k_fill_r(r_ro, h_rw, p):
    nom = 2*p*(((pi * r_ro)/(2 * p)) - gn.w_pole_frame - gn.r_bend_max)*(h_rw - 2*gn.h_pole_frame)
    den = h_rw * pi * r_ro
    return nom / den

gen.k_fill_s = get_k_fill_s(r_si= gen.dims.r_si, h_sw= gen.h_wndg_s, p= gen.p)
gen.k_fill_r = get_k_fill_r(r_ro= gen.dims.r_ro, h_rw= gen.h_wndg_r, p= gen.p)

gen.update_model_by_K(K_s= gen.k_fill_s * gen.J_e_s * gen.h_wndg_s, 
                      K_r= gen.k_fill_r * gen.J_e_r * gen.h_wndg_r, 
                      ks_d=0.866)

# gen.show_results("initialization")
# gen.fast_flux()

#%% Iterate over h_sw and h_rw
lst_stator = []
h_sw0 = h_pf*2.1
h_rw0 = h_pf*2.1
for iter_stator in np.linspace(0,3,9):
    gen.update_dimensions(h_wndg_s= h_sw0 + h_sw0*iter_stator)
    gen.k_fill_s = get_k_fill_s(r_si= gen.dims.r_si, h_sw= gen.h_wndg_s, p= gen.p)
    
    lst_rotor = []
    for iter_rotor in np.linspace(0,3,20):
        gen.update_dimensions(h_wndg_r= h_rw0 + h_rw0*iter_rotor)
        gen.k_fill_r = get_k_fill_r(r_ro= gen.dims.r_ro, h_rw= gen.h_wndg_r, p= gen.p)
        gen.update_model_by_K(K_s= gen.k_fill_s * gen.J_e_s * gen.h_wndg_s, 
                              K_r= gen.k_fill_r * gen.J_e_r * gen.h_wndg_r, 
                              ks_d=0.866)
        lst_rotor.append(copy.deepcopy(gen))
    lst_stator.append(lst_rotor)
        
#%% Plot iteration of h_sw and h_rw
fig = plt.figure(dpi=300, figsize=(6,7))
ax1 = plt.subplot()
plt.grid()

for lst in lst_stator:
    ax1.plot([g.h_wndg_r for g in lst], [g.P for g in lst], label=f"h_sw= {np.round(lst[0].h_wndg_s, 3)}")

ax1.set_xlabel("h_rw")
ax1.set_ylabel("P")
ax1.legend(loc="upper left")
plt.show()

#%% Apply Lift Factor 

gen.update_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_s= 0.03, h_wndg_r= 0.03)
gen.k_fill_s = get_k_fill_s(r_si= gen.dims.r_si, h_sw= gen.h_wndg_s, p= gen.p)
gen.k_fill_r = get_k_fill_r(r_ro= gen.dims.r_ro, h_rw= gen.h_wndg_r, p= gen.p)

gen.update_model_by_K(K_s= gen.k_fill_s * gen.J_e_s * gen.h_wndg_s, 
                      K_r= gen.k_fill_r * gen.J_e_r * gen.h_wndg_r, 
                      ks_d=0.866)

gen.apply_lift_factor()

#%% 
def adapt_yokes():
    # --- adapt stator yoke ---
    h_yoke_s = yoke_height(gen.p, gen.dims.r_si, gen.B_s, gen.B_yoke_max)    
    # --- adapt rotor yoke ---
    h_yoke_r = yoke_height(gen.p, gen.dims.r_ro, gen.B_r, gen.B_yoke_max)

    gen.update_dimensions(h_yoke_s=h_yoke_s, h_yoke_r=h_yoke_r)

adapt_yokes()

gen.k_fill_s = get_k_fill_s(r_si= gen.dims.r_si, h_sw= gen.h_wndg_s, p= gen.p)
gen.k_fill_r = get_k_fill_r(r_ro= gen.dims.r_ro, h_rw= gen.h_wndg_r, p= gen.p)

gen.update_model_by_K(K_s= gen.k_fill_s * gen.J_e_s * gen.h_wndg_s, 
                      K_r= gen.k_fill_r * gen.J_e_r * gen.h_wndg_r, 
                      ks_d=0.866)

gen.apply_lift_factor()


#%% Iterate over h_sw and h_rw with lift factor (and adapt yokes)
lst_stator = []
h_sw0 = h_pf*2.1 #0.3 
h_rw0 = h_pf*2.1 #0.3
for iter_stator in np.linspace(0,3,9):
    gen.update_dimensions(h_wndg_s= h_sw0 + h_sw0*iter_stator)
    gen.k_fill_s = get_k_fill_s(r_si= gen.dims.r_si, h_sw= gen.h_wndg_s, p= gen.p)
    
    lst_rotor = []
    for iter_rotor in np.linspace(0,3,20):
        gen.update_dimensions(h_wndg_r= h_rw0 + h_rw0*iter_rotor)
        gen.k_fill_r = get_k_fill_r(r_ro= gen.dims.r_ro, h_rw= gen.h_wndg_r, p= gen.p)
        gen.update_model_by_K(K_s= gen.k_fill_s * gen.J_e_s * gen.h_wndg_s, 
                              K_r= gen.k_fill_r * gen.J_e_r * gen.h_wndg_r, 
                              ks_d=0.866)
        gen.apply_lift_factor()
        adapt_yokes()
        gen.apply_lift_factor()
        lst_rotor.append(copy.deepcopy(gen))
    lst_stator.append(lst_rotor)

#%% Plot iteration of h_sw and h_rw with lift factor
fig = plt.figure(dpi=300, figsize=(6,7))
ax1 = plt.subplot()
plt.grid()

for lst in lst_stator:
    ax1.plot([g.h_wndg_r for g in lst], [g.P for g in lst], label=f"h_sw= {np.round(lst[0].h_wndg_s, 3)}")

ax1.set_xlabel("h_rw")
ax1.set_ylabel("P")
ax1.legend(loc="upper left")
plt.show()

#%% test configuration

gen.update_dimensions(h_yoke_s=h_pf, h_yoke_r=h_pf, h_wndg_s= 0.047, h_wndg_r= 0.055)
gen.k_fill_s = get_k_fill_s(r_si= gen.dims.r_si, h_sw= gen.h_wndg_s, p= gen.p)
gen.k_fill_r = get_k_fill_r(r_ro= gen.dims.r_ro, h_rw= gen.h_wndg_r, p= gen.p)
gen.update_model_by_K(K_s= gen.k_fill_s * gen.J_e_s * gen.h_wndg_s, 
                      K_r= gen.k_fill_r * gen.J_e_r * gen.h_wndg_r, 
                      ks_d=0.866)
gen.apply_lift_factor()
adapt_yokes()
gen.apply_lift_factor()

gen.show_results()
# gen.fast_flux()

