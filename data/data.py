# -*- coding: utf-8 -*-

from scipy.constants import pi

class Generator():
    
    m = 3
    Pel_out = 7 * 1e6  # == 7 [MW]
    n_syn = 8.33       # [rpm]
    w_syn = 0.8723     # [rad/s]        
    U_LL_N = 690       # [V] - nominal line to line voltage
    delta_mag = 0.04   # [m] == 40 [mm]
    mu_r_yoke = 1e5    # ferro magnetic
    mu_r_core = 1      # G-10 CR (Fiberglass Epoxy) or Stainless Steel
    
    
class FieldWinding():
    
    T_HTS = 30 # [K] - LNe - Liquid Neon
    J_e = 11.8e6 # [A/m2] == 20/1.7 ~ 11.8 [A/mm2] - engineering current density
    
    
class StatorWinding_Cu():
    
    # T_HTS = 77 # [K] - LN2 - Liquid Nitrogen
    J_e = 3.8e6 # [A/m2] == 3.8 [A/mm2] - current density
    C_e = 2262600.0 # [VAs/m3] == 37.71 [kVA min/m3] - Essons Number - electromagnetic utilization
    d2so_L = 22.28 # [m3]
    V_ges = pi/4 * d2so_L # [m3]

    
class StatorWinding_77K():
    
    T_HTS = 77 # [K] - LN2 - Liquid Nitrogen
    J_e = 11.8e6 # [A/m2] == 20/1.7 ~ 11.8 [A/mm2] - engineering current density
    C_e = 0.5*117*60*1e3 # [VAs/m3] == 0.5 * 117 [kVA min/m3] - Essons Number - electromagnetic utilization
    d2so_L = 11.28 # [m3]
    
    
class StatorWinding_30K():
    
    T_HTS = 30 # [K] - LNe - Liquid Neon
    J_e = 158.8e6 # [A/m2] == 270/1.7 ~ 158.8 [A/mm2] - engineering current density
    C_e = 0.5*1575.9*60*1e3 # [VAs/m3] == 0.5 * 1575.9 [kVA min/m3] - Essons Number - electromagnetic utilization
    d2so_L = 1 # [m3]