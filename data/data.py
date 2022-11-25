# -*- coding: utf-8 -*-
import numpy as np
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
    Ic_spec = 650      # [A]
    A_tape = (12*0.22) * 1e-6 # [m**2]
    h_pole_frame = 0.009 # [m] == 9 [mm] due to h_pressure_plate 18 [mm] / 2
    w_pole_frame = 0.009 # [m] same reasoning as above
    r_bend_max = 0.02 # [m] == 20 [mm] 
    HTS_weight_length_corr = 48 / 500 # [kg/m], guess from Bergen
    Iron_density = 7760 # [kg/m3], from isovac data sheet
    G10CR_density = 1850 # [kg/m3], from (insert source)
    HTS_density = 6423 # [kg/m3], from (insert source): 1/3 * (cu + hastelloy + kapton) = 1/3 * (8960 + 8890 + 1420)
    thickness_ratio = 0.25 # = l_e / d_so
    
    
class FieldWinding():
    
    T_HTS = 30 # [K] - LNe - Liquid Neon
    J_e = 154.5e6 # [A/m2] == 270/1.7 ~ 158.8 [A/mm2] - engineering current density
    
class StatorWinding():
    d2so_L = 5
    
    def l_e(self, r_so):
        return self.d2so_L/ (2*r_so)**2
    
    def r_so(self, l_e):
        return np.sqrt(self.d2so_L/l_e)/2
    
    
class StatorWinding_Cu(StatorWinding):
            
    # T_HTS = 77 # [K] - LN2 - Liquid Nitrogen
    J_e = 3.8e6 # [A/m2] == 3.8 [A/mm2] - current density
    C_e = 2262600.0 # [VAs/m3] == 37.71 [kVA min/m3] - Essons Number - electromagnetic utilization
    d2so_L = 22.28 # [m3]
    V_ges = pi/4 * d2so_L # [m3]

    
class StatorWinding_77K(StatorWinding):
            
    T_HTS = 77 # [K] - LN2 - Liquid Nitrogen
    J_e = 11.1e6 # [A/m2] == 20/1.7 ~ 11.8 [A/mm2] - engineering current density
    C_e = 55*60*1e3 # [VAs/m3] == 0.5 * 117 [kVA min/m3] - Essons Number - electromagnetic utilization
    d2so_L = 15 # [m3]
    
    
class StatorWinding_30K(StatorWinding):
        
    T_HTS = 30 # [K] - LNe - Liquid Neon
    J_e = 154.5e6 # [A/m2] == 270/1.7 ~ 158.8 [A/mm2] - engineering current density
    C_e = 766*60*1e3 # [VAs/m3] == 0.5 * 1575.9 [kVA min/m3] - Essons Number - electromagnetic utilization
    d2so_L = 1 # [m3]