# -*- coding: utf-8 -*-
import numpy as np

class n7_Dimensions():
    
    def __init__(self, r_so, r_si, r_sw, r_ag, r_rw, r_ro, r_ri):
        self.r_so = r_so
        self.r_si = r_si
        self.r_sw = r_sw
        self.r_ag = r_ag
        self.r_rw = r_rw
        self.r_ro = r_ro
        self.r_ri = r_ri
        
        
    def update_dimensions(self, h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r, delta_mag):
        self.r_si = self.r_so - h_yoke_s
        self.r_sw = self.r_si - h_wndg_s/2
        self.r_ag = self.r_sw - h_wndg_s/2 - delta_mag/2
        self.r_rw = self.r_ag - delta_mag/2 - h_wndg_r/2
        self.r_ro = self.r_rw - h_wndg_r/2
        self.r_ri = self.r_ro - h_yoke_r
        
     
    def show(self, header="n7_Dimensions"):
        print("---", header, "---")
        print(f"r_so = {np.round(self.r_so, 4)} [m]")
        print(f"r_si = {np.round(self.r_si, 4)} [m]")
        print(f"r_sw = {np.round(self.r_sw, 4)} [m]")
        print(f"r_ag = {np.round(self.r_ag, 4)} [m]")
        print(f"r_rw = {np.round(self.r_rw, 4)} [m]")
        print(f"r_ro = {np.round(self.r_ro, 4)} [m]")
        print(f"r_ri = {np.round(self.r_ri, 4)} [m]")