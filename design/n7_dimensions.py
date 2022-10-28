# -*- coding: utf-8 -*-
import numpy as np

class n7_Dimensions():
    
    def __init__(self, r_so, r_si, r_sA, r_rF, r_ro, r_ri):
        self.r_so = r_so
        self.r_si = r_si
        self.r_sA = r_sA
        self.r_rF = r_rF
        self.r_ro = r_ro
        self.r_ri = r_ri
     
    def show(self, header="n7_Dimensions"):
        print("---", header, "---")
        print(f"r_so = {np.round(self.r_so, 4)} [m]")
        print(f"r_si = {np.round(self.r_si, 4)} [m]")
        print(f"r_sA = {np.round(self.r_sA, 4)} [m]")
        print(f"r_rF = {np.round(self.r_rF, 4)} [m]")
        print(f"r_ro = {np.round(self.r_ro, 4)} [m]")
        print(f"r_ri = {np.round(self.r_ri, 4)} [m]")