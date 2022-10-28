# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi

from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot

from data import Generator, FieldWinding, StatorWinding



class n7_Model():
    
    def __init__(self, p: int, l_e: float, r_so: float, 
                 gn: Generator, fw: FieldWinding, sw: StatorWinding,
                 k_fill_s: float = 0.5, k_fill_r: float = 0.5,
                 B_yoke_max: float = 6.0):
        
        
        self.p = p
        self.l_e = l_e
        
        self.dims = n7_Dimensions(r_so)
        
        self.gn = gn
        self.fw = fw
        self.sw = sw
        
        self.J_e_s = sw.J_e
        self.J_e_r = fw.J_e
        
        self.k_fill_s = k_fill_s
        self.k_fill_r = k_fill_r
        self.B_yoke_max = B_yoke_max

        self.mdl = Model(self.p, self.l_e)
        self.plt = None
        self.res = None
        
        self.h_yoke_s = 0.0
        self.h_yoke_r = 0.0
        self.h_wndg_s = 0.0
        self.h_wndg_r = 0.0
        
        self.K_s = 0.0
        self.K_r = 0.0
        
        self.B_s_c = 0.0
        self.B_r_c = 0.0
        
        self.M = 0.0
        self.P = 0.0
        
    
              
    def update_dimensions(self, h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r):
        self.h_yoke_s = h_yoke_s
        self.h_yoke_r = h_yoke_r
        self.h_wndg_s = h_wndg_s
        self.h_wndg_r = h_wndg_r
        self.dims.update_dimensions(h_yoke_s, h_yoke_r, 
                                    h_wndg_s, h_wndg_r, 
                                    self.gn.delta_mag)
            
    
    def update_model(self, K_s, K_r, **kwargs):
        
        self.K_s = K_s
        self.K_r = K_r
                
        alpha_r = kwargs.get("alpha_r", pi/2)
        alpha_s = kwargs.get("alpha_s", 0.0)
        
        ks_d = kwargs.get("ks_d", 1)
        ks_p = kwargs.get("ks_p", 1)
        kr_b = kwargs.get("kr_b", 1)
        
        this_mdl = Model(self.p, self.l_e)
        
        this_mdl.add_layer(AirLayer(r=self.dims.r_ri))
        this_mdl.add_layer(MagneticLayer(r=self.dims.r_ro, 
                                    mu_r=self.gn.mu_r_yoke))
   
        this_mdl.add_layer(CurrentLoading(K=K_r*np.sqrt(2)*kr_b, 
                                    r=self.dims.r_rw, 
                                    alpha=alpha_r,
                                    mu_r=1.
                                    ))

        this_mdl.add_layer(CurrentLoading(K=K_s*np.sqrt(2)*ks_d*ks_p,
                                    r=self.dims.r_sw,
                                    alpha=alpha_s,
                                    mu_r=1.
                                    ))

        this_mdl.add_layer(AirLayer(r=self.dims.r_si))
        this_mdl.add_layer(MagneticLayer(r=self.dims.r_so, 
                                    mu_r=self.gn.mu_r_yoke))
        
        this_mdl.build()
        this_mdl.solve()
        this_mdl.total_torque()
        
        this_plt = PlanePlot(this_mdl, fgsz=100)
        
        flux_data = this_mdl.get_B_data(r=np.array([self.dims.r_rw, 
                                                    self.dims.r_sw]), 
                                        t=np.linspace(0, np.pi/self.p, 400)
                                        )
        
        self.B_r_c = np.max(np.abs(np.sin(np.pi/6) * flux_data.Br[:,0] + \
                              np.cos(np.pi/6) * flux_data.Bt[:, 0])
                            )
        self.B_s_c = np.max(np.abs(np.sin(np.pi/6) * flux_data.Br[:,1] + \
                              np.cos(np.pi/6) * flux_data.Bt[:, 1])
                            )
        
        
        self.M = this_mdl.Mpos
        self.P = self.gn.w_syn * self.M
        
        
        self.mdl = this_mdl
        self.plt = this_plt
        
        
    def show_results(self, header="n7 - Results"):
        print("---", header, "---")
        print(f"M = {np.round(self.M * 1e-6, 2)} [MNm]")
        print(f"P = {np.round(self.P * 1e-6, 2)} [MW]")
        print(f"K_r = {np.round(self.K_r * 1e-3, 2)} [kA/m]")
        print(f"K_s = {np.round(self.K_s * 1e-3, 2)} [kA/m]")
        print(f"h_wndg_r = {np.round(self.h_wndg_r, 4)} [m]")
        print(f"h_wndg_s = {np.round(self.h_wndg_s, 4)} [m]")
        print(f"h_yoke_r = {np.round(self.h_yoke_r, 4)} [m]")
        print(f"h_yoke_s = {np.round(self.h_yoke_s, 4)} [m]")
        print(f"B_r_c = {np.round(self.B_r_c, 4)} [T]")
        print(f"B_s_c = {np.round(self.B_s_c, 4)} [T]")
        print(f"J_e_r = {np.round(self.J_e_r * 1e-6, 4)} [A/mm2]")
        print(f"J_e_s = {np.round(self.J_e_s * 1e-6, 4)} [A/mm2]")
        

    
    
class n7_Dimensions():
    
    def __init__(self, r_so):
        self.r_so = r_so
        self.r_si = 0.0
        self.r_sw = 0.0
        self.r_ag = 0.0
        self.r_rw = 0.0
        self.r_ro = 0.0
        self.r_ri = 0.0
        
        
    def update_dimensions(self, h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r, delta_mag):
        self.r_si = self.r_so - h_yoke_s
        self.r_sw = self.r_si - h_wndg_s/2
        self.r_ag = self.r_sw - h_wndg_s/2 - delta_mag/2
        self.r_rw = self.r_ag - delta_mag/2 - h_wndg_r/2
        self.r_ro = self.r_rw - h_wndg_r/2
        self.r_ri = self.r_ro - h_yoke_r
        
     
    def show(self, header="n7 - Dimensions"):
        print("---", header, "---")
        print(f"r_so = {np.round(self.r_so, 4)} [m]")
        print(f"r_si = {np.round(self.r_si, 4)} [m]")
        print(f"r_sw = {np.round(self.r_sw, 4)} [m]")
        print(f"r_ag = {np.round(self.r_ag, 4)} [m]")
        print(f"r_rw = {np.round(self.r_rw, 4)} [m]")
        print(f"r_ro = {np.round(self.r_ro, 4)} [m]")
        print(f"r_ri = {np.round(self.r_ri, 4)} [m]")
        
    
    