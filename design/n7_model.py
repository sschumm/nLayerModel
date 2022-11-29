# -*- coding: utf-8 -*-
import sys
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as pyplot
from scipy.constants import pi

from analytics import yoke_height

from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot

from data import Generator, FieldWinding, StatorWinding

from design import get_L_TPL2100



class n7_Model():
    
    def __init__(self, p: int, l_e: float, r_so: float, 
                 gn: Generator, fw: FieldWinding, sw: StatorWinding,
                 k_fill_s: float = 0.5, k_fill_r: float = 0.5,
                 B_yoke_max: float = 6.0):
        
        
        self.p = p
        self.l_e = l_e
        self.q = 0.5
        
        self.dims = n7_Dimensions(r_so)
        
        self.gn = gn
        self.fw = fw
        self.sw = sw
        
        self.f = self.p * self.gn.n_syn / 60
        self.J_e_spec = gn.Ic_spec / (1.7 * gn.A_tape)
        
        self.J_e_s = sw.J_e
        self.J_e_r = fw.J_e
        
        self.k_fill_s = k_fill_s
        self.k_fill_r = k_fill_r
        self.B_yoke_max = B_yoke_max

        self.mdl = Model(self.p, self.l_e)
        self.plt = None
        
        self.h_yoke_s = 0.0
        self.h_yoke_r = 0.0
        self.h_wndg_s = 0.0
        self.h_wndg_r = 0.0
        
        self.K_s = 0.0
        self.K_r = 0.0
        
        self.ks_d = 1
        self.ks_p = 1
        self.kr_b = 1
        
        self.B_s = 0.0
        self.B_r = 0.0
        self.B_s_c = 0.0
        self.B_r_c = 0.0
        
        self.M = 0.0
        self.P = 0.0
        
        self.coil = None
        self.N_s = 0
        self.J_s_amplitude = 0.0
        
        self.runs = 0
        self.HTS_volume = 0
        self.HTS_volume_rotor = 0
        self.HTS_volume_stator = 0
        self.HTS_length = 0
        self.HTS_length_rotor = 0
        self.HTS_length_stator = 0
        self.weight = 0
        self.weight_rotor_yoke = 0
        self.weight_rotor_wndg = 0
        self.weight_stator_yoke = 0
        self.weight_stator_wndg = 0
        
        self.reference_weight = 0
        self.weight_rotor_kryo = 0
        self.weight_stator_kryo = 0
        
    
              
    def init_dimensions(self, h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r, **kwargs):
        self.h_yoke_s = h_yoke_s
        self.h_yoke_r = h_yoke_r
        self.h_wndg_s = h_wndg_s
        self.h_wndg_r = h_wndg_r
        self.dims.update_dimensions(h_yoke_s, h_yoke_r, 
                                    h_wndg_s, h_wndg_r, 
                                    self.gn.delta_mag)
        if kwargs.get("coil_shapes", True):
            self.coil_shapes()
        
        
    def update_dimensions(self, **kwargs):
        self.h_yoke_s = kwargs.get("h_yoke_s", self.h_yoke_s)
        self.h_yoke_r = kwargs.get("h_yoke_r", self.h_yoke_r)
        self.h_wndg_s = kwargs.get("h_wndg_s", self.h_wndg_s)
        self.h_wndg_r = kwargs.get("h_wndg_r", self.h_wndg_r)
        self.dims.update_dimensions(self.h_yoke_s, self.h_yoke_r, 
                                    self.h_wndg_s, self.h_wndg_r, 
                                    self.gn.delta_mag)
        if kwargs.get("coil_shapes", True):
            self.coil_shapes()
        if kwargs.get("keep_const_kr_b", False):
            self.keep_const_kr_b(**kwargs)
            
    
    def update_model_by_K(self, K_s, K_r, **kwargs):
        
        self.K_s = K_s
        self.K_r = K_r
                
        alpha_r = kwargs.get("alpha_r", pi/2)
        alpha_s = kwargs.get("alpha_s", 0.0)
        
        self.ks_d = kwargs.get("ks_d", self.ks_d)
        self.ks_p = kwargs.get("ks_p", self.ks_p)
        self.kr_b = kwargs.get("kr_b", self.kr_b)
        
        this_mdl = Model(self.p, self.l_e)
        
        this_mdl.add_layer(AirLayer(r=self.dims.r_ri))
        this_mdl.add_layer(MagneticLayer(r=self.dims.r_ro, 
                                    mu_r=self.gn.mu_r_yoke))
   
        this_mdl.add_layer(CurrentLoading(K=K_r*np.sqrt(2)*self.kr_b, 
                                    r=self.dims.r_rw, 
                                    alpha=alpha_r,
                                    mu_r=1.
                                    ))

        this_mdl.add_layer(CurrentLoading(K=K_s*np.sqrt(2)*self.ks_d*self.ks_p,
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
        self.runs += 1
        
        this_plt = PlanePlot(this_mdl, fgsz=100)
        
        flux_data = this_mdl.get_B_data(r=np.array([self.dims.r_rw, 
                                                    self.dims.r_sw]), 
                                        t=np.linspace(0, np.pi/self.p, 400)
                                        )
        
        self.B_r = np.max(np.abs(flux_data.Br[:, 0]))
        self.B_s = np.max(np.abs(flux_data.Br[:, 1]))
        
        # self.B_rt = np.max(np.abs(flux_data.Bt[:, 0]))
        # self.B_st = np.max(np.abs(flux_data.Bt[:, 1]))
        
        arctan_r = np.arctan2(flux_data.Br[:,0], flux_data.Bt[:,0])
        B_r_norm = np.sqrt(flux_data.Br[:,0]**2 + flux_data.Bt[:,0]**2)
        phi_r = np.abs(pi/6 - arctan_r)
        self.B_r_c = np.max(B_r_norm * np.cos(phi_r))
        
        arctan_s = np.arctan2(flux_data.Br[:,1], flux_data.Bt[:,1])
        B_s_norm = np.sqrt(flux_data.Br[:,1]**2 + flux_data.Bt[:,1]**2)
        phi_s = np.abs(pi/6 - arctan_s)
        self.B_s_c = np.max(B_s_norm * np.cos(phi_s))
        
        # self.B_r_c = np.max(np.abs(flux_data.Br[:,0]))
        # self.B_s_c = np.max(np.abs(flux_data.Br[:,1]))
        
        # self.B_r_c = np.max([np.abs(  np.sin(np.pi/6)*flux_data.Br[:,0] + np.cos(np.pi/6)*flux_data.Bt[:, 0]),
        #                      np.abs(  np.sin(np.pi/6)*flux_data.Br[:,0] - np.cos(np.pi/6)*flux_data.Bt[:, 0]),
        #                      np.abs(- np.sin(np.pi/6)*flux_data.Br[:,0] + np.cos(np.pi/6)*flux_data.Bt[:, 0]),
        #                      np.abs(- np.sin(np.pi/6)*flux_data.Br[:,0] - np.cos(np.pi/6)*flux_data.Bt[:, 0]),
        #                     ])
        # self.B_s_c = np.max([np.abs(  np.sin(np.pi/6)*flux_data.Br[:,1] + np.cos(np.pi/6)*flux_data.Bt[:, 1]),
        #                      np.abs(  np.sin(np.pi/6)*flux_data.Br[:,1] - np.cos(np.pi/6)*flux_data.Bt[:, 1]),
        #                      np.abs(- np.sin(np.pi/6)*flux_data.Br[:,1] + np.cos(np.pi/6)*flux_data.Bt[:, 1]),
        #                      np.abs(- np.sin(np.pi/6)*flux_data.Br[:,1] - np.cos(np.pi/6)*flux_data.Bt[:, 1])
        #                     ])
        
        
        self.M = this_mdl.Mpos
        self.P = self.gn.w_syn * self.M
        
        
        self.mdl = this_mdl
        self.plt = this_plt
    
        
    
    def update_model_by_J(self, J_e_s, J_e_r, **kwargs):
        self.J_e_s = J_e_s
        self.J_e_r = J_e_r
        K_s = self.k_fill_s * self.h_wndg_s * self.J_e_s
        K_r = self.k_fill_r * self.h_wndg_r * self.J_e_r 
        self.update_model_by_K(K_s, K_r, **kwargs)
                
        
    def apply_lift_factor(self, **kwargs):
        
        lift_factor_angle = kwargs.get("lift_factor_angle", 30)
        verbose = kwargs.get("verbose", False)
        pretty = kwargs.get("pretty", False)
        
        J_s_history = [self.J_e_s]
        J_r_history = [self.J_e_r]   
        
        K_s_history = [self.K_s]
        K_r_history = [self.K_r]
        
        if pretty:
            B_s_c_history = [self.B_s_c]
            B_r_c_history = [self.B_r_c]
        
        self.K_s = self.k_fill_s * self.h_wndg_s * self.J_e_s
        self.K_r = self.k_fill_r * self.h_wndg_r * self.J_e_r
        
        J_s_history.append(self.J_e_s)
        J_r_history.append(self.J_e_r) 
        
        K_s_history.append(self.K_s)
        K_r_history.append(self.K_r)
        
        if pretty:
            B_s_c_history.append(self.B_s_c)
            B_r_c_history.append(self.B_r_c)
        
        self.update_model_by_K(self.K_s, self.K_r)
        
        limit = 50
        for idx in range(1, limit):
            B_s_c = self.B_s_c
            B_r_c = self.B_r_c            
            
            J_e_s=get_L_TPL2100(T = self.sw.T_HTS, B = B_s_c, 
                                theta = lift_factor_angle) * self.J_e_spec
            J_e_r=get_L_TPL2100(T = self.fw.T_HTS, B = B_r_c, 
                                theta = lift_factor_angle) * self.J_e_spec
            
            while np.isnan(J_e_s):
                print("\nINFO: lift_factor overflow for J_e_s...")
                B_s_c *= 0.9
                J_e_s=get_L_TPL2100(T = self.sw.T_HTS, B = B_s_c, 
                                    theta = lift_factor_angle) * self.J_e_spec
                
            if np.isnan(J_e_r):
                print("\nINFO: lift_factor overflow for J_e_r...")
                B_r_c *= 0.9
                J_e_r=get_L_TPL2100(T = self.fw.T_HTS, B = B_r_c, 
                                    theta = lift_factor_angle) * self.J_e_spec
            
            J_s_history.append(J_e_s)
            J_r_history.append(J_e_r)
                    
            self.update_model_by_J(J_e_s=J_s_history[-2] * 0.3 + J_s_history[-1] * 0.7, 
                                   J_e_r=J_r_history[-2] * 0.3 + J_r_history[-1] * 0.7, 
                                   **kwargs)
       
            K_s_history.append(self.K_s)
            K_r_history.append(self.K_r)
            
            if pretty:
                B_s_c_history.append(self.B_s_c)
                B_r_c_history.append(self.B_r_c)
            
            if idx > 20:
                if np.allclose(np.array(J_s_history[-4:-1]), self.J_e_s) and np.allclose(np.array(J_r_history[-4:-1]), self.J_e_r):
                        break
            if idx == limit -1:
                print("\nINFO: no convergence while applying lift factor ...\n")
        
            
        # --- plot convergence ---
        if verbose:
            fig, axs = pyplot.subplots(2)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(6)
            fig.set_figwidth(5)
            
            ax1 = axs[0]
            ax1.ticklabel_format(useOffset=False)
            ax1.grid()
            ax2 = ax1.twinx()
            ax2.ticklabel_format(useOffset=False)
            ax3 = axs[1]
            ax3.ticklabel_format(useOffset=False)
            ax3.grid()
            ax4 = ax3.twinx()
            ax4.ticklabel_format(useOffset=False)
            
            ax1.set_xticks(range(1,len(K_s_history)+1))
            ax1.scatter([i+1 for i in range(len(K_s_history))], 
                        [y*1e-3 for y in K_s_history], c="b")
            ax1.set_ylabel('K_s [kA/m]', color="b")
            ax2.scatter([i+1 for i in range(len(K_r_history))], 
                        [y*1e-3 for y in K_r_history], c="r", marker="x")
            ax2.set_ylabel('K_r [kA/m]', color="r")
            
            ax3.set_xticks(range(1,len(J_s_history)+1))
            ax3.scatter([i+1 for i in range(len(J_s_history))], 
                        [y*1e-6 for y in J_s_history], c="b")
            ax3.set_ylabel('J_e_s [A/mm2]', color="b")
            ax4.scatter([i+1 for i in range(len(J_r_history))], 
                        [y*1e-6 for y in J_r_history], c="r", marker="x")
            ax4.set_ylabel('J_e_r [A/mm2]', color="r")
            
        # --- plot convergence in single plots ---
        if pretty:
            n = 10
            w = 4
            h = 6
            do_tikz = False
            
            fig, ax = pyplot.subplots(1)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(w)
            fig.set_figwidth(h)
            # ax.ticklabel_format(useOffset=False)
            ax.grid()
            # ax.set_xticks(range(1,n+1))

            ax.scatter([i+1 for i in range(n)], 
                        [y*1e-3 for y in K_s_history[:n]], c="b", marker="o")
            ax.plot([i+1 for i in range(n)], 
                        [y*1e-3 for y in K_s_history[:n]], c="b", marker="o")
            ax.set_ylabel('Ks [kA/m]', color="b")
            if do_tikz:
                tkz.clean_figure()
                tkz.save("aus_iterating_Ks.tex")
            
            fig, ax = pyplot.subplots(1)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(w)
            fig.set_figwidth(h)
            # ax.ticklabel_format(useOffset=False)
            ax.grid()
            # ax.set_xticks(range(1,n+1))
            
            ax.scatter([i+1 for i in range(n)], 
                        [y*1e-6 for y in J_s_history[:n]], c="r")
            ax.plot([i+1 for i in range(n)], 
                        [y*1e-6 for y in J_s_history[:n]], c="r")
            ax.set_ylabel('J_e_s [A/mm2]', color="r")
            if do_tikz:
                tkz.clean_figure()
                tkz.save("aus_iterating_Jes.tex")
            
            fig, ax = pyplot.subplots(1)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(w)
            fig.set_figwidth(h)
            # ax.ticklabel_format(useOffset=False)
            ax.grid()
            # ax.set_xticks(range(1,n+1))

            ax.scatter([i+1 for i in range(n)], 
                        [y*1 for y in B_s_c_history[:n]], c="g", marker="o")
            ax.plot([i+1 for i in range(n)], 
                        [y*1 for y in B_s_c_history[:n]], c="g", marker="o")
            ax.set_ylabel('B_s_c [T]', color="g")
            if do_tikz:
                tkz.clean_figure()
                tkz.save("aus_iterating_Bsc.tex")
            
            # ----------- Rotor -----------
            fig, ax = pyplot.subplots(1)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(w)
            fig.set_figwidth(h)
            # ax.ticklabel_format(useOffset=False)
            ax.grid()
            # ax.set_xticks(range(1,n+1))
            
            ax.scatter([i+1 for i in range(n)], 
                        [y*1e-3 for y in K_r_history[:n]], c="b", marker="o")
            ax.plot([i+1 for i in range(n)], 
                        [y*1e-3 for y in K_r_history[:n]], c="b", marker="o")
            ax.set_ylabel('Kr [kA/m]', color="b")
            if do_tikz:
                tkz.clean_figure()
                tkz.save("aus_iterating_Kr.tex")
            
            fig, ax = pyplot.subplots(1)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(w)
            fig.set_figwidth(h)
            # ax.ticklabel_format(useOffset=False)
            ax.grid()
            # ax.set_xticks(range(1,n+1))
            
            ax.scatter([i+1 for i in range(n)], 
                        [y*1e-6 for y in J_r_history[:n]], c="r")
            ax.plot([i+1 for i in range(n)], 
                        [y*1e-6 for y in J_r_history[:n]], c="r")
            ax.set_ylabel('J_e_r [A/mm2]', color="r")
            if do_tikz:
                tkz.clean_figure()
                tkz.save("aus_iterating_Jer.tex")
            
            fig, ax = pyplot.subplots(1)            
            fig.tight_layout() 
            fig.dpi=500
            fig.set_figheight(w)
            fig.set_figwidth(h)
            # ax.ticklabel_format(useOffset=False)
            ax.grid()
            # ax.set_xticks(range(1,n+1))
            
            ax.scatter([i+1 for i in range(n)], 
                        [y*1 for y in B_r_c_history[:n]], c="g", marker="o")
            ax.plot([i+1 for i in range(n)], 
                        [y*1 for y in B_r_c_history[:n]], c="g", marker="o")
            ax.set_ylabel('B_r_c [T]', color="g")
            if do_tikz:
                tkz.clean_figure()
                tkz.save("aus_iterating_Brc.tex")

            

        
        
    def fast_flux(self, **kwargs):
        
        export_pdf = kwargs.get("pdf", False)
        export_svg = kwargs.get("svg", False)
        
        r_max = self.dims.r_so * 1.05
        show_axis=True
            
        
        self.plt.fluxplot(dr=1000, dt=300, lvls=10, lw=None,
                          t_min=0, t_max=np.pi/2,
                          custom_dims=True,
                          y0=0, x0=0,
                          y1=r_max, x1=r_max,
                          show_cbar=True, 
                          show_axis=show_axis,
                          show_borders=True,
                          border_width=0.2,
                          transparent=True,
                          padding=0,
                          pdf=export_pdf, pdf_dpi=300, 
                          svg=export_svg, svg_dpi=300, **kwargs)
        
    
    def coil_shapes(self, **kwargs):
        verbose = kwargs.get("verbose", False)

        # --- rotor ---
        w_rp = pi * self.dims.r_ro/self.p
        
        A_rc = self.k_fill_r * (self.h_wndg_r * w_rp)/2
        h_rc = self.h_wndg_r - 2 * self.gn.h_pole_frame
        if h_rc <= 0:
            if verbose: print(f"\nConfiguration invalid due to {h_rc = } [m]")

            h_rc = 0.01
        
        w_rc = A_rc/h_rc
        if 2 * (self.gn.w_pole_frame + w_rc) > w_rp:
            if verbose: print(f"\nConfiguration invalid due to {w_rc = } [m], {w_rp = } [m]")
            w_rc = (w_rp/2) - self.gn.w_pole_frame
        
        r_r_bend = 0.5 * (w_rp - 2 * self.gn.w_pole_frame - 2 * w_rc)
        if r_r_bend < self.gn.r_bend_max:
            if verbose: print(f"\nConfiguration invalid due to {r_r_bend = } [m]")
        
        # --- stator ---
        w_sp = (2 * pi * (self.dims.r_si - self.h_wndg_s)) / (3 * self.p)
        
        A_sc = self.k_fill_s * (self.h_wndg_s * w_sp)/2
        h_sc = self.h_wndg_s - 2 * self.gn.h_pole_frame
        if h_sc <= 0:
            if verbose: print(f"\nConfiguration invalid due to {h_sc = } [m]")

            h_sc = 0.01
            
        w_sc = A_sc/h_sc
        if 2 * (self.gn.w_pole_frame + w_sc) > w_sp:
            if verbose: print(f"\nConfiguration invalid due to {w_sc = } [m], {w_sp = } [m]")
            w_sc = (w_sp/2) - self.gn.w_pole_frame
        
        r_s_bend = 0.5 * (w_sp - 2 * self.gn.w_pole_frame - 2 * w_sc)
        if r_s_bend < self.gn.r_bend_max:
            if verbose: print(f"\nConfiguration invalid due to {r_s_bend = } [m]")
        

        self.coil = coil_Dimensions(w_rp, A_rc, h_rc, w_rc, r_r_bend, 
                                    w_sp, A_sc, h_sc, w_sc, r_s_bend)

            
    
    def keep_const_kr_b(self, **kwargs):
        kr_b = kwargs.get("kr_b", self.kr_b)
        self.k_fill_r = kr_b * (1 - (2*self.gn.h_pole_frame) / self.h_wndg_r)
        if self.k_fill_r > 1 or self.k_fill_r < 0:
            print(f"INFO: {self.k_fill_r = }")
    
    
    def apply_coil_sizes(self, **kwargs):
        self.coil_shapes(**kwargs)
        kr_b = 2 * self.coil.w_rc / self.coil.w_rp
        self.update_model_by_K(K_s=self.K_s,
                               K_r=self.K_r,
                               kr_b=kr_b)
        
    
    def apply_coil_sizes_and_lift_factor(self, **kwargs):
        
        self.apply_coil_sizes(**kwargs)
        self.apply_lift_factor(**kwargs)
        
        
    def HTS_usage(self, **kwargs):
        self.coil_shapes()
        l_rotor = self.l_e + pi * self.coil.r_r_bend
        l_stator = self.l_e + pi * self.coil.r_s_bend
        
        self.HTS_volume_rotor = l_rotor * self.p * 4 * self.coil.A_rc
        self.HTS_volume_stator = l_stator * self.p * 6 * self.coil.A_sc
        self.HTS_volume = self.HTS_volume_rotor + self.HTS_volume_stator
        self.HTS_length = self.HTS_volume / self.gn.A_tape
        self.HTS_length_rotor = self.HTS_volume_rotor / self.gn.A_tape
        self.HTS_length_stator = self.HTS_volume_stator / self.gn.A_tape
        return self.HTS_volume, self.HTS_length
    
    
    def compute_weight(self, **kwargs):
        
        weight_HTS_r = self.gn.HTS_weight_length_corr * self.HTS_length_rotor
        weight_HTS_s = self.gn.HTS_weight_length_corr * self.HTS_length_stator
        weight_Iron_r = self.gn.Iron_density * pi * (self.dims.r_ro**2 - self.dims.r_ri**2)
        weight_Iron_s = self.gn.Iron_density * pi * (self.dims.r_so**2 - self.dims.r_si**2)
        self.reference_weight = weight_HTS_r+ weight_Iron_r + weight_HTS_s + weight_Iron_s
        
                        
        self.weight_rotor_wndg = self.gn.HTS_density * self.HTS_volume_rotor
        self.weight_stator_wndg= self.gn.HTS_density * self.HTS_volume_stator
        
        self.weight_rotor_yoke = self.gn.Iron_density * pi * (self.dims.r_ro**2 - self.dims.r_ri**2)
        self.weight_stator_yoke= self.gn.Iron_density * pi * (self.dims.r_so**2 - self.dims.r_si**2)
        
        self.weight_rotor_kryo = self.gn.G10CR_density * self.l_e * (pi * ((self.dims.r_ro+self.h_wndg_r)**2 - self.dims.r_ro**2) - (self.p * 4 * self.coil.A_rc))
        self.weight_stator_kryo = self.gn.G10CR_density * self.l_e * (pi *(self.dims.r_si**2 - (self.dims.r_si-self.h_wndg_s)**2) - (self.p * 6 * self.coil.A_sc))
        
        self.weight_rotor = self.weight_rotor_yoke + self.weight_rotor_wndg + self.weight_rotor_kryo
        self.weight_stator= self.weight_stator_yoke+ self.weight_stator_wndg+ self.weight_stator_kryo
        self.weight = self.weight_rotor + self.weight_stator
                
    
    def guess_Ns(self, **kwargs):
        U_i = kwargs.get("U_i", self.gn.U_LL_N / np.sqrt(3))
        tau_p = (self.dims.r_si * pi) / self.p
        flux_data = self.mdl.get_B_data(r=np.array([self.dims.r_ag]), 
                                        t=np.linspace(0, np.pi/self.p, 400)
                                        )
        B_delta_hat = np.max(np.abs(flux_data.Br[:, 0]))
        
        num = U_i
        den = np.sqrt(2)*2 * self.f * self.ks_d * self.ks_p * self.l_e * tau_p * B_delta_hat
        
        N_s = num / den
        return N_s
    
    
    def guess_N_c_div_a(self, **kwargs):
        N_s = kwargs.get("N_s", self.guess_Ns(**kwargs))
        N_c_div_a = N_s / (2*self.p*self.q)
        return N_c_div_a
    
    
    def maximize_k_fill_s(self):
        r_si = self.dims.r_si
        h_sw = self.h_wndg_s
        p = self.p
        gn = self.gn
        nom = (pi*(r_si-h_sw)/(3*p)-gn.w_pole_frame-gn.r_bend_max)*(h_sw-2*gn.h_pole_frame)*3*p
        den = h_sw * (r_si - h_sw) * pi
        self.k_fill_s = nom / den
        
        
    def maximize_k_fill_r(self):
        r_ro = self.dims.r_ro
        h_rw = self.h_wndg_r
        p = self.p
        gn = self.gn
        nom = 2*p*(((pi * r_ro)/(2 * p)) - gn.w_pole_frame - gn.r_bend_max)*(h_rw - 2*gn.h_pole_frame)
        den = h_rw * pi * r_ro
        self.k_fill_r = nom / den
        
    
    def maximize_k_fill(self):
        self.maximize_k_fill_s()
        self.maximize_k_fill_r()
        
        
    def adapt_yokes(self, **kwargs):
        # --- adapt stator yoke ---
        h_yoke_s = yoke_height(self.p, self.dims.r_si, self.B_s, self.B_yoke_max)    
        # --- adapt rotor yoke ---
        h_yoke_r = yoke_height(self.p, self.dims.r_ro, self.B_r, self.B_yoke_max)
        # --- apply ---
        if h_yoke_s >= self.dims.r_so * 0.6:
            h_yoke_s = self.dims.r__so * 0.39
            print("WARNING: Stator Yoke too big.")
        self.update_dimensions(h_yoke_s=h_yoke_s, **kwargs)
        if h_yoke_r >= self.dims.r_ro:
            h_yoke_r = self.dims.r_ro * 0.9
            print("WARNING: Rotor Yoke too big.")
        self.update_dimensions(h_yoke_r=h_yoke_r, **kwargs)
            
        
    def apply_coil_thickness_ratio(self):
        self.maximize_k_fill()
        self.coil_shapes()
        # ctr = w_c / h_c = (0.2, 6)
        ctr_min = 0.2
        ctr_max = 6
        
        ctr_rotor = self.coil.w_rc / self.coil.h_rc
        if ctr_rotor < ctr_min: self.coil.w_rc = ctr_min * self.coil.h_rc
        if ctr_rotor > ctr_max: self.coil.w_rc = ctr_max * self.coil.h_rc
            
        # self.coil.w_rc = max(ctr_min * self.coil.h_rc, self.coil.w_rc)
        # self.coil.w_rc = min(ctr_max * self.coil.h_rc, self.coil.w_rc)
        
        ctr_stator = self.coil.w_sc / self.coil.h_sc
        if ctr_stator < ctr_min: self.coil.w_sc = ctr_min * self.coil.h_sc
        if ctr_stator > ctr_max: self.coil.w_sc = ctr_max * self.coil.h_sc
        
        # self.coil.w_sc = max(ctr_min * self.coil.h_sc, self.coil.w_sc)
        # self.coil.w_sc = min(ctr_max * self.coil.h_sc, self.coil.w_sc)
        
        self.coil.A_rc = self.coil.h_rc * self.coil.w_rc
        self.coil.A_sc = self.coil.h_sc * self.coil.w_sc
        
        self.coil.r_r_bend = 0.5 * (self.coil.w_rp - 2 * self.gn.w_pole_frame - 2 * self.coil.w_rc)
        self.coil.r_s_bend = 0.5 * (self.coil.w_sp - 2 * self.gn.w_pole_frame - 2 * self.coil.w_sc)
        
        self.kr_b = 2 * self.coil.w_rc / self.coil.w_rp
        self.k_fill_r = 2 * self.coil.A_rc / (self.h_wndg_r * self.coil.w_rp)
        self.k_fill_s = 2 * self.coil.A_sc / (self.h_wndg_s * self.coil.w_sp)
    
    
    def show_results(self, header="n7 - Results"):
        print("\n---", header, "---")
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
        

class coil_Dimensions():
    
    def __init__(self, w_rp, A_rc, h_rc, w_rc, r_r_bend,
                       w_sp, A_sc, h_sc, w_sc, r_s_bend):
        
        self.w_rp = w_rp 
        self.A_rc = A_rc
        self.h_rc = h_rc
        self.w_rc = w_rc
        self.r_r_bend = r_r_bend
        self.w_sp = w_sp
        self.A_sc = A_sc
        self.h_sc = h_sc
        self.w_sc = w_sc
        self.r_s_bend = r_s_bend
        
    def show(self, header="Coil - Dimensions"):
        print("\n---", header, "---")
        print(f"w_rp = {self.w_rp} [m]")
        print(f"A_rc = {self.A_rc} [m]")
        print(f"h_rc = {self.h_rc} [m]")
        print(f"w_rc = {self.w_rc} [m]")
        print(f"r_r_bend = {self.r_r_bend} [m]")
        print(f"w_sp = {self.w_sp} [m]")
        print(f"A_sc = {self.A_sc} [m]")
        print(f"h_sc = {self.h_sc} [m]")
        print(f"w_sc = {self.w_sc} [m]")
        print(f"r_s_bend = {self.r_s_bend} [m]")
            
    
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
        print("\n---", header, "---")
        print(f"r_so = {np.round(self.r_so, 4)} [m]")
        print(f"r_si = {np.round(self.r_si, 4)} [m]")
        print(f"r_sw = {np.round(self.r_sw, 4)} [m]")
        print(f"r_ag = {np.round(self.r_ag, 4)} [m]")
        print(f"r_rw = {np.round(self.r_rw, 4)} [m]")
        print(f"r_ro = {np.round(self.r_ro, 4)} [m]")
        print(f"r_ri = {np.round(self.r_ri, 4)} [m]")
        
    
    