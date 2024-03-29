import numpy as np
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import PlanePlot
from analytics import taup
from data import Generator as gn



def create_n_Layer_model(dims, p, l, Ks=False, Kr=False, **kwargs):
    
    is_stator_aircore = kwargs.get("is_stator_aircore", True)
    ks_d = kwargs.get("ks_d", 1)
    ks_p = kwargs.get("ks_p", 1)
    kr_b = kwargs.get("kr_b", 1)
    J_e_r = kwargs.get("J_e_r", None)
    J_e_s = kwargs.get("J_e_s", None)

    r_so = dims.r_so
    r_si = dims.r_si
    r_sA = dims.r_sA
    r_rF = dims.r_rF
    r_ro = dims.r_ro
    r_ri = dims.r_ri

    mdl = Model(p=p, l=l)
    mdl.add_layer(AirLayer(r=r_ri))
    mdl.add_layer(MagneticLayer(r=r_ro, 
                                mu_r=gn.mu_r_yoke))
    if Kr:    
        mdl.add_layer(CurrentLoading(K=Kr*np.sqrt(2)*kr_b, 
                                    r=r_rF, 
                                    alpha=np.pi/2, 
                                    mu_r=1.
                                    ))
    else:
        Kr = 0
    if Ks:
        mdl.add_layer(CurrentLoading(K=Ks*np.sqrt(2)*ks_d*ks_p,
                                    r=r_sA,
                                    alpha=0.0,
                                    mu_r=1.
                                    ))
    else:
        Ks = 0
    if is_stator_aircore:
        mdl.add_layer(AirLayer(r=r_si))
    else:
        mdl.add_layer(AirLayer(r=r_sA))
    mdl.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl.build()
    mdl.solve()
    mdl.total_torque()
    p_plt = PlanePlot(mdl, fgsz=100)
    
    airgap = r_ro + 2*(r_rF-r_ro) + 0.5 * ((r_si - 2*(r_si - r_sA)) - (r_ro + 2*(r_rF-r_ro)))
    
    flux_data = mdl.get_B_data(r=np.array([r_rF, airgap, r_sA]), 
                                t=np.linspace(0,np.pi/p, 400))
    B_r = np.max(np.abs(flux_data.Br[:, 0]))
    B_airgap = np.max(np.abs(flux_data.Br[:, 1]))
    B_s = np.max(np.abs(flux_data.Br[:, 2]))
    
    B_r_c = np.max(np.abs(np.sin(np.pi/6) * flux_data.Br[:,0] + np.cos(np.pi/6) * flux_data.Bt[:, 0]))
    B_s_c = np.max(np.abs(np.sin(np.pi/6) * flux_data.Br[:,2] + np.cos(np.pi/6) * flux_data.Bt[:, 2]))
    
    res = Main_Results(mdl.Mpos, mdl.Mpos*gn.w_syn, Kr, Ks, 
                       h_wndg_r=2*(r_rF-r_ro), 
                       h_wndg_s=2*(r_si-r_sA),
                       h_yoke_r=r_ro-r_ri, 
                       h_yoke_s=r_so-r_si,
                       B_r=B_r,
                       B_airgap=B_airgap,
                       B_s=B_s,
                       B_r_c=B_r_c,
                       B_s_c=B_s_c,
                       J_e_r=J_e_r,
                       J_e_s=J_e_s)
    return mdl, p_plt, res


class Main_Params():
    
    def __init__(self, p, l_e, r_so, r_si, k_fill_r=0.5, k_fill_s=0.5, B_yoke_max=2.5):
        self.p = p
        self.l_e = l_e
        self.r_so = r_so
        self.r_si = r_si
        
        self.k_fill_r = k_fill_r
        self.k_fill_s = k_fill_s
        self.B_yoke_max = B_yoke_max
        self.pole_pitch = taup(d_si = 2 * self.r_si, p = self.p)
        
    def show(self, header="Main_Params"):
        print("---", header, "---")
        print(f"p   = {self.p} [-]")
        print(f"l_e = {np.round(self.l_e, 3)} [m]")
        print(f"r_so = {np.round(self.r_so, 3)} [m]")
        print(f"r_si = {np.round(self.r_si, 3)} [m]")
        print(f"k_fill_r = {np.round(self.k_fill_r, 3)} [-]")
        print(f"k_fill_s = {np.round(self.k_fill_s, 3)} [-]")
        print(f"B_yoke_max = {np.round(self.B_yoke_max, 3)} [T]")
        print(f"pole_pitch = {np.round(self.pole_pitch, 3)} [m]")
        
        

class Main_Dims():
    
    def __init__(self, r_so, r_si, r_sA, r_rF, r_ro, r_ri):
        self.r_so = r_so
        self.r_si = r_si
        self.r_sA = r_sA
        self.r_rF = r_rF
        self.r_ro = r_ro
        self.r_ri = r_ri
     
    def show(self, header="Main_Dims"):
        print("---", header, "---")
        print(f"r_so = {np.round(self.r_so, 4)} [m]")
        print(f"r_si = {np.round(self.r_si, 4)} [m]")
        print(f"r_sA = {np.round(self.r_sA, 4)} [m]")
        print(f"r_rF = {np.round(self.r_rF, 4)} [m]")
        print(f"r_ro = {np.round(self.r_ro, 4)} [m]")
        print(f"r_ri = {np.round(self.r_ri, 4)} [m]")
        
        
class Main_Results():
    
    def __init__(self, M, P, K_r, K_s, h_wndg_r, h_wndg_s, h_yoke_r, h_yoke_s, 
                 B_r=None, B_airgap=None, B_s=None, B_r_c=None, B_s_c=None,
                 J_e_r=None, J_e_s=None):
        self.M = np.round(M, 2)
        self.P = np.round(P, 2)
        self.K_r = np.round(K_r, 2)
        self.K_s = np.round(K_s, 2)
        self.h_wndg_r = np.round(h_wndg_r, 4)
        self.h_wndg_s = np.round(h_wndg_s, 4)                 
        self.h_yoke_r = np.round(h_yoke_r, 4)
        self.h_yoke_s = np.round(h_yoke_s, 4)
        self.B_r = B_r
        self.B_airgap = B_airgap
        self.B_s = B_s
        self.B_r_c = B_r_c
        self.B_s_c = B_s_c
        self.J_e_r = J_e_r
        self.J_e_s = J_e_s
    
         
    def show(self, header="Main_Results"):
        print("---", header, "---")
        # print(f"M = {np.round(self.M * 1e-6, 2)} [MNm]")
        print(f"P = {np.round(self.P * 1e-6, 2)} [MW]")
        print(f"K_r = {np.round(self.K_r * 1e-3, 2)} [kA/m]")
        print(f"K_s = {np.round(self.K_s * 1e-3, 2)} [kA/m]")
        print(f"h_wndg_r = {self.h_wndg_r} [m]")
        print(f"h_wndg_s = {self.h_wndg_s} [m]")
        print(f"h_yoke_r = {self.h_yoke_r} [m]")
        print(f"h_yoke_s = {self.h_yoke_s} [m]")
        # print(f"B_r = {self.B_r} [T]")
        # print(f"B_airgap = {self.B_airgap} [T]")
        # print(f"B_s = {self.B_s} [T]")
        print(f"B_r_c = {np.round(self.B_r_c, 4)} [T]")
        print(f"B_s_c = {np.round(self.B_s_c, 4)} [T]")
        if self.J_e_r is not None:
            print(f"J_e_r = {np.round(self.J_e_r * 1e-6, 4)} [A/mm2]")
        if self.J_e_s is not None:
            print(f"J_e_s = {np.round(self.J_e_s * 1e-6, 4)} [A/mm2]")
        
        
    def show_P(self, header="P"):
        print(f"--- {header} --> P = {self.P/1e6} [MW]")
        
    def show_B_r(self, header="B_r"):
        print(f"--- {header} --> B_r = {self.B_r} [T]")
        
    def show_B_s(self, header="B_s"):    
        print(f"--- {header} --> B_s = {self.B_s} [T]")
        
    def show_B_airgap(self, header="B_airgap"):        
        print(f"--- {header} --> B_airgap = {self.B_airgap} [T]")
        
    def get_some(self):
        data = ""
        data += f"P= {np.round(self.P * 1e-6, 2)} [MW]" + "\n"
        data += f"B= {self.B_airgap} [T]" + "\n"
        data += f"K_r= {np.round(self.K_r * 1e-3, 2)} [kA/m]" + "\n"
        data += f"K_s= {np.round(self.K_s * 1e-3, 2)} [kA/m]" + "\n"
        return data
    