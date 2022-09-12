# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True, linewidth=250, precision=2)

from scipy.constants import pi
from modules import Model, MagneticLayer, AirLayer, CurrentLoading
from modules.plot import RadialMultiPlot, PlanePlot
from analytics import B_from, Phi_from, taup, K_from
from data import Generator as gn
from data import FieldWinding as fw
# --------- using the standard copper specs -----------------
from data import StatorWinding_Cu as sw


def new_model(**kwargs):
    
    # ---------------- d2so_L ----------------------------
    if "r_so" in kwargs:
        r_so= kwargs.get("r_so", 2)
        l_e = sw.d2so_L/(2 * r_so)**2
    else:
        l_e = kwargs.get("l_e", 0.8)
        r_so = 0.5 * np.sqrt(sw.d2so_L/l_e)
        
        
    # ---------- first assumptions (based on song gen#1) --------------
    p = kwargs.get("p", 12)
    k_fill_r = kwargs.get("kf_r", 0.5)
    k_fill_s = kwargs.get("kf_s", 0.5)
    h_yoke_s = kwargs.get("hy_s", 0.08)
    h_winding_s = kwargs.get("hw_s", 0.055)
    B_yoke_max = kwargs.get("By_max", 4)

    
    # ------- first derivations --------------------
    r_si = r_so - h_yoke_s - h_winding_s # fill winding space with iron due to cu-winding
    r_a = r_si # - h_winding_s/2
    pole_pitch = taup(d_si=2*r_si, p=p)
    
    K_s = k_fill_s * h_winding_s * sw.J_e
    K_s_amplitude = np.sqrt(2) * K_s
    
    B_airgap_amplitude = B_from(K=K_s, C_e=sw.C_e)
    estimated_flux_per_pole = Phi_from(B_airgap_amplitude, 
                                       taup=pole_pitch, l_e=l_e)
    estimated_flux_in_yoke = estimated_flux_per_pole/2
    h_yoke_r1 = estimated_flux_in_yoke / (l_e * B_yoke_max)
    
    h_winding_r1 = h_winding_s
    r_ro1 = r_si - h_winding_s - h_winding_r1 - gn.delta_mag
    r_ri1 = r_ro1 - h_yoke_r1
    r_f1 = r_ro1 + h_winding_r1/2
    airgap1 = r_ro1 + (r_si - r_ro1)/2
    
    
    # ------ build model with only stator current loading -------
    mdl1 = Model(p=p, l=l_e)
    mdl1.add_layer(AirLayer(r=r_ri1))
    mdl1.add_layer(MagneticLayer(r=r_ro1, 
                                mu_r=gn.mu_r_yoke))
    # mdl1.add_layer(CurrentLoading(K=K_r_amplitude, 
    #                             r=r_f, 
    #                             alpha=0.0, 
    #                             mu_r=1
    #                             ))
    mdl1.add_layer(CurrentLoading(K=K_s_amplitude,
                                r=r_a,
                                alpha=0.0,
                                mu_r=1.
                                ))
    # mdl1.add_layer(AirLayer(r=r_si))
    mdl1.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl1.build()
    mdl1.solve()
    mdl1.total_torque()
    
    # ---------- get flux data from stator current loading model ----------
    Bdata = mdl1.get_B_data(r=np.array([airgap1]), 
                           t=np.linspace(0,2*pi,10000))
    B_airgap_amplitude_measured = np.sqrt(2) * np.max(np.sqrt(Bdata.Br**2 + Bdata.Bt**2))
    
    K_r = K_from(B_airgap_amplitude_measured, C_e=sw.C_e)
    K_r_amplitude = np.sqrt(2) * K_r
    
    # --------- redo the rotor dimensions with new data
    h_winding_r2 = K_r / (k_fill_r * fw.J_e)
    
    r_ro2 = r_si - h_winding_s - h_winding_r2 - gn.delta_mag
    r_ri2 = r_ro2 - h_yoke_r1
    r_f2 = r_ro2 + h_winding_r2/2
    # airgap = r_ro2 + (r_si - r_ro2)/2
    
    # ------------ build model with only rotor current loading ----------
    mdl2 = Model(p=p, l=l_e)
    mdl2.add_layer(AirLayer(r=r_ri2))
    mdl2.add_layer(MagneticLayer(r=r_ro2, 
                                mu_r=gn.mu_r_yoke))
    mdl2.add_layer(CurrentLoading(K=K_r_amplitude, 
                                r=r_f2, 
                                alpha=0.0, 
                                mu_r=1.
                                ))
    # mdl2.add_layer(CurrentLoading(K=K_s_amplitude,
    #                             r=r_a,
    #                             alpha=0.0,
    #                             mu_r=1.
    #                             ))
    mdl2.add_layer(AirLayer(r=r_si))
    mdl2.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl2.build()
    mdl2.solve()
    mdl2.total_torque()
    
    # --------- build reference model with only rotor current loading ------
    mdl3 = Model(p=p, l=l_e)
    mdl3.add_layer(AirLayer(r=r_ri1))
    mdl3.add_layer(MagneticLayer(r=r_ro1, 
                                mu_r=gn.mu_r_yoke))
    mdl3.add_layer(CurrentLoading(K=K_r_amplitude, 
                                r=r_f1, 
                                alpha=0.0, 
                                mu_r=1.
                                ))
    # mdl3.add_layer(CurrentLoading(K=K_s_amplitude,
    #                             r=r_a,
    #                             alpha=0.0,
    #                             mu_r=1.
    #                             ))
    mdl3.add_layer(AirLayer(r=r_si))
    mdl3.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl3.build()
    mdl3.solve()
    mdl3.total_torque()
    
    # ------------- compute full model -------------
    mdl4 = Model(p=p, l=l_e)
    mdl4.add_layer(AirLayer(r=r_ri2))
    mdl4.add_layer(MagneticLayer(r=r_ro2, 
                                mu_r=gn.mu_r_yoke))
    mdl4.add_layer(CurrentLoading(K=K_r_amplitude, 
                                r=r_f2, 
                                alpha=0.5, 
                                mu_r=1.
                                ))
    mdl4.add_layer(CurrentLoading(K=K_s_amplitude,
                                r=r_a,
                                alpha=0.0,
                                mu_r=1.
                                ))
    mdl4.add_layer(AirLayer(r=r_si))
    mdl4.add_layer(MagneticLayer(r=r_so, 
                                mu_r=gn.mu_r_yoke))
    mdl4.build()
    mdl4.solve()
    mdl4.total_torque()
    
    return mdl1, mdl2, mdl3, mdl4

dr, dt = 1000, 1000
models = new_model()

#%%
# p_plt = PlanePlot(models[0])
# p_plt.fluxplot(dr, dt, lvls=8)
# p_plt.quiver(dr=20, dt=200, scale=20, width=0.001)

# #%%
# p_plt = PlanePlot(models[1])
# p_plt.fluxplot(dr, dt, lvls=8)
# p_plt.quiver(dr=20, dt=200, scale=100, width=0.001)

# #%%
# p_plt = PlanePlot(models[2])
# p_plt.fluxplot(dr, dt, lvls=8)
# p_plt.quiver(dr=20, dt=200, scale=100, width=0.001)

#%%
p_plt = PlanePlot(models[3])
# p_plt.fluxplot(dr, dt, lvls=8)
# p_plt.quiver(dr=20, dt=200, scale=100, width=0.001)
M = models[3].Mpos
print(f"{M = } [Nm]")
P_out = gn.w_syn * M/1e6
print(f"{P_out = } [MW]")

