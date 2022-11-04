# -*- coding: utf-8 -*-
import pathlib
import datetime
from data import Generator as gn
from data import StatorWinding_77K
from design import Main_Params, Main_Dims, Main_Results
from design import n7_Model

def save_params(filename: str, model: n7_Model, **kwargs):
    
    now = datetime.datetime.now()
    
    if len(str(now.hour)) == 1:
        hour = "0"+str(now.hour)
    else:
        hour = now.hour

    if len(str(now.minute)) == 1:
        minute = "0"+str(now.minute)
    else:
        minute = now.minute
        
    if len(str(now.day)) == 1:
        day = "0"+str(now.day)
    else:
        day = now.day
    
    timestamp = f"{now.month}{day}_{hour}{minute}_"
    
    
    path = pathlib.Path(__file__).parent / "comsol" / (timestamp + filename + ".txt")
    sw = kwargs.get("sw", StatorWinding_77K)
    
    with open(path, "w") as text_file:
        print("PARAMS_PREDEFINED 1", file=text_file)
        
        print(f"mu_r {model.gn.mu_r_yoke} relative permeability of iron", file=text_file)
        print(f"d2so_L {model.sw.d2so_L}[m^3] generator volume", file=text_file)
        print(f"delta_mag {model.gn.delta_mag}[m] height magnetic airgap", file=text_file)
        print(f"n_syn {model.gn.n_syn}[rpm] generator speed", file=text_file)
        print(f"h_pole_frame {model.gn.h_pole_frame}[m] height pole frame", file=text_file)
        print(f"w_pole_frame {model.gn.w_pole_frame}[m] width pole frame", file=text_file)       
        
        
        print("PARAMS_FROM_PYTHON 1", file=text_file)
        
        print(f"p {model.p} pole pair count", file=text_file)
        print(f"l_e {model.l_e}[m] effective generator length", file=text_file)
        print(f"h_wdng_r {model.h_wndg_r}[m] height winding rotor", file=text_file)
        print(f"h_wdng_s {model.h_wndg_s}[m] height winding stator", file=text_file)
        print(f"h_yoke_r {model.h_yoke_r}[m] height yoke rotor", file=text_file)
        print(f"h_yoke_s {model.h_yoke_s}[m] height yoke stator", file=text_file)
        print(f"k_fill_r {model.k_fill_r} fill factor rotor", file=text_file)
        print(f"k_fill_s {model.k_fill_s} fill factor stator", file=text_file)
        print(f"J_e_r {model.J_e_r * 1e-6}[A/mm^2] Engineering Current Density Rotor", file=text_file)
        print(f"J_e_s {model.J_e_s * 1e-6}[A/mm^2] Engineering CUrrent Density Stator", file=text_file)
        
        print(f"r_so {model.dims.r_so}[m] radius stator outer", file=text_file)
        print(f"r_si {model.dims.r_si}[m] radius stator inner", file=text_file)
        print(f"r_sw {model.dims.r_sw}[m] radius stator armature winding", file=text_file)
        print(f"r_ag {model.dims.r_ag}[m] radius airgap", file=text_file)
        print(f"r_rw {model.dims.r_rw}[m] radius rotor field winding", file=text_file)
        print(f"r_ro {model.dims.r_ro}[m] radius rotor outer", file=text_file)
        print(f"r_ri {model.dims.r_ri}[m] radius rotor inner", file=text_file)
        

        print(f"w_rp {model.coil.w_rp}[m] width of the rotor pole", file=text_file)
        print(f"A_rc {model.coil.A_rc}[m^2] Area of one rotor coil per pole", file=text_file)
        print(f"h_rc {model.coil.h_rc}[m] height of one rotor coil per pole", file=text_file)
        print(f"w_rc {model.coil.w_rc}[m] width of one rotor coil per pole", file=text_file)
        print(f"w_sp {model.coil.w_sp}[m] width of the stator pole", file=text_file)
        print(f"A_sc {model.coil.A_sc}[m^2] Area of one stator coil per pole pair", file=text_file)
        print(f"h_sc {model.coil.h_sc}[m] height of one stator coil per pole pair", file=text_file)
        print(f"w_sc {model.coil.w_sc}[m] width of one stator coil per pole pair", file=text_file)
        
        print(f"N_s {model.N_s} number of turns per phase", file=text_file)
        print(f"J_s_amplitude {model.J_s_amplitude}[A/m^2] amplitude of the stator current density", file=text_file)
        
        
        
        print("PARAMS_DERIVED 1", file=text_file)
        
        print(f"omega 2*pi*n_syn*p electrical angular frequency", file=text_file)
        print(f"T_el 1/(n_syn*p) electrical period", file=text_file)
        print(f"T_el_noUnit T_el*1[1/s] electrical period without unit", file=text_file)
        print("tau_p (pi*r_si)/p pole pitch", file=text_file)
        print("phi_tau_p 360[°]/(2*p) pole pitch angle", file=text_file)
        print("phi_sector 2*phi_tau_p sector angle", file=text_file)
        print("phi_rw1 phi_tau_p/2 angle for lower rotor winding", file=text_file)
        print("phi_rw2 phi_tau_p*(3/2) angle for upper rotor winding", file=text_file)
        print("psi_sw_W phi_sector*(1/3) angle for stator winding W", file=text_file)
        print("psi_sw_V phi_sector*(2/3) angle for stator winding V", file=text_file)
        print("alpha 0 torque generating angle", file=text_file)
        print("TDC_steps 50 steps in time dependant coarse study", file=text_file)
        print("TDF_steps 300 steps in time dependant fine study", file=text_file)
        
    print("\nINFO: data export succesfull")


def save_params_as_txt(filename: str, dims: Main_Dims, 
                       params: Main_Params, res: Main_Results, **kwargs):
    
    path = pathlib.Path(__file__).parent / "comsol" / (filename + ".txt")
    sw = kwargs.get("sw", StatorWinding_77K)
    
    if dims is not None and params is not None and res is not None:
        
        with open(path, "w") as text_file:
        
            print("PARAMS_PREDEFINED 1", file=text_file)
            
            print(f"mu_r {gn.mu_r_yoke} relative permeability of iron", file=text_file)
            print(f"d2so_L {sw.d2so_L}[m^3] generator volume", file=text_file)
            print(f"delta_mag {gn.delta_mag}[m] height magnetic airgap", file=text_file)
            print(f"n_syn {gn.n_syn}[rpm] generator speed", file=text_file)
            print(f"n_syn {gn.w_syn}[1/s] generator speed", file=text_file)
            
            print("PARAMS_FROM_PYTHON 1", file=text_file)
            
            print(f"p {params.p} pole pair count", file=text_file)
            print(f"l_e {params.l_e}[m] effective generator length", file=text_file)
            print(f"h_wndg_r {res.h_wndg_r}[m] height winding rotor", file=text_file)
            print(f"h_wndg_s {res.h_wndg_s}[m] height winding stator", file=text_file)
            print(f"h_yoke_r {res.h_yoke_r}[m] height yoke rotor", file=text_file)
            print(f"h_yoke_s {res.h_yoke_s}[m] height yoke stator", file=text_file)
            
            print("PARAMS_DERIVED 1", file=text_file)
            
            print("r_so sqrt(d2so_L/l_e)/2 radius stator outer", file=text_file)
            print("r_si r_so-h_yoke_s radius stator inner", file=text_file)
            print("r_sA r_si-h_wndg_s/2 radius stator armature winding", file=text_file)
            print("r_rF r_sA-h_wndg_s/2-delta_mag-h_wndg_r/2 radius rotor field winding", file=text_file)
            print("r_ro r_rF-h_wndg_r/2 radius rotor outer", file=text_file)
            print("r_ri r_ro-h_yoke_r radius rotor inner", file=text_file)
            print("tau_p (pi*r_si)/p pole pitch", file=text_file)
            print("angle_tau_p 360[°]/(2*p) pole pitch angle", file=text_file)
            print("angle_sector 2*angle_tau_p sector angle", file=text_file)
            
            print("PARAMS_YET_TO_DEFINE 1", file=text_file)
            
            print("w_wdng_r 0.1[m] width winding rotor", file=text_file)
            print("w_wdng_r_pole 0.07[m] width winding rotor pole", file=text_file)
            print("h_wndg_r_margin 0.02[m] margin to yoke and airgap", file=text_file)
            print("w_wdng_s 0.07[m] width winding stator", file=text_file)
            print("w_wdng_s_pole 0.05[m] width winding stator pole", file=text_file)
            print("h_wndg_s_margin 0.02[m] margin to yoke and airgap", file=text_file)
            print("angle_wdng_r_shift (((r_rF*pi)/(2*p)-w_wdng_r_pole/2-w_wdng_r/2)/(2*r_rF*pi))*360[°] angle to shift the rotor windings", file=text_file)
            print("angle_wdng_s_shift ((w_wdng_s_pole/2+w_wdng_s/2)/(2*r_sA*pi))*360[°] angle to shift the stator winding", file=text_file)
            print("N_f 100 number of rotor winding turns", file=text_file)
            print("A_rotor_winding 1e-6[m^2]", file=text_file)
            
            
    else:
        with open(path, "w") as text_file:
        
            print("from main")

if __name__ == "__main__":
    save_params_as_txt(filename="from_main", dims=None, params=None, res=None)