# -*- coding: utf-8 -*-
import pathlib
from data import Generator as gn
from data import StatorWinding_77K
from design import Main_Params, Main_Dims, Main_Results

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
            print(f"h_wdng_r {res.h_wdng_r}[m] height winding rotor", file=text_file)
            print(f"h_wdng_s {res.h_wdng_s}[m] height winding stator", file=text_file)
            print(f"h_yoke_r {res.h_yoke_r}[m] height yoke rotor", file=text_file)
            print(f"h_yoke_s {res.h_yoke_s}[m] height yoke stator", file=text_file)
            
            print("PARAMS_DERIVED 1", file=text_file)
            
            print("r_so sqrt(d2so_L/l_e)/2 radius stator outer", file=text_file)
            print("r_si r_so-h_yoke_s radius stator inner", file=text_file)
            print("r_sA r_si-h_wdng_s/2 radius stator armature winding", file=text_file)
            print("r_rF r_sA-h_wdng_s/2-delta_mag-h_wdng_r/2 radius rotor field winding", file=text_file)
            print("r_ro r_rF-h_wdng_r/2 radius rotor outer", file=text_file)
            print("r_ri r_ro-h_yoke_r radius rotor inner", file=text_file)
            print("tau_p (pi*r_si)/p pole pitch", file=text_file)
            print("angle_tau_p 360[°]/(2*p) pole pitch angle", file=text_file)
            print("angle_sector 2*angle_tau_p sector angle", file=text_file)
            
            print("PARAMS_YET_TO_DEFINE 1", file=text_file)
            
            print("w_wdng_r 0.1[m] width winding rotor", file=text_file)
            print("w_wdng_r_pole 0.07[m] width winding rotor pole", file=text_file)
            print("h_wdng_r_margin 0.02[m] margin to yoke and airgap", file=text_file)
            print("w_wdng_s 0.07[m] width winding stator", file=text_file)
            print("w_wdng_s_pole 0.05[m] width winding stator pole", file=text_file)
            print("h_wdng_s_margin 0.02[m] margin to yoke and airgap", file=text_file)
            print("angle_wdng_r_shift (((r_rF*pi)/(2*p)-w_wdng_r_pole/2-w_wdng_r/2)/(2*r_rF*pi))*360[°] angle to shift the rotor windings", file=text_file)
            print("angle_wdng_s_shift ((w_wdng_s_pole/2+w_wdng_s/2)/(2*r_sA*pi))*360[°] angle to shift the stator winding", file=text_file)
            print("N_f 100 number of rotor winding turns", file=text_file)
            print("A_rotor_winding 1e-6[m^2]", file=text_file)
            
            
    else:
        with open(path, "w") as text_file:
        
            print("PARAMS_PREDEFINED 1", file=text_file)
            
            print(f"mu_r {gn.mu_r_yoke} relative permeability of iron", file=text_file)
            print("d2so_L 15[m^3] generator volume", file=text_file)
            print(f"delta_mag {gn.delta_mag}[m] height magnetic airgap", file=text_file)
            print(f"n_syn {gn.w_syn}[1/s] generator speed", file=text_file)
            print(f"n_syn {gn.n_syn}[rpm] generator speed", file=text_file)
            
            print("PARAMS_FROM_PYTHON 1", file=text_file)
            
            print("p 32 pole pair count", file=text_file)
            print("l_e 0.3[m] effective generator length", file=text_file)
            print("h_wdng_r 0.12[m] height winding rotor", file=text_file)
            print("h_wdng_s 0.068[m] height winding stator", file=text_file)
            print("h_yoke_r 0.1532[m] height yoke rotor", file=text_file)
            print("h_yoke_s 0.0843[m] height yoke stator", file=text_file)
            
            print("PARAMS_DERIVED 1", file=text_file)
            
            print("r_so sqrt(d2so_L/l_e)/2 radius stator outer", file=text_file)
            print("r_si r_so-h_yoke_s radius stator inner", file=text_file)
            print("r_sA r_si-h_wdng_s/2 radius stator armature winding", file=text_file)
            print("r_rF r_sA-h_wdng_s/2-delta_mag-h_wdng_r/2 radius rotor field winding", file=text_file)
            print("r_ro r_rF-h_wdng_r/2 radius rotor outer", file=text_file)
            print("r_ri r_ro-h_yoke_r radius rotor inner", file=text_file)
            print("tau_p (pi*r_si)/p pole pitch", file=text_file)
            print("angle_tau_p 360[°]/(2*p) pole pitch angle", file=text_file)
            print("angle_sector 2*angle_tau_p sector angle", file=text_file)
            
            print("PARAMS_YET_TO_DEFINE 1", file=text_file)
            
            print("w_wdng_r 0.1[m] width winding rotor", file=text_file)
            print("w_wdng_r_pole 0.07[m] width winding rotor pole", file=text_file)
            print("h_wdng_r_margin 0.02[m] margin to yoke and airgap", file=text_file)
            print("w_wdng_s 0.07[m] width winding stator", file=text_file)
            print("w_wdng_s_pole 0.05[m] width winding stator pole", file=text_file)
            print("h_wdng_s_margin 0.02[m] margin to yoke and airgap", file=text_file)
            print("angle_wdng_r_shift (((r_rF*pi)/(2*p)-w_wdng_r_pole/2-w_wdng_r/2)/(2*r_rF*pi))*360[°] angle to shift the rotor windings", file=text_file)
            print("angle_wdng_s_shift ((w_wdng_s_pole/2+w_wdng_s/2)/(2*r_sA*pi))*360[°] angle to shift the stator winding", file=text_file)
            print("N_f 100 number of rotor winding turns", file=text_file)
            print("A_rotor_winding 1e-6[m^2]", file=text_file)


if __name__ == "__main__":
    save_params_as_txt(filename="from_main", dims=None, params=None, res=None)