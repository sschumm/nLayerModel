# -*- coding: utf-8 -*-
import numpy as np

from data import Generator as gn
from data import FieldWinding as fw
from data import StatorWinding_77K as sw

from design import n7_Model

p = 32
l_e = 0.3
r_so = sw.r_so(sw, l_e)


h_yoke_s, h_yoke_r = 0.02 * r_so, 0.02 * r_so
h_wndg_s, h_wndg_r = gn.delta_mag, gn.delta_mag


generator = n7_Model(p, l_e, r_so, 
                     gn=gn, fw=fw, sw=sw,
                     k_fill_s=0.5,
                     k_fill_r=0.5,
                     B_yoke_max=6.0)
generator.update_dimensions(h_yoke_s, h_yoke_r, h_wndg_s, h_wndg_r)


K_s = generator.k_fill_s * generator.J_e_s * generator.h_wndg_s
K_r = generator.k_fill_r * generator.J_e_r * generator.h_wndg_r


generator.update_model(K_s = K_s, K_r = K_r, alpha_r=0.5)
# generator.plt.fluxplot(1000,1000, lvls=10)
generator.show_results()