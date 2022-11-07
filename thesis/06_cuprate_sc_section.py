# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib as tkz

from scipy.constants import pi
from design import get_L_TPL2100

I_c_spec = 300 # [A]

#%% plot 1 - I_c(theta) at constant B


B = 2 # [T]

lst_T = [10, 30, 77] # [K]
lst_theta = [i for i in range(181)] # [°]

lst_I_c = []


for T in lst_T:
    
    a = []
    for theta in lst_theta:    
        
        I_c = get_L_TPL2100(T, B, theta) * I_c_spec
        a.append(I_c)

    lst_I_c.append(a)
    

fig = plt.figure(dpi=300, figsize=(5,7))
ax = plt.subplot()
ax.set_xlim(0, 180)
ax.set_ylim(0, 2500)
ax.set_xticks([i for i in range(0, 181, 20)])
ax.set_yticks([i for i in range(0, 2501, 250)])
plt.grid()

for i in range(len(lst_I_c)):
    ax.plot(lst_theta, lst_I_c[i], label=f"T = {lst_T[i]} / K")

ax.legend(loc="upper left")


#%% plot 2 - I_c(T) at constant theta with some flux densities B


theta = 30 # [T]

lst_B = [1, 2, 3] # [°]
lst_T = [i for i in range(10, 81)] # [T]

lst_I_c = []


for B in lst_B:
    
    a = []
    for T in lst_T:  
        
        I_c = get_L_TPL2100(T, B, theta) * I_c_spec
        a.append(I_c)

    lst_I_c.append(a)
    

fig = plt.figure(dpi=300, figsize=(6,5))
ax = plt.subplot()
ax.set_xlabel("T / K")
ax.set_ylabel("I_c / A")
ax.set_xlim(10, 80)
ax.set_ylim(0, 1650)
# ax.set_xticks([i for i in range(0, 181, 20)])
# ax.set_yticks([i for i in range(0, 2501, 250)])
plt.grid()

for i in range(len(lst_I_c)):
    ax.plot(lst_T, lst_I_c[i], label=f"B = {lst_B[i]} / T")

ax.legend(loc="upper right")


tkz.clean_figure()
tkz.save("sc_I_c_ofT_at_const_theta_and_B.tex")