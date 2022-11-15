# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib as tkz
from design import get_L_TPL2100 as L

P = 7000
n = 8.33

I_c = 300
A_HTS = 2.64
J_base = I_c / (1.7 * A_HTS)

J_cu = 3.8
d2l_cu = 22.3
C_cu = 37.71

#%% 77 K

T_77 = 77
B_77 = 1.3

J_77 = J_base * L(T_77, B_77, 30)
print(f"{J_77 = }")
C_77 = C_cu * J_77/J_cu
print(f"{C_77 = }")
C_77_a = C_77/2
print(f"{C_77_a = }")
d2l_77 = P / (C_77_a * n)
print(f"{d2l_77 = }")

#%% 30 K

T_30 = 30
B_30 = 1.8

J_30 = J_base * L(T_30, B_30, 30)
print(f"{J_30 = }")
C_30 = C_cu * J_30/J_cu
print(f"{C_30 = }")
C_30_a = C_30/2
print(f"{C_30_a = }")
d2l_30 = P / (C_30_a * n)
print(f"{d2l_30 = }")

#%% visualize

theta = 30 # [°]

lst_T = [77, 30] # [T]
lst_B = np.linspace(0,4,200) # [B]

lst_J_e = []


for T in lst_T:
    
    a = []
    for B in lst_B:  
        
        J_e = L(T, B, theta) * J_base
        a.append(J_e)

    lst_J_e.append(a)
    
#%% plot 77K

fig = plt.figure(dpi=300, figsize=(6,5))
ax = plt.subplot()
ax.set_xlabel("B / T")
ax.set_ylabel("J_e / Amm2")
ax.set_xlim(0, 3)
ax.set_ylim(0, 25)
# ax.set_xticks([i for i in range(0, 181, 20)])
# ax.set_yticks([i for i in range(0, 2501, 250)])
plt.grid()

ax.plot(lst_B, lst_J_e[0], label=f"T = {lst_T[0]} / K")

# ax.legend(loc="upper right")

tkz.clean_figure()
# tkz.save("aus_J_e_over_B_77K.tex")

#%% plot 30K

fig = plt.figure(dpi=300, figsize=(6,5))
ax = plt.subplot()
ax.set_xlabel("B / T")
ax.set_ylabel("J_e / Amm2")
ax.set_xlim(0, 3)
ax.set_ylim(0, 400)
# ax.set_xticks([i for i in range(0, 181, 20)])
ax.set_yticks([i for i in range(0, 401, 80)])
plt.grid()

ax.plot(lst_B, lst_J_e[1], label=f"T = {lst_T[1]} / K")

# ax.legend(loc="upper right")

tkz.clean_figure()
# tkz.save("aus_J_e_over_B_30K.tex")