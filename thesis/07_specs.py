# -*- coding: utf-8 -*-
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

