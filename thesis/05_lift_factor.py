# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from data import Generator as gn
from design import get_L_TPL2100

#%% vary only I_c_spec -- 77K

T = 77 # [K]
B = 4.5 # [T]
theta = 30 # [°]
I_c_spec = 650 # [A]

x = np.linspace(300, 650,10)
y = []
    
for I_c_spec in x:
    J_e_c = get_L_TPL2100(T, B, theta) * I_c_spec/(1.7 * gn.A_tape * 1e6)
    y.append(J_e_c)

fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")
ax.set_xlabel("I_c_spec")
ax.set_ylabel("J_e_c")
ax.plot(x, y)

#%% vary only B -- 77K

T = 77 # [K]
B = 4.5 # [T]
theta = 30 # [°]
I_c_spec = 650 # [A]

x = np.linspace(2, 4.5, 100)
y = []
    
for B in x:
    J_e_c = get_L_TPL2100(T, B, theta) * I_c_spec/(1.7 * gn.A_tape * 1e6)
    y.append(J_e_c)

fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")
ax.set_xlabel("B")
ax.set_ylabel("J_e_c")
ax.plot(x, y)


#%% vary both B and I_c_spec -- 77K

T = 77 # [K]
B = 4.5 # [T]
theta = 30 # [°]
I_c_spec = 650 # [A]

x = np.linspace(2, 4.5, 100)
c = np.linspace(300, 650, 5)
y = []

fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")
ax.set_xlabel("B")
ax.set_ylabel("J_e_c")
    
for I_c_spec in c:
    y = []
    for B in x:
        J_e_c = get_L_TPL2100(T, B, theta) * I_c_spec/(1.7 * gn.A_tape * 1e6)
        y.append(J_e_c)
        
    ax.plot(x, y, label=f"{I_c_spec = } [A]")
ax.legend()
ax.plot([x[0], x[-1]], [11.8, 11.8], color="black")


#%% vary both B and I_c_spec -- 30K

T = 30 # [K]
B = 4.5 # [T]
theta = 30 # [°]
I_c_spec = 650 # [A]

x = np.linspace(2, 4.5, 100)
c = np.linspace(300, 650, 5)
y = []

fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")
ax.set_xlabel("B")
ax.set_ylabel("J_e_c")
    
for I_c_spec in c:
    y = []
    for B in x:
        J_e_c = get_L_TPL2100(T, B, theta) * I_c_spec/(1.7 * gn.A_tape * 1e6)
        y.append(J_e_c)
        
    ax.plot(x, y, label=f"{I_c_spec = } [A]")
ax.legend()
ax.plot([x[0], x[-1]], [158.8, 158.8], color="black")