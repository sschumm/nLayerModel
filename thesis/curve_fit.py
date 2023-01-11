# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import pandas as pd
import pathlib as pl
import tikzplotlib as tkz


#%%
f = np.array([50, 60, 200, 400, 500])

B_10 = np.array([1.74, 2.17, 13.94, 44.72, 66.26])
B_15 = np.array([3.88, 4.99, 36.67, 129.17, 195.33])


def P_fe(f, p_hy, p_ft):
    return p_hy * f/50 + p_ft * (f/50)**2


popt_10, pcov_10 = curve_fit(P_fe, f, B_10)
popt_15, pcov_15 = curve_fit(P_fe, f, B_15)


def res(f, B, popt):
    plt.figure(dpi=300)
    plt.scatter(f, B)
    rng = np.linspace(0, np.max(f), 1000)
    plt.plot(rng, P_fe(rng, p_hy = popt[0], p_ft = popt[1]), "red")
    plt.grid()
    # plt.show()
    
    tkz.clean_figure()
    tkz.save("aus_mag_losses.tex")
    
    print(f"{popt = }")

res(f, B_10, popt_10)
# res(f, B_15, popt_15)

#%%

directory = "1206_0720_data_77K_p30_r2_plots_second_guess"

filename = "B2_iron_losses"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["f", "B2"]
df = pd.read_csv(path, usecols=columns) 

#%%

gamma_yoke = 7760 # kg/m3
l_e = 0.94
r_sa = 2
r_si = r_sa - 0.0481
V_sy = pi * (r_sa**2 - r_si**2) * l_e

r_ra = r_si - 0.044 - 0.04 - 0.029
r_ri = r_ra - 0.0895
V_ry = pi * (r_ra**2 - r_ri**2) * l_e

V = V_sy + V_ry

m_sy = gamma_yoke * V

p = 30
k_Vy = 1.5

p_hy, p_ft = popt_10[0], popt_10[1]

f_arr = np.array(df.f)
B_arr = np.array(df.B2)


#%%

kf_arr = P_fe(f=f_arr, p_hy=p_hy, p_ft=p_ft)
B_arr_total = B_arr * p * l_e

P_Fey_arr = k_Vy * B_arr_total * kf_arr * gamma_yoke


P_total = np.sum(P_Fey_arr)
print(f"P_Fe = {P_total*1e-3} [kW]")

#%%
plt.figure(dpi=300)
plt.bar(f_arr, P_Fey_arr)
plt.show()

plt.figure(dpi=300)
plt.bar(f_arr[2:], P_Fey_arr[2:], color="red")
plt.show()






#%%


export = True


directory = "fit_B"

filename = "m_s"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")

columns = ["p", "m"]
df = pd.read_csv(path, usecols=columns) # df = dataframe


# fig = plt.figure(dpi=300, figsize=(10,7))
# ax1 = plt.subplot()
# plt.grid()

# ax1.set_xlabel("p")
# ax1.set_ylabel("m")
# ax1.scatter(df.p, df.m, color = "blue", label="m")


p_arr = np.array(df.p)
m_arr = np.array(df.m)

def m_of_p(p, x):
    return x / p

popt_pm, pcov_pm = curve_fit(m_of_p, p_arr, m_arr)

#%%
fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("p")
ax1.set_ylabel("m")
ax1.scatter(df.p, df.m, color = "blue", label="m")
p_extra = np.linspace(4, 50,100) 
# ax1.plot(p_extra, m_of_p(p_extra, popt_pm[0]), color = "red")
ax1.plot(p_arr, m_of_p(p_arr, popt_pm[0]), color = "red")

if export:
    tkz.clean_figure()
    tkz.save("num_m_of_p_fit_stator.tex")


#%%

filename = "m_r"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")
columns = ["p", "m"]
df = pd.read_csv(path, usecols=columns) # df = dataframe
p_arr = np.array(df.p)
mry_arr = np.array(df.m)

filename = "B_rc"
path = pl.Path(__file__).parent / "visualize" / directory / (filename + ".csv")
columns = ["p", "B"]
df = pd.read_csv(path, usecols=columns) # df = dataframe
By_arr = np.array(df.B)

X = np.array([df.p, df.B])

def m_of_p(p, x):
    return x / p

def fit_pB(X, x):
    p, B = X
    return x * B/p

# popt_pm_ry, _ = curve_fit(m_of_p, p_arr, mry_arr)
# popt_pB_rc, _ = curve_fit(m_of_B, p_arr, mry_arr)
popt_pB, _ = curve_fit(fit_pB, X, mry_arr)

fig = plt.figure(dpi=300, figsize=(10,7))
ax1 = plt.subplot()
plt.grid()

ax1.set_xlabel("p")
ax1.set_ylabel("m")
ax1.scatter(p_arr, mry_arr, color = "blue", label="m")
# ax1.scatter(p_arr, By_arr, color = "green", label="m")
p_extra = np.linspace(4, 50,100) 
B_extra = np.linspace(0, 3, 100)
ax1.plot(p_arr, fit_pB((p_arr, By_arr), popt_pB[0]), color = "red")
# ax1.plot(p_extra, m_of_B(p_extra, popt_pB_rc[0]), color = "yellow")
# ax1.plot(p_extra, m_of_p(p_extra, popt_pm_ry[0]) + m_of_B(p_extra, popt_pB_rc[0]), color = "purple")



if export:
    tkz.clean_figure()
    tkz.save("num_m_of_p_fit_rotor.tex")














