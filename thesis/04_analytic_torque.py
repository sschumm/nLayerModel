# -*- coding: utf-8 -*-
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as plt
from scipy.constants import pi
from modules import Model, CurrentLoading, MagneticLayer

from analytics.four_regions import analytic_solution_K1, analytic_solution_K2
from analytics.four_regions import analytic_torque_on_K1

np.set_printoptions(suppress=True, linewidth=250, precision=40)


def f(p, r_1, r_2, r_3, K_r, K_s, alpha_r, alpha_s):
    l=2
    model = Model(p, l)
    model.add_layer(CurrentLoading(K=K_r, r=r_1, mu_r=1e5, alpha=alpha_r))
    model.add_layer(CurrentLoading(K=K_s, r=r_2, mu_r=1, alpha=alpha_s))
    model.add_layer(MagneticLayer(r=r_3, mu_r=1e5))
    model.build()
    model.solve()    
    model.total_torque()
    
    x_analytic = analytic_solution_K2(p, K2=K_s, r1=r_1, r2=r_2, r3=r_3,  
                              mu_r1=1e5, mu_r3=1e5)    
        
    
    M_numeric = model.Mpos
    M_analytic = analytic_torque_on_K1(p, l=2, K1=K_r, r1=r_1, 
                                       a1_K2=x_analytic[0],
                                       alpha1=alpha_r, alpha2=alpha_s)
    
    return M_numeric, M_analytic
    
    
#%% test

f(2, 1., 2., 3., 2e6, 1e6, 0, pi/2)


#%% random
x = [i for i in range(1000)]
y = []
for i in x:
    p = np.random.randint(1, 100)
    r_1 = 0.3 * np.random.rand() + 1.2
    r_2 = r_1 + 0.3 * np.random.rand() + 0.2
    r_3 = r_2 + 0.3 * np.random.rand() + 0.2
    K_r = 2e6 + np.random.rand() * 1e6
    K_s = 1e6 + np.random.rand() * 1e6
    alpha_r = np.random.rand() * pi/2
    alpha_s = 0
    
    n, a = f(p, r_1, r_2, r_3, K_r, K_s, alpha_r, alpha_s)
    e = np.abs((n-a)/a)
    y.append(e)
    if e > 0.002:
        print(n,a)
    
    
fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")

markersize = 5
ax.scatter(x, y, s=markersize)

#%% over p

# m = [i for i in range(100)]
x, y = [], []

for p in range(1, 101, 2):
    
    for j in range(20):
        r_1 = 0.3 * np.random.rand() + 1.2
        r_2 = r_1 + 0.3 * np.random.rand() + 0.2
        r_3 = r_2 + 0.3 * np.random.rand() + 0.2
        K_r = 2e6 + np.random.rand() * 1e6
        K_s = 1e6 + np.random.rand() * 1e6
        alpha_r = np.random.rand() * pi/2
        alpha_s = 0
        
        n, a = f(p, r_1, r_2, r_3, K_r, K_s, alpha_r, alpha_s)
        e = np.abs((n-a))#/a)
        x.append(p)
        y.append(e)

#%%
    
fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")

markersize = 5
ax.scatter(x, y, s=markersize)
ax.set_yscale("log")



tkz.clean_figure()
tkz.save("msm_abs_error_torque.tex")


#%%

x, y, dy = [], [], []

for alpha_s in np.linspace(0, pi/2, 20):
    
    data = []
    
    for i in range(10000):
        p = np.random.randint(1, 100)
        r_1 = 0.3 * np.random.rand() + 1.2
        r_2 = r_1 + 0.3 * np.random.rand() + 0.2
        r_3 = r_2 + 0.3 * np.random.rand() + 0.2
        K_r = 2e6 + np.random.rand() * 1e6
        K_s = 1e6 + np.random.rand() * 1e6
        alpha_r = 0
        alpha_s = alpha_s
    
        n, a = f(p, r_1, r_2, r_3, K_r, K_s, alpha_r, alpha_s)
        e = np.abs((n-a))#/a)
        # x.append(alpha_s)
        data.append(e)
        
    median = np.median(data)
    top = np.max(data)
    bot = np.min(data)
    x.append(alpha_s)
    y.append(median)
    dy.append(top-bot)
    

#%%
fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid("both")
# plt.style.use('seaborn-whitegrid')

# markersize = 5
# ax.scatter(x, y, s=markersize)
# ax.set_yscale("log")
ax.errorbar(x, y, yerr=dy, fmt='o', markersize=1.5, linewidth=1, capsize=3)

# tkz.clean_figure()
tkz.save("msm_abs_error_torque.tex")















