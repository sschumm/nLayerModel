# -*- coding: utf-8 -*-
import numpy as np
import tikzplotlib as tkz
import matplotlib.pyplot as plt
from modules import Model, CurrentLoading
from analytics.two_regions import analytic_solution

np.set_printoptions(suppress=True, linewidth=250, precision=40)


def f(p, r, K, atol, rtol, only_b2=True):
    model = Model(p)
    model.add_layer(CurrentLoading(K, r, mu_r=1))
    model.build()
    
    x_numeric = model.solve()[0]
    x_analytic= analytic_solution(p, r, K)
    
    a1_n, b1_n, a2_n, b2_n = x_numeric
    a1_a, b1_a, a2_a, b2_a = x_analytic
    
    close = np.allclose(x_numeric, x_analytic, rtol=rtol, atol=atol)
    
    if only_b2:
        return b2_n, b2_a
    else:
        return close, x_numeric, x_analytic

rtol = 0
atol = 1e-15

#%% static print

f(p=1, r=2.0, K=5.31e6, rtol=rtol, atol=atol, only_b2=False)

#%% static work

r = 0.425
p = 100
c, x_num, x_an = f(p=p, r=r, K=5.31e6, rtol=rtol, atol=atol, only_b2=False)

A_jr_num = np.array([r**p, r**-p]*2) * x_num
A_jr_an  = np.array([r**p, r**-p]*2) * x_an

print(A_jr_num)
print(A_jr_an)


#%% only b2 plot

x = np.arange(1,100)
y1 = list()
y2 = list()
y3 = list()
y4 = list()

for p in x:
    n, a = f(p=p, r=2., K=5.31e6, rtol=rtol, atol=atol)
    y1.append(np.abs(n - a))

for p in x:
    n, a = f(p=p, r=1.5, K=5.31e6, rtol=rtol, atol=atol)
    y2.append(np.abs(n - a))
    
for p in x:
    n, a = f(p=p, r=0.425, K=5.31e6, rtol=rtol, atol=atol)
    y3.append(np.abs(n - a))

for p in x:
    n, a = f(p=p, r=0.8, K=5.31e6, rtol=rtol, atol=atol)
    y4.append(np.abs(n - a))
    
fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid(which="both")
# ax.plot(x, y)
# ax.semilogy(x, y)
markersize = 5
ax.scatter(x, y1, s=markersize)
ax.scatter(x, y2, s=markersize)
# ax.scatter(x, y3, s=markersize)
ax.scatter(x, y4, s=markersize)
ax.set_yscale("log")
tkz.clean_figure()
tkz.save("msm_abs_error_diverge.tex")

#%%

x = np.arange(1,100)
y1 = list()
y2 = list()
y3 = list()
y4 = list()

for p in x:
    r=2.
    c, x_num, x_an = f(p=p, r=r, K=5.31e6, rtol=rtol, atol=atol, only_b2=False)

    A_jr_num = np.array([r**p, r**-p]*2) * x_num
    A_jr_an  = np.array([r**p, r**-p]*2) * x_an
    
    n = A_jr_num[3]
    a = A_jr_an[3]
    
    y1.append(np.abs(n - a))
    
for p in x:
    r=1.5
    c, x_num, x_an = f(p=p, r=r, K=5.31e6, rtol=rtol, atol=atol, only_b2=False)

    A_jr_num = np.array([r**p, r**-p]*2) * x_num
    A_jr_an  = np.array([r**p, r**-p]*2) * x_an
    
    n = A_jr_num[3]
    a = A_jr_an[3]
    
    y2.append(np.abs(n - a))
    
for p in x:
    r=0.425
    c, x_num, x_an = f(p=p, r=r, K=5.31e6, rtol=rtol, atol=atol, only_b2=False)

    A_jr_num = np.array([r**p, r**-p]*2) * x_num
    A_jr_an  = np.array([r**p, r**-p]*2) * x_an
    
    n = A_jr_num[3]
    a = A_jr_an[3]
    y3.append(np.abs(n - a))

for p in x:
    r=0.8
    c, x_num, x_an = f(p=p, r=r, K=5.31e6, rtol=rtol, atol=atol, only_b2=False)

    A_jr_num = np.array([r**p, r**-p]*2) * x_num
    A_jr_an  = np.array([r**p, r**-p]*2) * x_an
    
    n = A_jr_num[3]
    a = A_jr_an[3]
    
    print(n)
    print(a)
    print()
    y4.append(np.abs(n - a))
    
fig = plt.figure(dpi=300)
ax = plt.subplot()
plt.grid(which="both")
# ax.plot(x, y)
# ax.semilogy(x, y)
markersize = 5
ax.scatter(x, y1, s=markersize)
ax.scatter(x, y2, s=markersize)
# ax.scatter(x, y3, s=markersize)
ax.scatter(x, y4, s=markersize)
ax.set_yscale("log")

tkz.clean_figure()
tkz.save("msm_abs_error_clean.tex")