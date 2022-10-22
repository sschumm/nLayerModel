# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from modules import Model, CurrentLoading
from analytics.two_regions import analytic_solution

np.set_printoptions(suppress=True, linewidth=250, precision=40)


def f(p, r, K, atol, rtol):
    model = Model(p)
    model.add_layer(CurrentLoading(K, r, mu_r=1))
    model.build()
    
    x_numeric = model.solve()[0]
    x_analytic= analytic_solution(p, r, K)
    
    a1_n, b1_n, a2_n, b2_n = x_numeric
    a1_a, b1_a, a2_a, b2_a = x_analytic
    
    close = np.allclose(x_numeric, x_analytic, rtol=rtol, atol=atol)
    
    #return close, x_numeric, x_analytic
    return (b2_n, b2_a)


#%% 

rtol = 0
atol = 1e-15


x = np.arange(1,100)
y1 = list()
y2 = list()

for p in x:
    n, a = f(p=p, r=2., K=1e6, rtol=rtol, atol=atol)
    y1.append(np.abs(n - a))

for p in x:
    n, a = f(p=p, r=1.5, K=1e6, rtol=rtol, atol=atol)
    y2.append(np.abs(n - a))

fig = plt.figure(dpi=300)
ax = plt.subplot()
# ax.plot(x, y)
# ax.semilogy(x, y)
markersize = 5
ax.scatter(x, y1, s=markersize)
ax.scatter(x, y2, s=markersize)
ax.set_yscale("log")
