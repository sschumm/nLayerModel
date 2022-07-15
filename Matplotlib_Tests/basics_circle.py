# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 22:12:44 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 
 
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
 
radius = 0.4
 
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 
 
figure, axes = plt.subplots( 1) 
 
axes.plot( x, y ) 
axes.set_aspect( 1 ) 

radius = 0.2
 
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 

axes.plot( x, y ) 
axes.set_aspect( 1 )
 
plt.title( 'Parametric Equation Circle' ) 
plt.show() 

#%%

import numpy as np 
import matplotlib.pyplot as plt 

cc = plt.Circle((0, 0), 0.4, fill=False, color="b", linestyle="-.")

figure, axes = plt.subplots( 1) 
lim = 0.5
axes.add_artist(cc)
axes.set_xlim(-lim, lim)
axes.set_ylim(-lim, lim)
axes.set_aspect( 1 )

plt.show()

