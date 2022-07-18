# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 12:26:16 2022

@author: svens
"""

import numpy as np 
import matplotlib.pyplot as plt 



def plot(X, Y, U, V):
    
    # r = np.array([1, 2, 3, 4])
    # c = ["black", "blue", "red", "green"]
    # lim = np.max(r) * 1.2
    # lim = 5
    
    plt.figure(figsize=(10, 10))
    ax = plt.subplot()
    ax.set_aspect( 1 )
    # ax.set_xlim(-lim, lim)
    # ax.set_ylim(-lim, lim)
    
    # =============================================================================
    # for i in range(3, -1, -1):
    #     ax.add_artist(plt.Circle((0, 0), 
    #                              radius = r[i], 
    #                              fill = True, 
    #                              edgecolor = "r",
    #                              facecolor = c[i], 
    #                              linestyle = "-",
    #                              linewidth = 3, 
    #                              alpha=0.1))
    #  
    # =============================================================================
       
        
    """ Adding some streamplot stuff """
    # 1D arrays
    # x = np.arange(-5,5,0.1)
    # y = np.arange(-5,5,0.1)
      
    # Meshgrid
    # X,Y = np.meshgrid(x,y)
      
    # Assign vector directions
    # Ex = (X + 1)/((X+1)**2 + Y**2) - (X - 1)/((X-1)**2 + Y**2)
    # Ey = Y/((X+1)**2 + Y**2) - Y/((X-1)**2 + Y**2)
    
    
    
    ax.streamplot(X, Y, U, V, density=1.4, linewidth=None, color="black")
