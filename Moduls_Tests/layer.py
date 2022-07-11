# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:58:47 2022

@author: svens
"""
import numpy as np



class Layer():
    """Header
    
    Description
    
    Args:
        
        
    Attributes:
        
        
    """
    
    mu_0 = 4 * np.pi * 10**(-7) # [N / A**2]
    
    def __init__(self, radius: float = None, mu: float = 1.0):
        
        if mu is None:
            mu = 1.0
            
        if isinstance(mu, int):
            mu = float(mu)
            
        if isinstance(radius, int):
            radius = float(radius)
            
        self.radius = radius
        self.mu = mu
            
    
    
    
    
    
    
    
    