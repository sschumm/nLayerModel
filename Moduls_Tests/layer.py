# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:58:47 2022

@author: svens
"""
import numpy as np

mu_0 = 4 * np.pi * 10**(-7) # [N / A**2]

class Layer():
    """Header
    
    Description
    
    Args:
        
        
    Attributes:
        
        
    """
    
   
    
    def __init__(self, r: float, mu_r: float = 1.0):
        
        if mu_r is None:
            mu_r = 1.0
            
        if isinstance(mu_r, int):
            mu_r = float(mu_r)
            
        if isinstance(r, int):
            r = float(r)
            
        self.r = r
        self.mu_r = mu_r
        self.mu = mu_0 * self.mu_r
        self.mu_inv = 1 / self.mu
        # self.index = None
        
        
    def add_index(self, index: int):
        
        self.index = index
    
    

    