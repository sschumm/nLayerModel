# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:12:38 2022

@author: svens
"""
import numpy as np

from .layer import Layer

mu_0 = 4 * np.pi * 10**(-7) # [N / A**2]

class Model():
    """Header
    
    Description
    
    Args:
        
        
    Attributes:
        
        
    """
    
    
    def __init__(self, p: int = 1):
        
        if p is None:
            p = 1
        
        self.p = p
        self.n_layers = 0
        self.layers = dict()
         
    
    def add_layer(self, layer: Layer):
        
        self.n_layers += 1
        
        layer.add_index(self.n_layers)
        self.layers[layer.index] = layer
                
        
    def build(self):
        
        n_vars = self.n_layers * 2 + 2
        
        self.M = np.zeros((n_vars, n_vars))
        self.y = np.zeros(n_vars)
        
        self.M[ 0,  1] = 1
        self.M[-1, -2] = 1
        
        
        for i, lay in self.layers.items():
        
            r = lay.r
            mu_inv1 = lay.mu_inv
            
            if i == self.n_layers:
                mu_inv2 = 1 / mu_0
            else:
                mu_inv2 = self.layers[i + 1].mu_inv
        
            # B_r = const.
            self.M[  i, i-1:i+3] = np.array([-r**self.p, -r**-self.p, 
                                              r**self.p,  r**-self.p])
            # H_0 = const.
            self.M[i+1, i-1:i+3] = np.array([-mu_inv1 * r**(self.p-1),  mu_inv1 * r**-(self.p+1),
                                              mu_inv2 * r**(self.p-1), -mu_inv2 * r**-(self.p+1)]) * self.p
            
        
            
        
        
        
        
        
    
    
    
    