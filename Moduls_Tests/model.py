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
        
        layer.add_index(self.n_layers)
        layer.sync_pole_pairs(self.p)
        self.layers[layer.index] = layer
        
        self.n_layers += 1
                
        
    def build(self):
        
        if self.n_layers == 0:
            raise ValueError("Model contains no layers. Can't build the model.")
        
        for i in range(self.n_layers-1):
            self.layers[i].sync_layer_info(self.layers[i+1].mu_inv)
        
        n_vars = self.n_layers * 2 + 2
        
        self.M = np.zeros((n_vars, n_vars))
        self.y = np.zeros(n_vars)
        
        self.M[ 0,  1] = 1
        self.M[-1, -2] = 1
        
        
        for i in range(self.n_layers):
            
            # j = i+1
            # r = self.layers[i].r
            # mu_inv1 = self.layers[i].mu_inv
            
# =============================================================================
#             if j == self.n_layers:
#                 mu_inv2 = 1 / mu_0
#             else:
#                 mu_inv2 = self.layers[j].mu_inv
# =============================================================================
        
            # B_r = const.
            self.M = self.layers[i].apply_boundaries(self.M, i)
# =============================================================================
#             self.M[i+1, i:i+4] = np.array([-r**self.p, -r**-self.p, 
#                                             r**self.p,  r**-self.p])
# =============================================================================
            # H_0 = const.
            
# =============================================================================
#             self.M[i+2, i:i+4] = np.array([-mu_inv1 * r**(self.p-1),  mu_inv1 * r**-(self.p+1),
#                                             mu_inv2 * r**(self.p-1), -mu_inv2 * r**-(self.p+1)]) * self.p
# =============================================================================
            
        
            
        
        
        
        
        
    
    
    
    