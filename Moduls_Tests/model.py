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
        
        n_bounds = sum([len(vals.bounds) for vals in self.layers.values()])
        n_vars = n_bounds + 2
        
        self.M = np.zeros((n_vars, n_vars))
        self.y = np.zeros(n_vars)
        
        self.M[ 0,  1] = 1
        self.M[-1, -2] = 1
        
        
        i = 0
        for layer in self.layers.values():
                      
            self.M = layer.apply_boundaries(self.M, i)
            i += len(layer.bounds)
        

        
            
        
        
        
        
        
    
    
    
    