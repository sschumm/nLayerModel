# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:12:38 2022

@author: svens
"""
import numpy as np

from .layer import Layer, CurrentLoading

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
        self.built = False
         
    
    def add_layer(self, layer: Layer):
        
        layer.add_index(self.n_layers)
        layer.sync_pole_pairs(self.p)
        self.layers[layer.index] = layer
        
        self.n_layers += 1
                
        
    def build(self):
        
        if self.n_layers == 0:
            raise Exception("Model contains no layers. Use add_layer() before trying to build the model.")
        
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
            
            if isinstance(layer,CurrentLoading):
                self.y = layer.apply_solution(self.y, i)
            
            self.M = layer.apply_boundaries(self.M, i)
            i += len(layer.bounds)
            
        self.built = True
            
            
    def solve(self, allclose_check: bool= False):
        if self.built:    
            x = np.linalg.solve(a = self.M, b = self.y)
            
            if allclose_check:
                return x, np.allclose(np.dot(self.M, x), self.y)    
            else:
                return x
        else:
            raise Exception("Model has not been built yet. Use build() before trying to solve the model.")
        
        
    def vary_param(self, layer: Layer, param: str, val: float):
        if self.built:    
           if hasattr(layer, param):
               exec(f"layer.{param} = {val}")
               layer.update_boundaries()
               if isinstance(layer,CurrentLoading):
                   self.y = layer.apply_solution(self.y, layer.index)
                   
               self.M = layer.apply_boundaries(self.M, layer.index)
               

               print(self.M)
               print(self.y)
               print(self.solve())
               # return self.solve()
           else:
               raise Exception(f"The layer '{layer}' has no parameter '{param}'.")
        else:
            raise Exception("Model has not been built yet. Use build() before trying to solve the model.")
        
        

        
            
        
        
        
        
        
    
    
    
    