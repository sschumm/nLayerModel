# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:12:38 2022

@author: svens
"""
from .layer import Layer



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
                
        pass
        
        
    
    
    
    