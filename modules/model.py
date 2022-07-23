# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:41:40 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0
from .layer import Layer


class Model():
    
    def __init__(self, p: int):
        
        self.p = p
        
        
        self.sysA = None
        self.sysx = None
        self.sysb = None
        
        self.layers = list()
        
        self.radii = list()
        self.mu_i_inv = list()
        self.mu_o_inv = list()
        
        
    def add_layer(self, layer: Layer):
        self.layers.append(layer)
        
    def build(self):
        # reorder the list of layers by radius
        self.layers.sort(key= lambda lay: lay.r)
        # extract and store the layer data 
        self.radii = [i.r for i in self.layers]
        self.mu_i_inv = [i.mu_inv for i in self.layers]
        self.mu_o_inv = self.mu_i_inv[1:] + [1 / mu_0]
        