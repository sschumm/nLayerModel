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
        self.mu_i = list()
        self.mu_o = list()
        
        
    def add_layer(self, layer: Layer):
        self.layers.append(layer)
        
    def build_model(self):
        pass