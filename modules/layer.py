# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:34:28 2022

@author: svens
"""

# import numpy as np

from scipy.constants import mu_0



class Layer():
    
    
    def __init__(self, r: float, mu_r: float):
        
        self.r = r
        self.mu = mu_0 * mu_r
        self.mu_inv = 1 / self.mu
        
        self.idx = None
        
        
    def set_index(self, index: int):
        self.idx = index
        
        

class CurrentLoading(Layer):
    
    def __init__(self, K: float, r: float, mu_r: float):
        super().__init__(r, mu_r)
        
        self.K = K
        
class AirLayer(Layer):
    def __init__(self, r: float):
        super().__init__(r, mu_r=1)
        
    