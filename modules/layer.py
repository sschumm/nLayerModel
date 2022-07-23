# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:34:28 2022

@author: svens
"""

# import numpy as np

from scipy.constants import mu_0



class Layer():
    
    
    def __init__(self, r, mu_r):
        
        self.r = r
        self.mu = mu_0 * mu_r
        self.mu_inv = 1 / self.mu
        
        self.idx = None
        
        
    def set_index(self, index):
        self.idx = index
        
        

        
    