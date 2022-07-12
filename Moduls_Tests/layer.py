# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:58:47 2022

@author: svens
"""
import numpy as np

from .boundary import Boundary, B_r, H_0

mu_0 = 4 * np.pi * 10**(-7) # [N / A**2]

class Layer():
    """Header
    
    Description
    
    Args:
        
        
    Attributes:
        
        
    """
    
   
    
    def __init__(self, r: float, mu_r: float = 1.0, p: int = 1):
        
        if mu_r is None:
            mu_r = 1.0
            
        if p is None:
            p = 1
        
        if isinstance(mu_r, int):
            mu_r = float(mu_r)
            
        if isinstance(r, int):
            r = float(r)
        
        self.r = r
        self.mu_r = mu_r
        self.mu = mu_0 * self.mu_r
        self.mu_inv = 1 / self.mu
        self.mu_inv_outer = self.mu_inv
        self.p = p
        
        
        self.bounds = []
        
        self.add_boundary(B_r(r = self.r, p = self.p))
        self.add_boundary(H_0(r = self.r, p = self.p, mu_inv = self.mu_inv))
        
        
    def add_index(self, index: int):
        
        self.index = index
    
    
    def add_boundary(self, boundary: Boundary):
        self.bounds.append(boundary)
        
        
    def apply_boundaries(self, matrix: np.ndarray, index: int):
               
        for a, bound in enumerate(self.bounds):
            matrix[index+a+1, index:index+4] = bound.get_boundary()
        return matrix    
        
        
    def sync_pole_pairs(self, n_pole_pairs: int):
        self.p = n_pole_pairs
        for bound in self.bounds:
            bound.sync_pole_pairs(self.p)
            
            
    def sync_layer_info(self, mu_inv_out: float):
        self.mu_inv_outer = mu_inv_out
        for bound in self.bounds:
            if bound is H_0:
                bound.set_outer_mu_inv(mu_inv_out)
        
        
        
    
    
    
    
    
    
    
    
    
    
        
        
    