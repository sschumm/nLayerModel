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
        self.mu_inv_out = self.mu_inv
        self.p = p
        
        self.bounds = []
        self._init_boundaries()
        

    def _init_boundaries(self):
        self.add_boundary(Boundary(r = self.r, p = self.p))
        self.add_boundary(Boundary(r = self.r, p = self.p))

        
    def add_index(self, index: int):
        
        self.index = index
    
    
    def add_boundary(self, boundary: Boundary):
        self.bounds.append(boundary)
        
        
    def update_boundaries(self):
        self.bounds.clear()
        self._init_boundaries()
        
        
    def apply_boundaries(self, matrix: np.ndarray, index: int):
               
        for a, bound in enumerate(self.bounds):
            matrix[index+a+1, index:index+4] = bound.get_boundary()
        return matrix    
        
        
    def sync_pole_pairs(self, n_pole_pairs: int):
        self.p = n_pole_pairs
        for bound in self.bounds:
            bound.sync_pole_pairs(self.p)
            
            
    def sync_layer_info(self, mu_inv_out: float):
        self.mu_inv_out = mu_inv_out
        for bound in self.bounds:
            if bound is H_0:
                bound.set_outer_mu_inv(mu_inv_out)
        
        
        


class CurrentLoading(Layer):

    
    def __init__(self, r: float, K: float, mu_r: float = 1.0, p: int = 1):
        super().__init__(r, mu_r, p)
        
        if isinstance(K, int):
            r = float(K)
        
        self.K = K
    
    
    def _init_boundaries(self):
        self.add_boundary(B_r(r = self.r, p = self.p))
        self.add_boundary(H_0(r = self.r, p = self.p, mu_inv = self.mu_inv))
        
    
    def apply_solution(self, vector: np.ndarray, index: int):
        
        for a, bound in enumerate(self.bounds):
            if isinstance(bound, H_0):
                vector[index+a+1] = self.K
        return vector 


class MagneticLayer(Layer):

    
    def __init__(self, r: float, mu_r: float, p: int = 1):
        self.mu_r = mu_r
        super().__init__(r, mu_r, p)
        
       
    def _init_boundaries(self):
        self.add_boundary(B_r(r = self.r, p = self.p))
        self.add_boundary(H_0(r = self.r, p = self.p, mu_inv = self.mu_inv))
        
        
class AirLayer(Layer):

    
    def __init__(self, r: float, p: int = 1):
        mu_r = 1.0
        super().__init__(r, mu_r, p)
        
        
    def _init_boundaries(self):
        self.add_boundary(B_r(r = self.r, p = self.p))
        self.add_boundary(H_0(r = self.r, p = self.p, mu_inv = self.mu_inv))
        
    

    
        
        
    