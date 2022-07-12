# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 07:20:03 2022

@author: svens
"""
import numpy as np

mu_0 = 4 * np.pi * 10**(-7) # [N / A**2]

class Boundary():
    
    def __init__(self, r: float, p: int = 1, mu_inv: float = 1 / mu_0):
        
        if p is None:
            p = 1
        
        if mu_inv is None:
            mu_inv = 1 / mu_0
            
        if isinstance(mu_inv, int):
            mu_inv = float(mu_inv)
        
        self.p = p
        self.r = r
        self.mu_inv = mu_inv
        self.mu_inv_out = self.mu_inv
    
    
    def get_boundary(self):
        return np.ones() * 999.0 # for tests
         
    
    def sync_pole_pairs(self, n_pole_pairs: int):
        self.p = n_pole_pairs

    
    
class B_r(Boundary):

    def __init__(self, r: float, p: int = 1):
        super().__init__(r, p)
        
    
    def get_boundary(self):
        return np.array([-self.r** self.p, 
                         -self.r**-self.p, 
                          self.r** self.p,  
                          self.r**-self.p])
    

class B_0(Boundary):

    def __init__(self, r: float, p: int = 1):
        super().__init__(r, p)        


class H_r(Boundary):

    def __init__(self, r: float, p: int = 1):
        super().__init__(r, p)
        
        
class H_0(Boundary):
    
    def __init__(self, r: float, p: int = 1, mu_inv: float = 1.0):
        super().__init__(r, p, mu_inv)
        
        
    def set_outer_mu_inv(self, mu_inv_out: float):
        self.mu_inv_out = mu_inv_out
        
        
    def get_boundary(self):
        
        return np.array([-self.mu_inv     * self.r** (self.p-1), 
                          self.mu_inv     * self.r**-(self.p+1),
                          self.mu_inv_out * self.r** (self.p-1), 
                         -self.mu_inv_out * self.r**-(self.p+1)]) * self.p
        
