# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:34:28 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0



class Layer():
    
    
    def __init__(self, r: float, mu_r: float):
        
        self.r = r
        self.mu_r = mu_r
        self.mu = mu_0 * self.mu_r
        self.mu_inv = 1 / self.mu
        
        self.idx = None
        self.alpha = 0
        
        
    def set_index(self, index: int):
        self.idx = index
        
        
    # ------ computation of vector potential ------
    
    def A(self, p, R, a_j, b_j):
        if not np.all(R <= self.r):
            print("Warning: computed field lies outside of its corresponding layer")
        
        if self.idx == 0:
            return a_j * R**p
        else:
            return a_j * R**p  +  b_j * R**-p
        
    
    def dA(self, p, R, a_j, b_j):
        if not np.all(R <= self.r):
            print("Warning: computed field lies outside of its corresponding layer")
            
        if self.idx == 0:
            return a_j * p * R**(p-1)
        else:
            return a_j * p * R**(p-1)  -  b_j * p * R**-(p+1)
    
    
    def Az(self, p, R, T, a_j, b_j):    
        return np.sin(p*T + self.alpha) * self.A(p, R, a_j, b_j)
    
    
    # ------ computation of flux densities ------
    
    def Br(self, p, R, T, a_j, b_j):
        if not np.all(R <= self.r):
            print("Warning: computed field lies outside of its corresponding layer")
            
        if self.idx == 0:
            return np.cos(p*T + self.alpha) * p * (a_j * R**(p-1))
        else: 
            return np.cos(p*T + self.alpha) * p * (a_j * R**(p-1) + b_j * R**-(p+1))
        
    
    def Bt(self, p, R, T, a_j, b_j):           
        return np.sin(p*T + self.alpha) * -self.dA(p, R, a_j, b_j)
    
    
    # ------ computation of field strength ------
    
    def Hr(self, p, R, T, a_j, b_j):
        return self.mu_inv * self.Br(p, R, T, a_j, b_j)
    
    
    def Ht(self, p, R, T, a_j, b_j):
        return self.mu_inv * self.Bt(p, R, T, a_j, b_j)
        




class CurrentLoading(Layer):
    
    def __init__(self, K: float, r: float, alpha: float = 0., mu_r: float = 1.0):
        super().__init__(r, mu_r)
        
        self.K = K
        self.alpha = alpha
        self.tangential_force = None
        
    
    def Kt(self, p, T):
        # positive angle -> current loading leading clockwise
        # negative angle -> current loading lagging clockwise
        return self.K * np.sin(p * T + self.alpha)

        
class AirLayer(Layer):

    def __init__(self, r: float):
        super().__init__(r, mu_r=1.0)
        
        
class MagneticLayer(Layer):
    
    def __init__(self, r: float, mu_r: float):
        super().__init__(r, mu_r)
        

class Environment(Layer):

    def __init__(self, alpha: float = 0):
        super().__init__(r=np.inf, mu_r=1)
        self.alpha = alpha
    