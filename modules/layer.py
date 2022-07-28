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
        self.mu = mu_0 * mu_r
        self.mu_inv = 1 / self.mu
        
        self.idx = None
        
        
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
    
    
    def Az(self, p, R, T, a_j, b_j, alpha = 0.0):    
        return np.sin(p*T + alpha) * self.A(p, R, a_j, b_j)
    
    
    # ------ computation of flux densities ------
    
    def Br(self, p, R, T, a_j, b_j, alpha = 0.0):
        if not np.all(R <= self.r):
            print("Warning: computed field lies outside of its corresponding layer")
            
        if self.idx == 0:
            return np.cos(p*T + alpha) * p * (a_j * R**(p-1))
        else: 
            return np.cos(p*T + alpha) * p * (a_j * R**(p-1) + b_j * R**-(p+1))
        
    
    def Bt(self, p, R, T, a_j, b_j, alpha = 0.0):           
        return np.sin(p*T + alpha) * -self.dA(p, R, a_j, b_j)
    
    
    # ------ computation of field strength ------
    
    def Hr(self, p, R, T, a_j, b_j, alpha = 0.0):
        return self.mu_inv * self.Br(p, R, T, a_j, b_j, alpha)
    
    
    def Ht(self, p, R, T, a_j, b_j, alpha = 0.0):
        return self.mu_inv * self.Bt(p, R, T, a_j, b_j, alpha)
        




class CurrentLoading(Layer):
    
    def __init__(self, K: float, r: float, alpha: float = 0., mu_r: float = 1.0):
        super().__init__(r, mu_r)
        
        self.K = K
        self.alpha = alpha

        
class AirLayer(Layer):

    def __init__(self, r: float):
        super().__init__(r, mu_r=1.0)
        
        
class MagneticLayer(Layer):
    
    def __init__(self, r: float, mu_r: float):
        super().__init__(r, mu_r)
        

class Environment(Layer):

    def __init__(self):
        super().__init__(r=np.inf, mu_r=1)        
    