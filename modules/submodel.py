# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:54:41 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0

from .continuities import c_Br, c_Ht
from .layer import CurrentLoading

class SubModel():
    
    def __init__(self, p: int, layers: list, alpha: float):
        
        self.p = p
        self.alpha = alpha
        
        self.sysA = None
        self.sysx = None
        self.sysb = None
        
        # assume that the list of layers is already reordered by radius
        # assume that all layers got their indices
        
        self.layers = layers        
        
        # extract and store the layer data 
        self.radii = [i.r for i in self.layers]
        self.mu_i_inv = [i.mu_inv for i in self.layers]
        self.mu_o_inv = self.mu_i_inv[1:] + [1 / mu_0]
        
        self.x = None
        
        self._build_sysA()
        self._build_sysb()
        
        
    def _build_sysA(self):
        
        sysA = []
        n = (len(self.layers) - 1)
        
        # add inner boundary condition
        sysA.append([0, 1, 0, 0] + [0]*2*n)
        
        # add continuity conditions
        for i, layer in enumerate(self.layers):
            Br = c_Br(r = self.radii[i], 
                      p = self.p)
            Ht = c_Ht(r = self.radii[i], 
                      p = self.p, 
                      mu_i_inv=self.mu_i_inv[i], 
                      mu_o_inv = self.mu_o_inv[i])
            
            sysA.append([0]*2*i + Br + [0]*2*(n-i))
            sysA.append([0]*2*i + Ht + [0]*2*(n-i))
        
        # add outer boundary condition
        sysA.append([0]*2*n + [0, 0, 1, 0])
        
        self.sysA = np.asarray(sysA)

        
    def _build_sysb(self):
        
        sysb = [0.]
        
        for layer in self.layers:
            if isinstance(layer, CurrentLoading):
                sysb.append(0.)
                sysb.append(layer.K)
            else:
                sysb.append(0.)
                sysb.append(0.)
        sysb.append(0.)
        self.sysb = np.asarray(sysb)
            
        
    def solve(self):
        x = np.linalg.solve(self.sysA, self.sysb)
        
        if np.allclose(np.dot(self.sysA, x), self.sysb):
            self.x = x
            return x
        else:
            raise Exception("System solution x is NOT np.allclose().")




