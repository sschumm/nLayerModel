# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:54:41 2022

@author: svens
"""

import numpy as np

from .continuities import c_Br, c_Ht
from .layer import CurrentLoading

class SubModel():
    
    def __init__(self, p: int):
        
        self.p = p
        
        self.sysA = None
        self.sysx = None
        self.sysb = None
        
        self.layers = list()
        self.current_loadings = list()
        
        self.radii = list()
        self.mu_i_inv = list()
        self.mu_o_inv = list()
        
        self.x = None
        
        
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
            




