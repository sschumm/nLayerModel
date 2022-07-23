# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:41:40 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0
from .continuities import c_Br, c_Ht
from .layer import Layer, CurrentLoading


class Model():
    
    def __init__(self, p: int):
        
        self.p = p
        
        
        self.sysA = None
        self.sysx = None
        self.sysb = None
        
        self.layers = list()
        
        self.radii = list()
        self.mu_i_inv = list()
        self.mu_o_inv = list()
        
        self.x = None
        
        
    def add_layer(self, layer: Layer):
        self.layers.append(layer)

        
    def build(self):
        # reorder the list of layers by radius
        self.layers.sort(key= lambda lay: lay.r)
        
        # extract and store the layer data 
        self.radii = [i.r for i in self.layers]
        self.mu_i_inv = [i.mu_inv for i in self.layers]
        self.mu_o_inv = self.mu_i_inv[1:] + [1 / mu_0]
        
        # set layer indices
        for idx, layer in enumerate(self.layers):
            layer.idx = idx
        
        self._build_sysA()
        self._build_sysb()
        
    
    def solve(self):
        x = np.linalg.solve(self.sysA, self.sysb)
        
        if np.allclose(np.dot(self.sysA, x), self.sysb):
            self.x = x
            return x
        else:
            raise Exception("(sschumm) solution seems to be no good...")
        
        
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
            

        