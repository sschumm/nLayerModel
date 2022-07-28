# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:53:17 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0
from .continuities import c_Br, c_Ht
from .layer import Layer, CurrentLoading, Environment
from .utils import rt_to_xy, BrBt_to_UV
from .submodel import SubModel


class Model():
    def __init__(self, p: int):
        
        self.p = p
        
        self.layers = list()
        self.current_loadings = list()
        
        self.submodels = list()
        self.x = None
        
        
        
    def add_layer(self, layer: Layer):
        self.layers.append(layer)

        
    def build(self):
        # reorder the list of layers by radius
        self.layers.sort(key= lambda lay: lay.r)
        
        # add all current loading layers to a separate list
        for layer in self.layers:
            if isinstance(layer, CurrentLoading):
                self.current_loadings.append(layer)
        
        # extract and store the layer data 
        self.radii = [i.r for i in self.layers]
        self.mu_i_inv = [i.mu_inv for i in self.layers]
        self.mu_o_inv = self.mu_i_inv[1:] + [1 / mu_0]
        
        # set layer indices
        for idx, layer in enumerate(self.layers):
            layer.idx = idx
            
        # create a submodel for each current loading
        for idx, layer in enumerate(self.current_loadings):
            layers = self.layers.copy()
            
            for idx, lay in enumerate(layers):
                if isinstance(lay, CurrentLoading) and (lay.idx != layer.idx):
                    layers[idx] = CurrentLoading(K=0., r=lay.r, 
                                                 alpha=lay.alpha, mu_r=lay.mu_r)
            new_submodel = SubModel(self.p, layers)
            self.submodels.append(new_submodel)
            
        sysA = []
        sysb = []
        
        for subm in self.submodels:
            sysA.append(subm.sysA)
            sysb.append(subm.sysb)
                    
           
            
    def solve(self):
        
        x = []
        for subm in self.submodels:
            x.append(subm.solve())
        
        self.x = np.asarray(x)
        return x


















        