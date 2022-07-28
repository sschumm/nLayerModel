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
                    new_cl = CurrentLoading(K=0., r=lay.r, 
                                            alpha=0, mu_r=lay.mu_r)
                    
                    new_cl.set_index(lay.idx)
                    layers[idx] = new_cl
            new_submodel = SubModel(self.p, layers)
            self.submodels.append(new_submodel)                   
           
            
    def solve(self):
        
        x = []
        for subm in self.submodels:
            x.append(subm.solve())
        
        self.x = np.asarray(x)
        return x


    def get_A_data(self, r, t):
        R_tuple, T_tuple = tuple(), tuple()
        Az_tuple = tuple()
        
        
        # computes the vector potential for the i-th layer per iteration
        for i, j in enumerate(range(0, self.x.shape[-1], 2)):
            if i == 0:
                r_i = 0.
            else:
                r_i = self.layers[i - 1].r
                
            if j == (self.x.shape[-1] - 2):
                r_a = np.inf
            else:
                r_a = self.layers[i].r
            
            # create mesh for the i-th layer
            this_r = r[np.argwhere((r >= r_i) & (r < r_a)).flatten()]
            this_R, this_T = np.meshgrid(this_r, t)
            
            
            # start computing the superposition of the vector potential
            this_Az = np.zeros(this_R.shape)
            
            for subm in self.submodels:
                
                # this does not add a new layer to the model but is used to compute
                # the field for the environment, otherwise a_n & b_n would be unused
                plot_layers = subm.layers + [Environment()]
                
                # sums up the vector potential for all current loadings
                this_Az += plot_layers[i].Az(self.p, this_R, this_T, 
                                             a_j = subm.x[j], 
                                             b_j = subm.x[j+1])
            
 
            # store the result for the i-th layer
            R_tuple += (this_R, )
            T_tuple += (this_T, )
            Az_tuple += (this_Az, )
        
        R, T = np.hstack(R_tuple), np.hstack(T_tuple)
        Az = np.hstack(Az_tuple)
        X, Y = rt_to_xy(R, T)
        # data = (X, Y, Az, R, Ts)
        return X, Y, Az, R, T















        