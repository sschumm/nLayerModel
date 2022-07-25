# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 14:41:40 2022

@author: svens
"""

import numpy as np

from scipy.constants import mu_0
from .continuities import c_Br, c_Ht
from .layer import Layer, CurrentLoading, Environment
from .utils import rt_to_xy, BrBt_to_UV


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
            
            
    def get_B_data(self, r, t):
        R_tuple, T_tuple = tuple(), tuple()
        Br_tuple, Bt_tuple = tuple(), tuple()
        
        plot_layers = self.layers + [Environment()]
        for i, j in enumerate(range(0, len(self.x), 2)):
            if i == 0:
                r_i = 0.
            else:
                r_i = self.layers[i - 1].r
            
            this_r = r[np.argwhere((r >= r_i) & (r < plot_layers[i].r)).flatten()]
            this_R, this_T = np.meshgrid(this_r, t)
            
            this_Br= plot_layers[i].Br(self.p, this_R, this_T,
                                       a_j = self.x[j],
                                       b_j = self.x[j+1])
            this_Bt= plot_layers[i].Bt(self.p, this_R, this_T,
                                       a_j = self.x[j],
                                       b_j = self.x[j+1])

            R_tuple += (this_R, )
            T_tuple += (this_T, )
            Br_tuple += (this_Br, )
            Bt_tuple += (this_Bt, )

        R, T = np.hstack(R_tuple), np.hstack(T_tuple)
        Br, Bt = np.hstack(Br_tuple), np.hstack(Bt_tuple)
        
        X, Y = rt_to_xy(R, T)
        U, V = BrBt_to_UV(Br, Bt, T)
        return X, Y, U, V
    
    
    def get_A_xy_data(self, r, t):
        R_tuple, T_tuple = tuple(), tuple()
        Az_tuple = tuple()
        
        # this does not add a new layer to the model but is used to compute
        # the field for the environment, otherwise a_n & b_n would be unused
        plot_layers = self.layers + [Environment()]
        for i, j in enumerate(range(0, len(self.x), 2)):
            if i == 0:
                r_i = 0.
            else:
                r_i = self.layers[i - 1].r
            
            # create mesh for every layer
            this_r = r[np.argwhere((r >= r_i) & (r < plot_layers[i].r)).flatten()]
            this_R, this_T = np.meshgrid(this_r, t)
            
            # compute the vector potential for the i-th layer
            this_Az= plot_layers[i].Az(self.p, this_R, this_T, 
                                       a_j = self.x[j], 
                                       b_j = self.x[j+1])
            
            # store the result for the i-th layer
            R_tuple += (this_R, )
            T_tuple += (this_T, )
            Az_tuple += (this_Az, )
        
        R, T = np.hstack(R_tuple), np.hstack(T_tuple)
        Az = np.hstack(Az_tuple)
        X, Y = rt_to_xy(R, T)
        return X, Y, Az
            
        
        
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
            

        