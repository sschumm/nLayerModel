# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 11:53:17 2022

@author: svens
"""

import copy
import numpy as np

from scipy.constants import mu_0, pi
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
            layers = copy.deepcopy(self.layers)
            
            for idx, lay in enumerate(layers):
                if isinstance(lay, CurrentLoading) and (lay.idx != layer.idx):
                    lay.alpha = layer.alpha
                    lay.K = 0.
                else:
                    lay.alpha = layer.alpha
                    
            new_submodel = SubModel(self.p, layers)
            self.submodels.append(new_submodel)                   
           
            
    def solve(self):
        
        x = []
        for subm in self.submodels:
            x.append(subm.solve())
        
        self.x = np.asarray(x)
        return x
    
    
    def tangential_forces(self, dt=1000):        
        F = list()
        r = list()         
        t=np.linspace(0, 2*pi, dt, endpoint=True)
        
        # compute the tangential force on one current loading layer ...
        for i, cl in enumerate(self.current_loadings):
            r.append(cl.r)
            j = cl.idx * 2
            Br = 0.
            
            # ... by multiplying the radial flux created by 
            # all other current loadings at its radius ...
            for sm in self.submodels[:i] + self.submodels[i+1:]:
                Br += sm.layers[cl.idx].Br(self.p, r[i], t,
                                           a_j = sm.x[j],
                                           b_j = sm.x[j+1])
            
            # ... with it's own current loading ...
            Kt = cl.Kt(self.p, t)
            
            # ... and integrating over it using np.trapz().
            f = np.trapz(Kt * Br * r[i], t)
            cl.tangential_force = f
            F.append(f)
        return F, r
        
    
    def total_torque(self, dt=1000):
        self.tangential_forces(dt)        
        pos_torque = []
        neg_torque = []
        
        # store the acting forces in either direction multiplied
        # with the respective radius in 2 lists
        for cl in self.current_loadings:
            if cl.tangential_force >= 0.:
                pos_torque.append(cl.tangential_force * cl.r)
            else:
                neg_torque.append(cl.tangential_force * cl.r)
        
        # check if the torque in both directions is equal
        M = sum(pos_torque)
        if np.allclose(M, -sum(neg_torque)):
            return M
        else:
            raise Exception("Pos. and neg. torque are NOT np.allclose().")


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
        # data = (X, Y, Az, R, T)
        return X, Y, Az, R, T
    
    
    def get_B_data(self, r, t):
        R_tuple, T_tuple = tuple(), tuple()
        Br_tuple, Bt_tuple = tuple(), tuple()
        
        # computes the flux densities for the i-th layer per iteration
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
                      
            # start computing the superposition of the flux densities
            this_Br = np.zeros(this_R.shape)
            this_Bt = np.zeros(this_R.shape)
            
            for subm in self.submodels:
                
                # this does not add a new layer to the model but is used to compute
                # the field for the environment, otherwise a_n & b_n would be unused
                plot_layers = subm.layers + [Environment()]
                
                # sums up the radial flux density for all current loadings
                this_Br += plot_layers[i].Br(self.p, this_R, this_T, 
                                             a_j = subm.x[j], 
                                             b_j = subm.x[j+1])
                # sums up the tangential flux density for all current loadings
                this_Bt += plot_layers[i].Bt(self.p, this_R, this_T, 
                                             a_j = subm.x[j], 
                                             b_j = subm.x[j+1])
            
            # store the result for the i-th layer
            R_tuple += (this_R, )
            T_tuple += (this_T, )
            Br_tuple += (this_Br, )
            Bt_tuple += (this_Bt, )

        R, T = np.hstack(R_tuple), np.hstack(T_tuple)
        Br, Bt = np.hstack(Br_tuple), np.hstack(Bt_tuple)
        
        X, Y = rt_to_xy(R, T)
        U, V = BrBt_to_UV(Br, Bt, T)
        # data = (X, Y, U, V, R, T, Br, Bt)
        return X, Y, U, V, R, T, Br, Bt


    def get_H_data(self, r, t):
        R_tuple, T_tuple = tuple(), tuple()
        Hr_tuple, Ht_tuple = tuple(), tuple()
    
        # computes the field strength for the i-th layer per iteration
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
                      
            # start computing the superposition of the flux densities
            this_Hr = np.zeros(this_R.shape)
            this_Ht = np.zeros(this_R.shape)
            
            for subm in self.submodels:
                
                # this does not add a new layer to the model but is used to compute
                # the field for the environment, otherwise a_n & b_n would be unused
                plot_layers = subm.layers + [Environment()]
                
                # sums up the radial flux density for all current loadings
                this_Hr += plot_layers[i].Hr(self.p, this_R, this_T, 
                                             a_j = subm.x[j], 
                                             b_j = subm.x[j+1])
                # sums up the tangential flux density for all current loadings
                this_Ht += plot_layers[i].Ht(self.p, this_R, this_T, 
                                             a_j = subm.x[j], 
                                             b_j = subm.x[j+1])
    
    
    
            # store the result for the i-th layer
            R_tuple += (this_R, )
            T_tuple += (this_T, )
            Hr_tuple += (this_Hr, )
            Ht_tuple += (this_Ht, )
    
        R, T = np.hstack(R_tuple), np.hstack(T_tuple)
        Hr, Ht = np.hstack(Hr_tuple), np.hstack(Ht_tuple)
        
        X, Y = rt_to_xy(R, T)
        U, V = BrBt_to_UV(Hr, Ht, T)
        # data = (X, Y, U, V, R, T, Hr, Ht)
        return X, Y, U, V, R, T, Hr, Ht
    
