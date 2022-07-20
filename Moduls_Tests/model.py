# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:12:38 2022

@author: svens
"""
import numpy as np
import Moduls_Tests.somemath as sm

from .layer import Layer, CurrentLoading

mu_0 = 4 * np.pi * 10**(-7) # [N / A**2]

class Model():
    """Header
    
    Description
    
    Args:
        
        
    Attributes:
        
        
    """
    
    
    def __init__(self, p: int = 1):
        
        if p is None:
            p = 1
        
        self.p = p
        self.n_layers = 0
        self.layers = dict()
        self.built = False
        self.solution = None
         
    
    def add_layer(self, layer: Layer):
        
        layer.add_index(self.n_layers)
        layer.sync_pole_pairs(self.p)
        self.layers[layer.index] = layer
        
        self.n_layers += 1
                
        
    def build(self):
        
        if self.n_layers == 0:
            raise Exception("Model contains no layers. Use add_layer() before trying to build the model.")
        
        for i in range(self.n_layers-1):
            self.layers[i].sync_layer_info(self.layers[i+1].mu_inv)
        
        n_bounds = sum([len(vals.bounds) for vals in self.layers.values()])
        n_vars = n_bounds + 2
        
        self.M = np.zeros((n_vars, n_vars))
        self.y = np.zeros(n_vars)
        
        self.M[ 0,  1] = 1
        self.M[-1, -2] = 1
        
        
        i = 0
        for layer in self.layers.values():
            
            if isinstance(layer,CurrentLoading):
                self.y = layer.apply_solution(self.y, i)
            
            self.M = layer.apply_boundaries(self.M, i)
            i += len(layer.bounds)
            
        self.built = True
            
            
    def solve(self, allclose_check: bool= False):
        if self.built:    
            x = np.linalg.solve(a = self.M, b = self.y)
            self.solution = x
            if allclose_check:
                return x, np.allclose(np.dot(self.M, x), self.y)    
            else:
                return x
        else:
            raise Exception("Model has not been built yet. Use build() before trying to solve the model.")
                    
    
    def get_radii_data(self):
        radii = []
        for lay in self.layers.values():
            radii.append(lay.r)
        return radii
    
    
    
    def get_B_plot(self, theta):
        
        r_i, r_o = 0.01, 3
        r = np.linspace(r_i, r_o, 30)
        t = np.linspace(0, 2*np.pi, 40)
        # R, T = np.meshgrid(r, t)
        R_tuple, T_tuple = tuple(), tuple()
        Br_tuple = tuple()
        Bt_tuple = tuple()
        
        for lay in self.layers.values():
            if lay.index == 0:
                r_i = 0.01
            else:
                r_i = self.layers[lay.index -1].r
            
            r_o = lay.r
            
            lower = r >= r_i
            upper = r <   r_o
            idxs = np.argwhere(lower & upper).flatten()
            r_lay = r[idxs]
            
            R_lay, T_lay = np.meshgrid(r_lay, t)
            
            Br_lay = sm.Br_no_k(self.p,
                                R_lay, 
                                T_lay,
                                aj = self.solution[lay.index * 2],
                                bj = self.solution[lay.index * 2 + 1]
                                )
            
            Bt_lay = sm.B0_no_k(self.p,
                                R_lay, 
                                T_lay,
                                aj = self.solution[lay.index * 2],
                                bj = self.solution[lay.index * 2 + 1]
                                )
            
            R_tuple += (R_lay, )
            T_tuple += (T_lay, )
            Br_tuple += (Br_lay, )
            Bt_tuple += (Bt_lay, )
            
        R = np.hstack(R_tuple)
        T = np.hstack(T_tuple)        
        Br = np.hstack(Br_tuple)
        Bt = np.hstack(Bt_tuple)
        X, Y = sm.r0_to_xy(R, T)
        U, V = sm.BrB0_to_UV(Br, Bt, T)
            
        return X, Y, U, V
    
    
    def get_A_plot(self):
        r_i, r_o = 0.01, 3
        r = np.linspace(r_i, r_o, 30)
        t = np.linspace(0, 2*np.pi, 40)
        
        R_tuple, T_tuple = tuple(), tuple()
        Az_tuple = tuple()
        
        for lay in self.layers.values():
            if lay.index == 0:
                r_i = 0.01
            else:
                r_i = self.layers[lay.index -1].r
            
            r_o = lay.r
            
            lower = r >= r_i
            upper = r <   r_o
            idxs = np.argwhere(lower & upper).flatten()
            r_lay = r[idxs]
            
            R_lay, T_lay = np.meshgrid(r_lay, t)
            
            Az_lay = sm.Az_no_k(self.p,
                                R_lay, 
                                T_lay,
                                aj = self.solution[lay.index * 2],
                                bj = self.solution[lay.index * 2 + 1]
                                )
            
            R_tuple += (R_lay, )
            T_tuple += (T_lay, )
            Az_tuple += (Az_lay, )

        R = np.hstack(R_tuple)
        T = np.hstack(T_tuple)        
        Az = np.hstack(Az_tuple)
        X, Y = sm.r0_to_xy(R, T)
        return X, Y, Az

        
        
        
        
        
    
    
    
    