# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:12:38 2022

@author: svens
"""
import numpy as np
import Moduls_Tests.somemath as sm

from .somemath import Az_no_k, Br_no_k, B0_no_k
#from .somemath import *
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
        
        
    
    def test(self, r, theta):
        if not isinstance(r, np.ndarray):
            r = np.array([float(r)])
            
        Az = []
        for radius in r:
            Az.append(Az_no_k(self.p, 
                              radius, 
                              theta,
                              aj = self.solution[0], 
                              bj = self.solution[1]))
            
        return np.asarray(Az)
    
    
    def get_Br_plot(self, theta):
        
        # Br = []

        spacing = self.layers[list(self.layers)[-1]].r / 20.
                
        for lay in self.layers.values():
            if lay.index == 0:
                r_i = 0.01
            else:
                r_i = self.layers[lay.index -1].r
            radius = np.arange(r_i, lay.r, spacing)
            
            X, Y = [], []
            U, V = [], []
            Br, B0 = [], []
            for r in radius:
                X.append(sm.r0_to_x(r, theta))
                Y.append(sm.r0_to_y(r, theta))
                # Br.append(f"Radius = {r}, Br = ")
                Br.append(Br_no_k(self.p,
                                  r, 
                                  theta,
                                  aj = self.solution[lay.index * 2],
                                  bj = self.solution[lay.index * 2 + 1]
                                  )
                          )
                B0.append(B0_no_k(self.p,
                                  r, 
                                  theta,
                                  aj = self.solution[lay.index * 2],
                                  bj = self.solution[lay.index * 2 + 1]
                                  )
                          )
                
            lay.X = np.asarray(X)
            lay.Y = np.asarray(Y)
            U, V = sm.BrB0_to_UV(Br, B0, theta)
            lay.U = U
            lay.V = V
            lay.Br = Br
            lay.B0 = B0
        return"Job finished"
            # np.set_printoptions(suppress=True, linewidth=200, precision=4)
            # print(Br, "\n")
        



# =============================================================================
#         for lay in self.layers.values():
#             if lay.index == 0:
#                 r_i = 0
#             else:
#                 r_i = self.layers[lay.index -1].r
#             radius = np.linspace(r_i, lay.r, 10)
#             print(lay.index, "  ", radius)
# =============================================================================
        
# =============================================================================
#         if not isinstance(r, np.ndarray):
#             r = np.array([float(r)])
#             
#         Az = []
#         for radius in r:
#             Az.append(Az_no_k(self.p, 
#                               radius, 
#                               theta,
#                               aj = self.solution[0], 
#                               bj = self.solution[1]))
#             
#         return np.asarray(Az)
# =============================================================================
         

        
    
        
        
        
        
        
    
    
    
    