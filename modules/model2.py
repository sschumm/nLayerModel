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


class Model():
    def __init(self, p: int):
        
        self.p = p
        
        
    
        
        