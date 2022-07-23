# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 17:24:t2 2t22

@author: svens
"""

import numpy as np



xy_to_r = lambda x, y: np.sqrt(x**2, y**2)
xy_to_t = lambda x, y: np.arctan2(y, x)


def xy_to_rt(x, y):
    return xy_to_r(x,y), xy_to_t(x,y)

rt_to_x = lambda r, t: r * np.cos(t)
rt_to_y = lambda r, t: r * np.sin(t)


def rt_to_xy(r, t):
    return rt_to_x(r, t), rt_to_y(r, t) 


BrBt_to_U = lambda Br, Bt, t: Br * np.cos(t) - Bt * np.sin(t)
BrBt_to_V = lambda Br, Bt, t: Br * np.sin(t) + Bt * np.cos(t)


def BrBt_to_UV(Br, Bt, t):
    return BrBt_to_U(Br, Bt, t), BrBt_to_V(Br, Bt, t)