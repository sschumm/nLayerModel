# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:38:57 2022

@author: svens
"""
import numpy as np


def r0_to_x(r, theta):
    x = r * np.cos(theta)
    return x


def r0_to_y(r, theta):
    y = r * np.sin(theta)
    return y


def r0_to_xy(r, theta):
    x = r0_to_x(r, theta)
    y = r0_to_y(r, theta)    
    return x, y



