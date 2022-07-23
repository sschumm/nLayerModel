# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 16:11:09 2022

@author: svens
"""



def c_Br(r, p):
    return [-r**(p - 1), 
            -r**-(p + 1),
             r**(p - 1),
             r**-(p + 1)]


def c_Ht(r, p, mu_i_inv, mu_o_inv):
    return [ mu_i_inv * p * r**(p - 1),
            -mu_i_inv * p * r**-(p+ 1),
            -mu_o_inv * p * r**(p - 1),
             mu_o_inv * p * r**-(p+ 1)]