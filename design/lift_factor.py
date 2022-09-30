# -*- coding: utf-8 -*-
import numpy as np
from scipy.constants import pi

wt = 12
bt = 0.22
bSC = 3e-3
A = wt * bt
ASC = wt * bSC

Ic_0 = 300

linear_fun = lambda a, x: a[0] * x + a[1]
exponential_decay = lambda a, x: a[0] * np.exp(-(x/a[1])**a[2])
parabolic_increase = lambda a, x: a[0] * x**a[1] + a[2]
angle_dependency = lambda a, x: a[0] / (1 + a[1] * np.abs(np.sin(x - (a[3]/360) * 2*pi))**a[2]) + a[4]

Ic_fit_final  = np.array([8.600606578930105, 77.929267595050604, 0.614684024969411])
a_fit_final   = np.array([2.425185286960588, 0.003257562217890, 0.000341639734764])*1e3
b_fit_final   = np.array([-27.358801161089112, 4.373218208831309, 0.314268659635138])
exp_fit_final = 0.515800000000000

a1_fit_final = np.array([-0.088088589245889, 4.682645887291001])
a2_fit_final = np.array([ 0.000000000004849, 8.194684157196530, 9.990790629451066])
a3_fit_final = np.array([ 0.051247695994478, 0.170662898258885])
a4_fit_final = 1.216070888376809*1e2
a5_fit_final = np.array([-0.034902355980912, 2.623684547437763])

def get_L_TPL2100(T, B, theta):
    
    T_final = T
    B_final = B
    theta_final = (theta/180) * pi
    theta_critical = (30/180) * pi
    
    params_final = [exponential_decay(Ic_fit_final, B_final)/A, 
                    exponential_decay(a_fit_final, B_final),
                    exponential_decay(b_fit_final, B_final),
                    exp_fit_final]
      
    params_angle_final = [linear_fun(a1_fit_final,T_final),
                          parabolic_increase(a2_fit_final,T_final),
                          linear_fun(a3_fit_final,T_final),
                          a4_fit_final,
                          linear_fun(a5_fit_final,T_final)]
    
    num = angle_dependency(params_angle_final, theta_final)
    den = angle_dependency(params_angle_final, theta_critical)
    L_angle_normalized =  num / den
    
    f = params_final[0] * np.exp(- (T_final / (params_final[1] + params_final[2] * T_final))**params_final[3])
    
    L = 10**f * (A/Ic_0) * L_angle_normalized
    
    return L