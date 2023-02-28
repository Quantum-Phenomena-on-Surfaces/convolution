#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 10:42:33 2022

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
import cunv_func as cv

plt.rcParams.update({'font.size': 18})
plt.rcParams['figure.figsize'] = [20, 6]
#from scipy.signal import find_peaks


'''Parameters'''
tip = 'SC'####SC of Metal
Delta_tip = 1.0
Delta_subs = 2.0
Dynes = 0.05
T_tip = 4.0
T_sub = 4.0 
N_subs = 4.5E-4
N_tip = 4.5E-4
V_i = 8.0

###energy vector
V_init = -V_i/27211.6
V_final = V_i/27211.6
N = 2001
V = np.linspace(V_init,V_final,N)
bias = V*27211.6

#####gaussian noise
sigma = 2.7E-6

'''load data'''
#data = np.loadtxt('SC_pb_tip', skiprows = 1)

########convoluted dI/dV
(dIdV, dIdV_gauss) = cv.convolution(V, Delta_tip, Delta_subs, Dynes, T_tip, T_sub, N_subs, N_tip, N, tip, sigma)
factor = np.average(dIdV)/np.average(dIdV_gauss)
factor_data = np.average(dIdV)/np.average(data[:,1])


if (tip == 'SC'):
    plt.plot(bias, dIdV, color='C0', label = 'dI/dV')
    plt.plot(bias, dIdV_gauss*factor, color='C1', label = 'Gaussian noise')
    
    #plt.plot(data[:,0]*1.0E+3, data[:,1]*factor_data, '.', color = 'C2', label = 'data')
    
    plt.xlim([-5, 5])
    plt.xlabel('Bias (mV)')
    plt.ylabel('dI/dV (a.u.)')
    plt.legend()
    
elif (tip == 'Metal'):
    plt.plot(bias, dIdV, color='C0', label = 'Metal tip')
    plt.plot(bias, dIdV_gauss*factor, color='C1', label = 'Gaussian noise')
    plt.xlim([-5, 5])
    plt.xlabel('Bias (mV)')
    plt.ylabel('dI/dV (a.u.)')    
    plt.legend()
    
    
    
    
    
    