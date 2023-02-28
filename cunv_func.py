#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 10:48:43 2022

@author: cristina
"""

import numpy as np

def convolution(V, Delta_tip, Delta_subs, Dynes, T_tip, T_sub, N_subs, N_tip, N, tip, sigma):    

    import numpy as np
    ###common functions
    exp = np.exp
    sqrt = np.sqrt

    ###parameters
    Delta_tip = Delta_tip/27211.6 #SC gap
    Delta_subs = Delta_subs/27211.6 #SC gap
    Dynes = Dynes/27211.6
    T_tip = T_tip*25.86/300/27211.6
    T_sub = T_sub*25.86/300/27211.6    

    
    ####define Fermi functions
    def fermi(E,T):
        y=1.0/(1+exp(E/T))
        return(y)
    
    def dfermi(E,T):
        y = (-1.0/T)*exp(E/T)/(1.0+exp(E/T))**2
        return(y)
    
    ####SC density of states
    def rho0(omega, Dynes, Delta, N0):
    
        y = N0*np.abs(omega*np.real(1/sqrt(omega**2 + 2*1j*omega*Dynes - (Dynes + Delta)**2)))        
        return(y)

    def rhod(omega, Dynes, Delta, N0):    
        y1 = N0*np.sign(omega)*np.real(1 / sqrt(omega ** 2 + 1j*omega*Dynes - (Delta)**2))
        y2 = - N0*np.sign(omega)*np.real(omega**2/(omega**2 + 1j*omega*Dynes - Delta**2)**(3.0/2.0))    
        return(y1 + y2)
    
    ####constant density of states    
    def rho_cte(N0):    
        y = N0        
        return(y)  
        
    ######Gaussian noise    
    def gauss_noise(dIdV, omega, sigma, N):
        u = np.exp(-1/2*(omega**2/sigma**2))/(sigma*np.sqrt(np.pi))        
        dIdV_new = np.convolve(u, dIdV)/N
        dIdV_3 = dIdV_new[int((N-1)/2) : int((N-1)/2 + N)] 
        
        return(dIdV_3)
    

    #####conductance calculation
    def dIdV_1(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip):
        Vr =  -V 
        u = rho0(V,Dynes,Delta_subs, N_subs) 
        v = rhod(Vr,Dynes,Delta_tip, N_tip)*fermi(Vr,T_tip)
        dIdV = np.convolve(u,v)
        return(dIdV)
    
    def dIdV_2(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip):
        Vr =  -V 
        u = (rho0(V,Dynes,Delta_subs, N_subs))*fermi(V,T_sub)
        v = rhod(Vr,Dynes,Delta_tip, N_tip)
        dIdV = np.convolve(u,v)
        return(dIdV)
    
    def dIdV_3(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip):
        Vr =  -V
        u = rho0(V,Dynes,Delta_subs, N_subs)  
        
        if (tip == 'SC'):
            v = rho0(Vr,Dynes, Delta_tip, N_tip)*dfermi(Vr,T_tip)            
        if (tip == 'Metal'):
            v = rho_cte(N_tip)
            
        dIdV = np.convolve(u,v)
        return (dIdV)    

             
    if (tip == 'SC'):
        dIdV = (- dIdV_1(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip) + \
             dIdV_2(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip) - \
             dIdV_3(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip))/N
                
        dIdV_2 = dIdV[int((N-1)/2) : int((N-1)/2 + N)] 
                
    elif (tip == 'Metal'):  
        
        dIdV = dIdV_3(V, Dynes, Delta_tip, Delta_subs, N_subs, N_tip)/N
        dIdV_2 = dIdV
        
    ###add gaussian noise    
    dIdV_3 =  gauss_noise(dIdV_2, V, sigma, N)

    
    
    return(dIdV_2, dIdV_3)  

    
    
    
    
    
    
    
    
    