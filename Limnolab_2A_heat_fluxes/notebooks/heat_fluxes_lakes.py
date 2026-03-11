#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:14:47 2024

@author: Damien Bouffard
@email: damien.bouffard@eawag.ch
"""

import numpy as np

# Function in Python for short wave absorption

def SW(G, C, Adir):
    Adiff = 0.066
    Adir = 0.1 # should become time dependant
    
    Fdir = (1-C)/((1-C)+ 0.5*C)
    Fdiff = (0.5*C)/((1-C)+0.5*C)
    

    
    sw = ( G * Fdir * (1 - Adir) + G * Fdiff * (1 - Adiff) )
    return sw

# Function in Python for Absorption of atmospheric longwave radiation

def LW_in(Ta,Ea):
    
    # parameters
    A_L = 0.03
    sigma = 5.67 * 10**-8
    
    # Absorption of atmospheric longwave radiation
    lwi =  (1-A_L)*Ea*sigma*(Ta+273.25)**4
    
    return lwi

#def vapour_pressure(Ta, rh):
#
#    # water vapour pressure (and water vapour saturation )
#    es = 6.11*np.exp( (17.27*Ta)/(237.3+Ta) )
#    ea = rh*es/100.
    
#    return ea
    

def emissivity_(Ta, ea, C):
    #parameters
    a = 1.09
    
    #Emissivity
    Ea = a * (1 + 0.17 * C**2) * 1.24  * (ea / (Ta+273.25))**(1/7)

    return Ea

def emissivity_old(Ta, rh, ea, C):
    #parameters
    a = 1.09
    
    #Emissivity
    Ea = a * (1 + 0.17 * C**2) * 1.24  * (ea / (Ta+273.25))**(1/7)

    return Ea

# Function in Python for Longwave emission from water surfaces

def LW_out(Tw):
    
    # parameters
    E_w = 0.972
    sigma = 5.67 * 10**-8
    
    # Longwave emission from water surfaces
    lwo = - E_w * sigma * (Tw+273.25)**4
    
    return lwo

# Function in Python for Evaporation and condensation

def He(Ta, Tw, WS, ea):
    
    es = 6.112 * np.exp(17.62 * Tw / (243.12 + Tw))
    f = 4.8 + 1.98 * WS + 0.28 * (Tw - Ta)
    
    # Evaporation and condensation
    he = - f * (es - ea)  
    
    return he

def Hc(Ta, Tw, WS):
    
    cp = 1005
    p = 1015
    Lv = 2.47 * 10**6
    gamma = cp * p / (0.622 * Lv)

    f = 4.8 + 1.98 * WS + 0.28 * (Tw - Ta)
    
    # Evaporation and condensation
    hc = - gamma * f * (Tw - Ta) 
    
    return hc

def Hf(Q,Tr,Tl, A0):
    cp= 4180
    rho=1000
    
    hf = -cp * rho * Q * (Tl - Tr) / A0
    
    return hf