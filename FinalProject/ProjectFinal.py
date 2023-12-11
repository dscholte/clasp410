# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 16:38:03 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt

L = 2264000 #j/kg
Ps = 101325 #Pa 
Cp = 1004 #J/kg/K
Rv = 461 #J/kg/K

#Solution Part 1 
#-----------------------------------------------------------------------------

# Function returns the saturation vapor pressure e* in [Pa], given air
# temperature (in degrees C) #PROVEN 
def e_air(t_air):
    a =6.107799961
    b =4.436518521E-1
    c =1.428945805E-2
    d =2.650648471E-4
    e =3.031240396E-6
    f =2.034080948E-8
    g =6.136820929E-11
    vappress = 100*(a + b*t_air + c*t_air**2 + d*t_air**3 + e*t_air**4 + \
                    f*t_air**5 + g*t_air**6);
    
    return vappress

# Run loop to solve for saturated mixing ratio in  range in degree celsius 
# (set numbers) #PROVEN
def calcqstarval(templow,temphigh):
    
    
    sat_vparray = []
    for i in range(templow,temphigh+1,1):
        sat_vp = e_air(i)
        sat_vparray.append(sat_vp)

    #q* = 0.622 * e*/p
    qstarvals = []
    for i in range(len(sat_vparray)):
        qstarvalshere = ( sat_vparray[i] / Ps ) * 0.622
        qstarvals.append(qstarvalshere)
    
    return qstarvals, qstarvalshere, templow, temphigh+1

#Sat. Mix Ratio (q*) from (-40) to (40) degree celsius
qstarvals, qstarvalshere, templow, temphigh = calcqstarval(-40,40)
temp = np.arange(templow,temphigh,1)
#Solution Part 2 / 3 
#-----------------------------------------------------------------------------
#Given a constant relative humidity = 0.70
#Given constant air-sea temperature difference = 2
#Calculate q values
#RH = q/q*
def bowenratio(qstarvals=qstarvals, humid=0.70, airtosea=2, equilibrium=False):
    
    #Calculate mass mixing ratio q, given a humidity and RH = q/q*
    if equilibrium==True:
        humid=1
#Consider boundary of temp range for qvals
    specialqstarval2,specval2 , no, yes = calcqstarval(templow-(airtosea),templow-(airtosea))
    specialqstarval1,specval1 , no, yes = calcqstarval(templow-(airtosea-1),templow-(airtosea-1))
    
    qvals = []
    for i in range(len(qstarvals)):
        
        if i==0:
            vals = specval2 * humid
            
        if i==1:        
            vals = specval1 * humid
            
        else:
            vals = qstarvals[i] * humid
            
        
        qvals.append(vals)    


#Calculate Bowen Ratio given air-sea temp diff
    bowenratiovals=[]
    for i in range(len(qvals)):
        if equilibrium==True:
            
            #If calculating Be, Be^-1 = (L/Cp) * dq*/dT
            boweninv = (L/Cp) * qstarvals[i]
            bowen = (1/boweninv)*10
            bowenratiovals.append(bowen)
        
        else:
            #Otherwise calculating Bo values for temprange
            # Bo = (Cp*air-sea-temp-diff) / (L* )
            bowen = (Cp*airtosea) / (L*( qstarvals[i] - qvals[i]))
            bowenratiovals.append(bowen)
    
    
    return bowenratiovals


#Recreating Figure given
brvals = bowenratio(equilibrium=True)
# #g/kg units
qstarvals = np.array(qstarvals)*1000 
plt.semilogy(temp,qstarvals,'b',label='q* (g/kg)')
plt.semilogy(temp,brvals,'r',label="$B_{e}$", linestyle='--')
plt.xlim(templow,temphigh)
plt.legend()



