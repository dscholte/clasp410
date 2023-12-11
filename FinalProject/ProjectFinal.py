# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 16:38:03 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt

L = 2264000 #j/kg #latent heat of vaporization
Ps = 101325 #Pa #surface pressure
Cp = 1004 #J/kg/K #specific heat
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

'''
-e_air-
function given to use for the clasius-clapeyron relation

    inputs:
        t_air:    
            takes in a temperature in degrees celsius
        
    outputs:
        vappress:
            saturation vapor pressure at a certain temperature

'''




# Run loop to solve for saturated mixing ratio in range in degree celsius 
# (set numbers) #PROVEN
def calcqstarval(templow,temphigh):
    
    #Saturation Vapor Pressure for given range
    sat_vparray = []
    for i in range(templow,temphigh+1,1):
        sat_vp = e_air(i)
        sat_vparray.append(sat_vp)

    #Saturated mixing ratio for given range
    #q* = 0.622 * e*/p
    qstarvals = []
    for i in range(len(sat_vparray)):
        qstarvalshere = ( sat_vparray[i] / Ps ) * 0.622
        qstarvals.append(qstarvalshere)
    
    return qstarvals, qstarvalshere, templow, temphigh+1

'''
-calcqstarval-

inputs:
    templow and temphigh:
        temperatures inputted in degrees celsius
        
outputs:
    qstarvals:
        this creates the list of saturated mixing ratio values at each temp
        outputted with units of kg/kg
    
    qstarvalshere: 
        this second output for the boundary condition used next,
        needed to remove the value from appended list to use as float.
        I don't know how else to do this but it seems to work fine.
        
    templow and temphigh+1
        outputted for use in plots
'''




#Running calcqstarval for -40 to 40
#Sat. Mix Ratio (q*) from (-40) to (40) degree celsius
qstarvals, qstarvalshere, templow, temphigh = calcqstarval(-40,40)
#Create temp range from inputted range for quick use in plots
temp = np.arange(templow,temphigh,1)






def boundary(airtosea, humid):
    
    boundaryqvallist=[]
    boundaryvallist=[]
    
    #calc q* vals at boundary depending on air to sea diff subtraction
    for i in range(1, airtosea+1):
    #for i in range(airtosea+1, 1):
        
        notused, boundaryval , notused2, notused3 = calcqstarval(templow-(i), \
                                                            templow-(i))
        boundaryvallist.append(boundaryval)
    
    
    #multiply values by humidity to get q val at boundary
    for i in range(len(boundaryvallist)):
        boundaryqval = boundaryvallist[i] * humid
        boundaryqvallist.append(boundaryqval)
    
    
    boundaryqvallist.reverse()
    
    return boundaryqvallist

'''
-boundary-
    inputs:
        
        
    function:
        
    
    outputs:
        

'''


#Solution Part 2 / 3 
#-----------------------------------------------------------------------------
def bowenratio(qstarvals=qstarvals, humid=0.70, airtosea=2, equilibrium=False):
    
    #Calculate mass mixing ratio q, given a humidity and RH = q/q*
    
    #Consider boundary of temp range for qvals
    #Use to solve for q(T-airtoseatempdiff)
    
    #q* value for i[0-airtosea] boundary
    specialqstarval2,specval2 , no, yes = calcqstarval(templow-(airtosea), \
                                                       templow-(airtosea))
    #q* value for i[1-airtosea] boundary
    specialqstarval1,specval1 , no, yes = calcqstarval((templow+1)-(airtosea), \
                                                       (templow+1)-(airtosea))
    
    
    qvals = []
    # RH = q/q* 
    # Solving for q(temp - airtoseatempdiff) for use in bowen ratio
    for i in range(len(qstarvals)):          
        vals = qstarvals[i-airtosea] * humid
            
        qvals.append(vals)
        
    #q value for i[0-airtosea] boundary
    qvals[0] = specval2 * humid
    #q value for i[1-airtosea] boundary
    qvals[1] = specval1 * humid
    
    
#Calculate Bowen Ratio given air-sea temp diff
    bowenratiovals=[]
    for i in range(len(qvals)):
        if equilibrium==True:
            
            #Calculate Equilibrium Bowen Ratio
            #If calculating Be, Be^-1 = (L/Cp) * dq*/dT
            boweninv = (L/Cp) * qstarvals[i]
            bowen = (1/boweninv)*10
            bowenratiovals.append(bowen)
        
        else:
            #Otherwise calculating Bo values for temprange
            # Bo = (Cp*air-sea-temp-diff) / (L* (qsurface-qair) )
            bowen = (Cp*airtosea) / (L*( qstarvals[i] - qvals[i]))
            bowenratiovals.append(bowen)
    
    
    return bowenratiovals

'''
-bowenratio-

    inputs:
        qstarvals:
            input range of saturated mixing ratio values to use in bowen
            ratio equation.
            
        humid:
            constant relative humidity
        
        airtosea:
            temperature difference from the sea to the air
            
        equilibrium:
            if solving for Be, equilibrium bowen ratio, equilibrium=True
            when True:
                calculated bowen ratios use different equation:
                Be^-1 = (L/Cp) * dq*/dT  
                where dq*/dT is = clasius-clapeyron relation used for q* vals
    
    function:
        Because mass mixing ratio depends on temp, q (q-air) value used in bowen ratio
        has to take into account the air to sea temperature difference. At the
        lower boundary it's necessary to add '
        
        
        
        
        Bo = (Cp*air-sea-temp-diff) / (L* (qsurface-qair) )
        
        using inputted q* values and boundary q* values, a list of q values
        used in the bowen ratio equation is created using RH = q/q* :
            
            Bo = (Cp*air-sea-temp-diff) / (L* (qsurface-qair) )
            where qsurface = q*
            and   qair     = q
        
        
        
    outputs:
        bowenratiovals:
            Values of calculated Bowen Ratio


'''





brvals = bowenratio()
plt.semilogy(temp, brvals)
plt.show()

#Recreating Figure given
eqbrvals = bowenratio(equilibrium=True)
# #g/kg units
qstarvals = np.array(qstarvals)*1000 
plt.semilogy(temp,qstarvals,'b',label='q* (g/kg)')
plt.semilogy(temp,eqbrvals,'r',label="$B_{e}$", linestyle='--')
plt.xlim(templow,temphigh)
plt.title('Recreated q* and Equilibrium Bowen Ratio')
plt.xlabel('Temperature (degC)')
plt.legend()
plt.show()


