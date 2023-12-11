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
# temperature (in degrees C) #PROVEN #Clausisu-Clapeyron
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




#
#------------------------------------------------------------------------------
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
        airtosea:
            temperature difference from sea to air
            
        humid:
            relative humidity
        
    function:
        for a given difference in air to sea temperature we have to consider
        and include values for mass mixing ratio of the air at temperatures 
        below the give low temperature in the range.
        This function creates a list with length of values required which
        depends on the air-sea temp diff value.
    
    outputs:
        list of mass mixing ratio values below the low temperature in the
        range. For use in the bowen ratio equation.

'''


#
#-----------------------------------------------------------------------------
def bowenratio(qstarvals=qstarvals, humid=0.70, airtosea=2, equilibrium=False):
    
    #Calculate mass mixing ratio q, given a humidity and RH = q/q*
    qvals = [] 
    # Solving for q values at q(temp - airtoseatempdiff) for use in bowen ratio
    for i in range(len(qstarvals)):          
        vals = qstarvals[i-airtosea] * humid
            
        qvals.append(vals)
        
    #Adding in boundary q values
    #q value for i[0-airtosea] boundary
    qvals[0:airtosea] = boundary(airtosea,humid)
    #COMMENT OUT THE LINE ABOVE AND RUN TO SEE WHY IT'S IMPORTANT!
    
    
#Calculate Bowen Ratio given air-sea temp diff
    bowenratiovals=[]
    for i in range(len(qvals)):
        if equilibrium==True:
            
            #Calculate Equilibrium Bowen Ratio
            #If calculating Be, Be^-1 = (L/Cp) * dq*/dT
            boweninv = (L/Cp) * qstarvals[i]
            #units issue below so multiply by 10
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
                (our calcqstarval function)
    
    function:
        Bowen Ratio equation:
        Bo = (Cp * airseatempdiff) / (L* (qsurface - qair) )
        where qsurface = q*(T)
        and   qair     = q(T - airseatempdiff)
        
        
        First For loop calculates qair values:
            Because mass mixing ratio depends on temp, q-air value used in 
            bowen ratio has to take into account the air to sea temperature 
            difference.
        Here we calculate q-air(T - airseatempdiff) values from: 
            q-air(T - airseatempdiff) = qsurface(T - airseatempdiff) * humidity

        As described in the -boundary- docstring we must work to include values
        of qair below the temp range. Calling the boundary function here makes
        sure our bowen ratio line is complete.
        
        
        To calculate the bowen ratio at each temperature as shown in the 
        equation above 
        the second For loop runs through every temperature, 
        inputting q*(qsurface) and q(qair) values at each temp as well as the
        inputted airtosea difference.
        
    outputs:
        bowenratiovals:
            Values of calculated Bowen Ratio


'''



humiddiffarray=[]
humidrange = np.arange(0,1,0.1)
for i in humidrange:
      
    brvals = bowenratio(humid=i)
    avgbr = np.mean(brvals)
    humiddiffarray.append(avgbr)
    
    
plt.semilogy(humidrange, humiddiffarray)
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


