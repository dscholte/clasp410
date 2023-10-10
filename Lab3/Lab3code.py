# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 10:21:54 2023

@author: dscho, manishrv
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

plt.style.use("fivethirtyeight")

@dataclass
class modeldata:
    '''
    Holds the information for all variables in the model
    
        temperatures
            outputted temperature of layer
        layers
            List of each layer 
        nlayers
            integer of number of layers in model
        solar
            solar irradiance value (w/m^2)
        emissivity
            layer emission value between 0 and 1
        albedo
            reflectivity value between 0 and 1
        nuclear_winter
            
    '''
    
    temperatures: list[float]
    layers:list[float]
    nlayers:int
    solar:float
    emissivity:float
    albedo:float
    nuclear_winter: bool
    
def run_atmmodel(nlayers, solar, emissivity, albedo, nuclear_winter=False):
    
    
    
    nlayers = nlayers
    e_atm = emissivity
    e_ground = 1
    s_0 = solar
    if nuclear_winter==True:
        albedo=0
    s = s_0 * (1/4) * (1 - albedo)
    sigma = 5.67*(10**-8)
    
    temperatures = np.zeros(nlayers+1)
    
    
    A = np.zeros([nlayers+1, nlayers+1])
    b = np.zeros(nlayers+1)
    
    for i in range(nlayers+1):
        for j in range(nlayers+1):
            value = -9999
            if i == 0:
                e_layer = e_ground
            else:
                e_layer = e_atm
                
            if i==0 and j==0:
                value = -1
            elif i == j:
                value = -2
            else:
                value = e_layer* (1-e_atm)**(abs(i-j)-1)
            A[i,j] = value
         
            
    #fill b
    if not nuclear_winter:    
        b[0] = -s
    else:
        b[-1] = -s
        
        
        
    # Invert matrix:
    Ainv = np.linalg.inv(A)
    # Get solution of Fluxes with matrix multiplication:
    fluxes = np.matmul(Ainv, b) # Note our use of matrix multiplication!
    
    
    #Convert fluxes to temperature
    
    #Atmospheres
    temperatures[1:] = (fluxes[1:] / (sigma*e_atm)) **(1/4)
    
    #Ground
    temperatures[0] = (fluxes[0]/(sigma*e_ground)) **(1/4)
    
    
    results = modeldata(
        temperatures, 
        range(0,nlayers+1),
        nlayers, 
        solar,
        emissivity,
        albedo,
        nuclear_winter
        
        
        )
    return results

    
def main():
    
    #normaltemp = run_atmmodel(2, 1370, 0.3, 0.3, False)
    #nucleartemp = run_atmmodel(2, 1370, 0.3, 0.3, True)
    #plt.plot(normaltemp.layers, normaltemp.temperatures)
    #plt.plot(nucleartemp.layers, nucleartemp.temperatures)
    
#Question 3    
    emissivities = []
    surfacetemp = []
    for emissivity in np.arange(0.01, 1, 0.01):
        results = run_atmmodel(1, 1350, emissivity, 0.33, False)
        emissivities.append(emissivity)
        surfacetemp.append(results.temperatures[0])
    
    plt.plot(emissivities, surfacetemp)
    plt.title('Surface Temp related to Emissivity Values')
    plt.xlabel('Emissivity')
    plt.ylabel('Surface Temp(K)')
    plt.show()
    
    
    nlayers1 = []
    surfacetemp = []
    for nlayers in np.arange(2, 50, 1):
        results = run_atmmodel(nlayers, 1350, 0.255, 0.33, False)
        nlayers1.append(nlayers)
        surfacetemp.append(results.temperatures[0])
        
    plt.plot(surfacetemp, nlayers1)
    plt.title('Surface Temp related to number of layers')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Number of atm.layers')
    plt.show()
    
    
#Question 4 
    # number of layers is independent variable changed in for loop
    # solar irradiance set to 2600 W/m^-2
    # emissivity set to 1
    # albedo set to 0, no shortwave radiation reflected
    
    nlayers1 = []
    surfacetemp1 = []
    for nlayers in np.arange(2, 50, 1):
        results = run_atmmodel(nlayers, 2600, 1, 0.33, False)
        nlayers1.append(nlayers)
        surfacetemp1.append(results.temperatures[0])
    
    plt.plot(nlayers1, surfacetemp1)
    plt.axhline(700, color='r', linestyle='--')
    plt.title('Venus surface temp related to number of layers')
    plt.ylabel('Temperature (K)')
    plt.xlabel('Number of atm.layers')
    plt.show()

    #Find numbers of layers that brings us closest to Venus Temp.
    layergoal = (min(range(len(surfacetemp1)),
              key=lambda i: abs(surfacetemp1[i]-700)))+2
    
    print(layergoal)
   
    



    #Question 5 - Nuclear Winter = True
        # number of layers set to 5
        # solar irradiance set to 1370 W/m^-2
        # emissivity set to 0.5
        # albedo = 0 as nuclear_winter=True
    
    nuclear_results = run_atmmodel(5, 1350, 0.5, 0.33, True)
    
    plt.plot(nuclear_results.temperatures, nuclear_results.layers)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Number of atm.layers')
    plt.show()
    

    
if __name__ == "__main__":
    main()












