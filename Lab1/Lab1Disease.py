# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 16:05:02 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

area_cmap = ListedColormap(['blue', 'green', 'red', 'black'])

#4 is dead
#1 is had it and immune
#2 is healthy
#3 is sick

nx, ny, numiters = 30, 20, 300 # Number of cells in X and Y direction and # of interations
prob_spread = 0.65 # OVERALL CHANCE OF SPREAD FROM ONE SICK PERSON TO A HEALTHY PERSON



def disease(prob_fatality,prob_immune):
    # Create an initial grid, set all values to "2"(healthy). dtype sets the value
    # type in our array to integers only.
    area  =  np.zeros([numiters,ny,nx], dtype=int) + 2
    numberimmune = 0
    numberdead = 0
    # chance of any cell starting as immune person
    # Loop in the "x" direction:
    for i in range(nx):
        # Loop in the "y" direction:
        for j in range(ny):
            ## Perform logic here:
            if np.random.rand() < prob_immune:
                area[0,j, i] = 1
                numberimmune = numberimmune + 1
                
    # Set the center cell to "sick":
    area [0,ny//2,nx//2] = 3
    
    fig, ax = plt.subplots(1,1)
    ax.pcolor(area[0,:,:], cmap=area_cmap, vmin=1, vmax=4)
    plt.close()
    
    #disease spread loop
    for k in range(1, numiters):
        
        area[k,:,:] = area[k-1,:,:]
        
        for i in range(nx):
        
            for j in range(ny):
                
                
                # if person was already sick, move on to possibility of spread
                if area[k-1,j,i]==3:
            
                    #going left, if not an edge and cell is healthy, spread has a chance
                    if i!=0 and area[k-1,j,i-1]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            area[k,j, i-1] = 3
    
                    #going down, if not an edge and cell is healthy, spread has a chance
                    if j!=0 and area[k-1,j-1,i]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            area[k,j-1,i] = 3
                        
                    #going right, if not an edge and cell is healthy, spread has a chance
                    if i!=nx-1 and area[k-1,j,i+1]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            area[k,j, i+1] = 3
                        
                    #going up, if not an edge and cell is healthy, spread has a chance
                    if j!=ny-1 and area[k-1,j+1,i]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            area[k,j+1, i] = 3
                            
                    # current sick cell becomes immune or dead
                    deathchance = np.random.rand()
                    if deathchance < prob_fatality:
                        area[k,j,i] = 4
                        numberdead = numberdead + 1
                    else: 
                        area[k,j,i] = 1 
                        numberimmune = numberimmune + 1
                    
        #plot current conditions
        # fig, ax = plt.subplots(1,1)
        # ax.pcolor(area[k,:,:], cmap=area_cmap, vmin=1, vmax=4)
        # plt.close()
        
        # end loop if there are no more cells are sick and output # of iterations
        # it took to get there and number/percentage of cells that are dead
        if 3 not in area[k,:,:]:
            time = k
            percentdead = numberdead/(nx*ny)
                    
            break

    return time+1, numberdead, percentdead
        
#initialize a range of probimmune and probfatality values
immunerange = np.arange(0,1.1,0.05)
fatalityrange = np.arange(0,1.1,0.05)

#initialize array of outputs from change in probfatality values
timearray = []
numberdeadarray = []
fatalityarray = []

    #loop range of probfatality values through disease spread function
for i in fatalityrange:
    timeend, numberofdead, percentagedead = disease(i,0.01)
    timearray.append(timeend)
    numberdeadarray.append(numberofdead)
    fatalityarray.append(percentagedead)
    
    
plt.figure(figsize=(12,8))
plt.plot(fatalityrange,fatalityarray)
plt.title('Fraction of cells that are dead with 1% of starting cells immune', fontsize=18)
plt.xlabel('Probability of Fatality', fontsize=20)
plt.ylabel('Fraction of cells that are dead after spread', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


#initialize array of outputs from change in probimmune values
timearray2 = []
numberimmunearray = []
immunearray = []

    #loop range of probimmune values through disease spread function
for i in immunerange:
    timeend2, numberofdead, percentagedead = disease(0.1,i)
    timearray2.append(timeend2)
    numberimmunearray.append(numberofdead)
    immunearray.append(percentagedead)
    
    
plt.figure(figsize=(12,8))
plt.plot(immunerange,immunearray)
plt.title('Fraction of cells that are dead with a 10% chance of death from the disease', fontsize=20)
plt.xlabel('Fraction of cells that begin with immunity from death (Vaccine)', fontsize=18)
plt.ylabel('Fraction of cells that are dead after spread', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

