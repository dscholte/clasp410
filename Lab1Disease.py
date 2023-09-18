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

nx, ny, numiters = 30, 20, 30 # Number of cells in X and Y direction and # of interations
prob_spread = 0.6 # Chance to spread to adjacent people.
prob_immune = 0.01 # Chance of starting with immunity with vaccine
prob_start = 0.001 # Chance of cell randomly becoming infected.
#prob_fatality = 0.1 # Chance that an infected person dies (3 becomes 4)


def diseasespr(prob_fatality):
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
    
    #disease spread loop
    for k in range(1, numiters):
        
        area[k,:,:] = area[k-1,:,:]
        
        for i in range(nx):
        
            for j in range(ny):
                
                # random chance of a person getting sick 
                if area[k-1,j,i]==2:
                    firestart = np.random.rand()
                    if firestart <= prob_start:
                        area[k,j,i] = 3
                
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
        fig, ax = plt.subplots(1,1)
        ax.pcolor(area[k,:,:], cmap=area_cmap, vmin=1, vmax=4)
        
        # end loop if there are no more cells on fire and output # of iterations
        # it took to get there and number/percentage of cells that are bare
        if 3 not in area[k,:,:]:
            time = k
            percentdead = numberdead/(nx*ny)
                    
            break

    return time+1, numberdead, percentdead
        
        