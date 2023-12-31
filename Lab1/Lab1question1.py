# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 20:40:18 2023

@author: dscho

This code is the same as Lab1Fire.py
It is exclusively used to illustrate question #1 on lab 1.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

forest_cmap = ListedColormap(['tan', 'darkgreen', 'crimson'])

#1 is bare
#2 is forest
#3 is on fire

nx, ny, numiters = 11, 5, 300 # Number of cells in X and Y direction and # of interations

prob_start = 0.005 # Chance of cell being on fire at the beginning.

def fire(prob_spread,prob_bare):

    # Create an initial grid, set all values to "2". dtype sets the value
    # type in our array to integers only.
    forest  =  np.zeros([numiters,ny,nx], dtype=int) + 2
    #initialize number of bare cells
    numberbare = 0;
    # chance of any cell starting bare
    # Loop in the "x" direction:
    for i in range(nx):
        # Loop in the "y" direction:
        for j in range(ny):
            ## Perform logic here:
            if np.random.rand() < prob_bare:
                forest[0,j, i] = 1
                numberbare = numberbare+1
                
    # chance of a cell beginning on fire
    # Loop in the "x" direction:
    for i in range(nx):
        # Loop in the "y" direction:
        for j in range(ny):
            ## Perform logic here:
            if np.random.rand() < prob_start:
                forest[0,j, i] = 3
                
    # Set the center cell to "burning":
    forest [0,ny//2,nx//2] = 3
    
    
    
    #originate plot before fire spreads
    fig, ax = plt.subplots(1,1)
    ax.pcolor(forest[0,:,:], cmap=forest_cmap, vmin=1, vmax=3)
    
    #fire spread loop through number of iterations
    for k in range(1, numiters):
        
        forest[k,:,:] = forest[k-1,:,:]
        
        for i in range(nx):
        
            for j in range(ny):
                
                
                # if the forest was already on fire, move on to possibility of spread
                if forest[k-1,j,i]==3:
            
                    #going left, if not an edge and cell is a forest, spread has a chance
                    if i!=0 and forest[k-1,j,i-1]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            forest[k,j, i-1] = 3
    
                    #going down, if not an edge and cell is a forest, spread has a chance
                    if j!=0 and forest[k-1,j-1,i]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            forest[k,j-1,i] = 3
                        
                    #going right, if not an edge and cell is a forest, spread has a chance
                    if i!=nx-1 and forest[k-1,j,i+1]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            forest[k,j, i+1] = 3
                        
                    #going up, if not an edge and cell is a forest, spread has a chance
                    if j!=ny-1 and forest[k-1,j+1,i]==2:
                        firespread =np.random.rand()
                        if firespread <= prob_spread:
                            forest[k,j+1, i] = 3
                            
                    # current on fire cell becomes bare
                    forest[k,j,i] = 1
                    numberbare = numberbare +1
                    
        
        #plot current conditions
        fig, ax = plt.subplots(1,1)
        ax.pcolor(forest[k,:,:], cmap=forest_cmap, vmin=1, vmax=3)
        
        # end loop if there are no more cells on fire and output # of iterations
        # it took to get there and number/percentage of cells that are bare
        if 3 not in forest[k,:,:]:
            time = k
            percentbare = numberbare/(nx*ny)
                    
            break

    return time+1, numberbare, percentbare

fire(1,0.01)