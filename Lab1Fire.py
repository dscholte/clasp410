# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 10:08:13 2023

@author: dscho
"""

#!/usr/bin/env python3
'''
This file contains tools and scripts for completing Lab 1 for CLaSP 410.
To reproduce the plots shown in the lab report, do this...
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

forest_cmap = ListedColormap(['tan', 'darkgreen', 'crimson'])

#1 is bare
#2 is forest
#3 is on fire

nx, ny, numiters = 3, 3, 4 # Number of cells in X and Y direction.
prob_spread = 1.0 # Chance to spread to adjacent cells.
prob_bare = 0.0 # Chance of cell to start as bare patch.
prob_start = 0.0 # Chance of cell to start on fire.


# Create an initial grid, set all values to "2". dtype sets the value
# type in our array to integers only.
forest  =  np.zeros([numiters,ny,nx], dtype=int) + 2

 
# Loop in the "x" direction:
for i in range(nx):
    # Loop in the "y" direction:
    for j in range(ny):
        ## Perform logic here:
        if np.random.rand() < prob_bare:
            forest[0,i, j] = 1
            
# Set the center cell to "burning":
forest [0,1,1] = 3

fig, ax = plt.subplots(1,1)
ax.pcolor(forest[0,:,:], cmap=forest_cmap, vmin=1, vmax=3)

#fire spread loop
for k in range(1, numiters):
    
    forest[k,:,:] = forest[k-1,:,:]
    for i in range(nx):
        
        
        for j in range(ny):
            
            
            if forest[k-1,j,i]==3:
        
                
                #going left
                if i!=0 and forest[k-1,j,i-1]==2:
                    
                    forest[k,j, i-1] = 3

                #going down
                if j!=0 and forest[k-1,j-1,i]==2:
                    
                    forest[k,j-1,i] = 3
                    
                #going right 
                if i!=nx-1 and forest[k-1,j,i+1]==2:
                    
                    forest[k,j, i+1] = 3
                    
                #going up
                if j!=ny-1 and forest[k-1,j+1,i]==2:
                    
                    forest[k,j+1, i] = 3

                forest[k,j,i] = 1
    fig, ax = plt.subplots(1,1)
    ax.pcolor(forest[k,:,:], cmap=forest_cmap, vmin=1, vmax=3)






