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

#1 is bare
#2 is forest
#3 is on fire

nx, ny, numiters = 3, 3, 3 # Number of cells in X and Y direction.
prob_spread = 1.0 # Chance to spread to adjacent cells.
prob_bare = 0.0 # Chance of cell to start as bare patch.
prob_start = 0.0 # Chance of cell to start on fire.


# Create an initial grid, set all values to "2". dtype sets the value
# type in our array to integers only.
forest  =  np.zeros([ny,nx,numiters], dtype=int) + 2
forest([numiters], dtype=int) ==0

# Set the center cell to "burning":
forest [1,1] = 3

#fire spread loop
for k in range(1, numiters):
    for i in range(nx):
        for j in range(ny):
            # Roll our "dice" to see if we get a bare spot:
            if forest(k-1,j,i)==3:
                
                if forest(k-1,j,i-1)==2:
                    
                    forest[j, i-1] = 3
            
            
            
            
            
isbare = np.random.rand(nx, ny)
isbare = isbare < prob_bare
forest[isbare] = 1


# Loop in the "x" direction:
for i in range(nx):
    # Loop in the "y" direction:
    for j in range(ny):
        ## Perform logic here:
        forest[i, j] = 5







