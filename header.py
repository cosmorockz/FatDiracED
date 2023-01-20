'''This script contains all variables for the 
exact diagonalization code'''

import gc
import os
import shutil
import numpy as np
import concurrent.futures


cutoff = 1 # Landau level index cutoff
Q = 3 # Dirac magnetic Monopole charge
StrengthTol = 1e-8

l_max = Q - 0.5 + cutoff # Maximum Landau level for the given Magnetic
# field and "cutoff"

TotalSPStates = int((4*cutoff + 2) * Q + 2 * cutoff \
    * (cutoff + 1)) # Total number of single particle states

LargestNumber = 2**TotalSPStates

def Bset(Q,cutoff,m):
    # Set of quantum numbers if magentic 
    # field and cutoff is specified for a given m
    levels = []
    if ((Q*2) % 2) < 1e-8:
#         M = max(int(abs(m)+0.5),)
        for i in range(-cutoff,cutoff+1):
            l = Q + np.abs(i) - 1/2
            if ((np.abs(m)*2) % 2) > 0.999:
                if np.abs(m) <= np.abs(l):
                    levels.append(i)
    else:
        for i in range(-cutoff,cutoff+1):
            l = Q + np.abs(i) - 1/2
            if ((np.abs(m)*2) % 2) < 1e-8:
                if np.abs(m) <= np.abs(l):
                    levels.append(i)
    return levels




