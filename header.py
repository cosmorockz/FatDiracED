'''This script contains all variables for the 
exact diagonalization code'''

import gc
import os
import shutil
import threading
import numpy as np
import concurrent.futures

cutoff = 1 # Landau level index cutoff
Q = 1 # Dirac magnetic Monopole charge
StrengthTol = 1e-8
NCPU = os.cpu_count() # Number of CPUs in the current computer
pFactor = 4 # Parallization Factor; Basis States are divided into
# pFactor * NCPU groups

kinConst = 2 # Strength of the kinetic part
intConst = 4 # Strength of the interacting part

l_max = Q - 0.5 + cutoff # Maximum Landau level for the given Magnetic
# field and "cutoff"

TotalSPStates = int((4*cutoff + 2) * Q + 2 * cutoff \
    * (cutoff + 1)) # Total number of single particle states

NParticles = 2 # Total number of electrons in the system

LargestNumber = 2**TotalSPStates

lock = threading.Lock()

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




