'''This code evaluates the matrix elements for two-particle 
Coulomb interactions. Given any Dirac cutoff (specifies )

'''

import numpy as np

cutoff = 1 # Landau level index cutoff
Q = 4 # Dirac magnetic Monopole charge

l_max = Q - 0.5 + cutoff # Maximum Landau level for the given Magnetic
# field and "cutoff"

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


MM = [-l_max + i for i in range(int(2*l_max + 1))] # All the allowed m's


# Making a Dictionary for all the allowed states

StatesAll = {} # Dictionary of all the allowed states

k = 0 # Key for the dictionary


for m in MM:
    AllowedLevels = Bset(Q,cutoff,m)
    for level in AllowedLevels:
        la = np.sign(level) # Lambda
        n = np.abs(level) # Landau level index
        StatesAll[k] = {}
        StatesAll[k]['m'] = m 
        StatesAll[k]['n'] = n 
        StatesAll[k]['lambda'] = la
        k += 1








