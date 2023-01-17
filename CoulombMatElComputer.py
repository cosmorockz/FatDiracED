'''This code evaluates the matrix elements for two-particle 
Coulomb interactions. Given any Dirac cutoff (specifies )

'''

import numpy as np
from sympy.physics.wigner import wigner_3j

cutoff = 1 # Landau level index cutoff
Q = 2 # Dirac magnetic Monopole charge

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

TotalStates = k + 1 # Total number of states

# Finding the allowed combinations for Coulomb interaction

CoulombIntIndices = [] # All the allowed four-state 
# indices for Coulomb interaction

for key1, value1 in StatesAll.items():
    for key2, value2 in StatesAll.items():
        for key3, value3 in StatesAll.items():
            for key4, value4 in StatesAll.items():
                m1 = value1['m']
                m2 = value2['m']
                m3 = value3['m']
                m4 = value4['m']
                if abs(m1+m2+m3+m4) <= 1e-8:
                    CoulombIntIndices.append((key1,key2,key3,key4))


def CoulombStrength(Ind):

    ind1,ind2, ind3, ind4 = Ind

    state1 = StatesAll[ind1]
    state2 = StatesAll[ind2]
    state3 = StatesAll[ind3]
    state4 = StatesAll[ind4]

    n1 = state1['n']
    n2 = state2['n']
    n3 = state3['n']
    n4 = state4['n']

    la1 = state1['lambda']
    la2 = state2['lambda']
    la3 = state3['lambda']
    la4 = state4['lambda']

    m1 = round(state1['m'],1)
    m2 = round(state2['m'],1)
    m3 = round(state3['m'],1)
    m4 = round(state4['m'],1)

    l1 = round(n1 + Q - 0.5,1)
    l2 = round(n2 + Q - 0.5,1)
    l3 = round(n3 + Q - 0.5,1)
    l4 = round(n4 + Q - 0.5,1)

    IntStr = 0 # Interaction Strength

    for l in range(int(min(l1+l4,l2+l3))):
        # print(l)
        for m in range(-l,l+1):
        
            UpUp = (-1)**(int(Q+0.5-m1-m+l1+l+l4+Q+0.5-m2+l2+l+l3)) \
                * np.sqrt((2*l1+1)*(2*l+1)*(2*l4+1)/(4*np.pi)) \
                    * np.sqrt((2*l2+1)*(2*l+1)*(2*l3+1)/(4*np.pi)) \
                        * wigner_3j(l1,l,l4,m1,m,-m4) \
                            * wigner_3j(l1,l,l4,round(-Q-0.5,1),0,round(Q+0.5,1)) \
                                * wigner_3j(l2,l,l3,m2,-m,-m3) \
                                    * wigner_3j(l2,l,l3,round(-Q-0.5,1),0,round(Q+0.5,1))
            
            # print(UpUp)
            
            IntStr += UpUp

    return IntStr









