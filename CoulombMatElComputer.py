'''This code evaluates the matrix elements for two-particle 
Coulomb interactions. Given any Dirac cutoff (specifies )

'''

from header import *
from sympy.physics.wigner import wigner_3j
from header import intConst

print(intConst)


MM = [-l_max + i for i in range(int(2*l_max + 1))] # All the allowed m's


# Making a Dictionary for all the allowed states

StatesAll = {} # Dictionary of all the allowed states

k = 0 # Key for the dictionary


for m in MM:
    AllowedLevels = Bset(Q,cutoff,m)
    for level in AllowedLevels:
        la = np.sign(level) # Lambda
        if la == 0:
            la = 1
        n = np.abs(level) # Landau level index
        StatesAll[k] = {}
        StatesAll[k]['m'] = m 
        StatesAll[k]['n'] = n 
        StatesAll[k]['lambda'] = la
        k += 1

TotalStates = k # Total number of states

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

    # Given the Indices for four operators, this function
    # computes the Coulomb interaction strength between them

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

    # print(n1,n2,n3,n4)
    # print(m1,m2,m3,m4)
    # print(la1,la2,la3,la4)
    # print(l1,l2,l3,l4)

    # IntStr = 0 # Interaction Strength

    # print(la1,la2,la3,la4)

    UpUpInt = 0
    UpDownInt = 0
    DownUpInt = 0
    DownDownInt = 0

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
            
            UpDown = (-1)**(int(Q+0.5-m1-m+l1+l+l4+Q-0.5-m2+l2+l+l3)) \
                * np.sqrt((2*l1+1)*(2*l+1)*(2*l4+1)/(4*np.pi)) \
                    * np.sqrt((2*l2+1)*(2*l+1)*(2*l3+1)/(4*np.pi)) \
                        * wigner_3j(l1,l,l4,m1,m,m4) \
                            * wigner_3j(l1,l,l4,round(-Q-0.5,1),0,round(Q+0.5,1)) \
                                * wigner_3j(l2,l,l3,m2,-m,-m3) \
                                    * wigner_3j(l2,l,l3,round(-Q+0.5,1),0,round(Q-0.5,1))
            
            DownUp = (-1)**(int(Q-0.5-m1-m+l1+l+l4+Q+0.5-m2+l2+l+l3)) \
                * np.sqrt((2*l1+1)*(2*l+1)*(2*l4+1)/(4*np.pi)) \
                    * np.sqrt((2*l2+1)*(2*l+1)*(2*l3+1)/(4*np.pi)) \
                        * wigner_3j(l1,l,l4,m1,m,-m4) \
                            * wigner_3j(l1,l,l4,round(-Q+0.5,1),0,round(Q-0.5,1)) \
                                * wigner_3j(l2,l,l3,m2,-m,-m3) \
                                    * wigner_3j(l2,l,l3,round(-Q-0.5,1),0,round(Q+0.5,1))
            
            DownDown = (-1)**(int(Q-0.5-m1-m+l1+l+l4+Q-0.5-m2+l2+l+l3)) \
                * np.sqrt((2*l1+1)*(2*l+1)*(2*l4+1)/(4*np.pi)) \
                    * np.sqrt((2*l2+1)*(2*l+1)*(2*l3+1)/(4*np.pi)) \
                        * wigner_3j(l1,l,l4,m1,m,-m4) \
                            * wigner_3j(l1,l,l4,round(-Q+0.5,1),0,round(Q-0.5,1)) \
                                * wigner_3j(l2,l,l3,m2,-m,-m3) \
                                    * wigner_3j(l2,l,l3,round(-Q+0.5,1),0,round(Q-0.5,1))
            
            
            UpUpInt += UpUp
            UpDownInt += UpDown
            DownUpInt += DownUp
            DownDownInt += DownDown

    IntStr = (UpUpInt + la2*la3*UpDownInt + la1*la4*DownUpInt + la1*la2*la3*la4*DownDownInt)

    return Ind,float(IntStr)

NZIntIndices = [] # Indices for Non-zero Coulomb interaction Strength
IntStrength = [] # The corresponding interaction strength

# for Ind in CoulombIntIndices:
#     ind, strength = CoulombStrength(Ind)
#     if abs(strength) > 1e-8:
#         NZIntIndices.append(ind)
#         IntStrength.append(strength)

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(CoulombStrength,CoulombIntIndices)

    for result in results:
        ind, strength = result
        if abs(float(strength)) > StrengthTol:
            NZIntIndices.append(ind)
            IntStrength.append(strength)


np.savetxt("CoulombIntMatElemQ"+str(Q)+"cutoff"+str(cutoff)+".dat",\
    np.c_[NZIntIndices,IntStrength],fmt="%d    %d    %d    %d    %.8f")










