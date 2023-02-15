from header import *
from CommonFunctions import *

Lz = 0.0

BS_Folder = "./BasisStates"+str(cutoff)+"Q"+str(Q)+"/NParticles" \
    +str(NParticles)

AM_File = BS_Folder + "/AM" + str(Lz) + ".dat"

AllInds = np.loadtxt(AM_File,unpack=True,dtype=int)

def LpLn(state, n1, n2):
    # Acts as l_n2^+ l_n1^-

    level1, m1 = AllStates[n1]
    level2, m2 = AllStates[n2]

    BNumber = ContainerArray(state, TotalSPStates)

    if (BNumber[n1] == 0) or (BNumber[n2] == 0):
        return (-1,0)
    else:
        m1new = round(m1-1,1)
        m2new = round(m2+1,1)
        
        try:
            n1new = AllStatesInv[(level1,m1new)]
        except KeyError:
            return (-1,0)

        try:
            n2new = AllStatesInv[(level2,m2new)]
        except KeyError:
            return (-1,0)

        if BNumber[n1new] == 0:
            BNumber[n1] = 0
            BNumber[n1new] = 1
            if BNumber[n2] == 0:
                return (-1,0)
            else:
                if BNumber[n2new] == 0:
                    BNumber[n2] = 0
                    BNumber[n2new] = 1
                    StateF = FockBasisToDecimal(BNumber)
                    l1 = round(Q-0.5+abs(level1),1)
                    l2 = round(Q-0.5+abs(level2),1)
                    ValueF = np.sqrt(l1*(l1+1) - m1*(m1-1)) \
                        * np.sqrt(l2*(l2+1) - m2*(m2+1))
                    return (StateF,ValueF)

                else:
                    return (-1,0)
        
        else:
            return (-1,0)

def LnLp(state,n1,n2):
    # Acts as l_n2^- l_n1^+

    level1, m1 = AllStates[n1]
    level2, m2 = AllStates[n2]

    BNumber = ContainerArray(state, TotalSPStates)

    if (BNumber[n1] == 0) or (BNumber[n2] == 0):
        return (-1,0)
    else:
        m1new = round(m1+1,1)
        m2new = round(m2-1,1)

        try:
            n1new = AllStatesInv[(level1,m1new)]
        except KeyError:
            return (-1,0)

        try:
            n2new = AllStatesInv[(level2,m2new)]
        except KeyError:
            return (-1,0)

        if BNumber[n1new] == 0:
            BNumber[n1] = 0
            BNumber[n1new] = 1
            if BNumber[n2] == 0:
                return (-1,0)
            else:
                if BNumber[n2new] == 0:
                    BNumber[n2] = 0
                    BNumber[n2new] = 1
                    StateF = FockBasisToDecimal(BNumber)
                    l1 = round(Q-0.5+abs(level1),1)
                    l2 = round(Q-0.5+abs(level2),1)
                    ValueF = np.sqrt(l1*(l1+1) - m1*(m1+1)) \
                        * np.sqrt(l2*(l2+1) - m2*(m2-1))
                    return (StateF,ValueF)

                else:
                    return (-1,0)
        
        else:
            return (-1,0)

def LzLz(state,n1,n2):
    # Acts as l_n2^z l_n1^z

    BNumber = ContainerArray(state,TotalSPStates)

    if (BNumber[n1] == 0) or (BNumber[n2] == 0):
        return (-1,0)
    else:
        level1, m1 = AllStates[n1]
        level2, m2 = AllStates[n2]
        return (state,m1+m2)

L2OpRows = []
L2OpCols = []

L2OpData = []

for ind in AllInds:
    LzLzValue = 0
    for i in range(TotalSPStates):
        for j in range(TotalSPStates):
            state1, value1 = LpLn(ind,j,i) # l_i^+ l_j^-
            state2, value2 = LnLp(ind,j,i) # l_i^- l_j^+
            state3, value3 = LzLz(ind,j,i) # l_i^z l_j^z

            if state1 != -1:
                L2OpRows.append(ind)
                L2OpCols.append(state1)
                L2OpData.append(value1/2)
            
            if state2 != -1:
                L2OpRows.append(ind)
                L2OpCols.append(state2)
                L2OpData.append(value2/2)

                print(ind, i, j, state2)
            
            if state3 != -1:
                LzLzValue += value3
    if abs(LzLzValue) > StrengthTol:
        L2OpRows.append(ind)
        L2OpCols.append(ind)
        L2OpData.append(LzLzValue)







