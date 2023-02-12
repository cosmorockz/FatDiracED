from header import *
from BasicOperators import *

CoulombFile = "CoulombIntMatElemQ" + str(Q) + "cutoff" + str(cutoff) + ".dat"

if os.path.exists(CoulombFile):
    pass
else:
    os.system("python CoulombMatElComputer.py")

Ind1, Ind2, Ind3, Ind4, CoulStrs = np.loadtxt(CoulombFile, usecols=(0,1,2,3,4), unpack=True)

Ind1 = Ind1.astype(int)
Ind2 = Ind2.astype(int)
Ind3 = Ind3.astype(int)
Ind4 = Ind4.astype(int)

def CoulombOperator(state):
    # 1, 2, 4, 3
    MatrixElements = []
    for i in range(len(Ind1)):
        sign3, state3 = AnnihilationOperator(state,Ind3[i])
        if state3 == -1:
            continue
        else:
            sign4, state4 = AnnihilationOperator(state3,Ind4[i])
            if state4 == -1:
                continue
            else:
                sign2, state2 = CreationOperator(state4,Ind2[i])
                if state2 == -1:
                    continue
                else:
                    sign1, stateF = CreationOperator(state2,Ind1[i])
                    if stateF == -1:
                        continue
                    else:
                        MatrixElements.append((stateF, intConst * sign1 * sign2 \
                            * sign3 * sign4 * CoulStrs[i]))
    return MatrixElements



