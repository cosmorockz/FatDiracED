from BasicOperators import *

def EnergyN(n):
    return kinConst * np.sign(n) * np.sqrt((Q+abs(n))**2 - Q**2)

KE = [EnergyN(n) for n,m in AllStates]

def KineticOperator(state):
    MatrixElement = 0
    for i in range(len(AllStates)):
        signA, stateA = AnnihilationOperator(state,i)
        if stateA == -1:
            continue
        else:
            signAC, stateAC = CreationOperator(stateA,i)
            signT = signA * signAC
            MatrixElement += signT * KE[i]
    return MatrixElement



