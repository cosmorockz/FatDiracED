from header import *
from CommonFunctions import *

def CreationOperator(state,Op):
    BNumber = ContainerArray(state,TotalSPStates)
    if BNumber[Op] == 1:
        return (0,-1)
    else:
        Sign1 = (-1)**(sum(BNumber[:Op]))
        BNumber[Op] = 1
        return (Sign1,FockBasisToDecimal(BNumber))

def AnnihilationOperator(state,Op):
    BNumber = ContainerArray(state,TotalSPStates)
    if BNumber[Op] == 0:
        return (0,-1)
    else:
        Sign1 = (-1)**(sum(BNumber[:Op]))
        BNumber[Op] = 0
        return (Sign1,FockBasisToDecimal(BNumber))



