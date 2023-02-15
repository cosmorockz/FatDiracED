import re
import glob
from header import *
from scipy import sparse
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
from KineticPart import KineticOperator
from CoulombPart import CoulombOperator


BS_Folder = "./BasisStates"+str(cutoff)+"Q"+str(Q)+"/NParticles" \
    +str(NParticles)

if os.path.isdir(BS_Folder):
    pass
else:
    os.system("python BasisStates.py")

AllFiles = glob.glob(BS_Folder+"/*.dat")
LTotal = []

marker1 = "AM"
marker2 = ".dat"
for ff in AllFiles:
    regexPattern = marker1 + '(.+?)' + marker2
    am = float(re.search(regexPattern,ff).group(1))
    LTotal.append(am)

def IndexMapping(inds,sortedInds):
    indDict = {}
    # try:
    for i,ind in enumerate(inds):
        indDict[ind] = sortedInds[i]
    # except TypeError:
    #     indDict[inds.tolist()] = sortedInds[0]
    return indDict

AllEigvals = []

for l in LTotal:

    AM_File = BS_Folder + "/AM" + str(l) + ".dat"
    # print(AM_File)
    AllInds = np.loadtxt(AM_File,unpack=True,dtype=int)
    try:
        ll = len(AllInds)
    except TypeError:
        AllInds = [AllInds.tolist()]
    sIndices = np.argsort(AllInds)
    oIndices = np.argsort(sIndices)
    # print(oIndices)
    indDict = IndexMapping(AllInds, oIndices)
    rows = []
    cols = []
    values = []

    # # Kinetic part
    for ind in AllInds:
        ke = KineticOperator(ind)
        rowInd = indDict[ind]
        rows.append(rowInd)
        cols.append(rowInd)
        values.append(ke)
    
    # Coulomb part
    for ind in AllInds:
        MatrixElements = CoulombOperator(ind)
        rowInd = indDict[ind]
        for (stateF,element) in MatrixElements:
            rows.append(rowInd)
            cols.append(indDict[stateF])
            values.append(element)
    
    MatrixSize = len(AllInds)

    # print(H)

    # print(values)

    try:
        H = sparse.csr_matrix((values,(rows,cols)))
        if MatrixSize <= nEigvals + 1:
            eigvals, eigvecs = np.linalg.eigh(H.todense())
        else:
            eigvals, eigvecs = eigsh(H,k=nEigvals,which="SA")
    except ValueError:
        eigvals = [0]
        

    # print(H)
    AllEigvals.append(eigvals)
    
LL = []
EE = []
for i in range(len(AllEigvals)):
    for eigval in AllEigvals[i]:
        LL.append(LTotal[i])
        EE.append(eigval)

plt.scatter(LL,EE)
plt.show()








