import re
import copy
import glob

import numpy as np

from header import *
from scipy import sparse
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
from KineticPart import KineticOperator
from CoulombPart import CoulombOperator

from header import intConst

# print(intConst)

# Your main code here
# Ensure that you use intConst wherever it is needed in your main code




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

# plt.scatter(LL,EE)
# plt.show()


def organize_spectrum(spec_dict, tolerance=1e-8):
    spec_dict = copy.deepcopy(spec_dict)

    LzTots = sorted(spec_dict.keys())

    for i, LzTot in enumerate(LzTots):
        for eigval in spec_dict[LzTot]:
            if eigval is not None:
                for j in range(i + 1, len(LzTots)):
                    for k, eigval_inner in enumerate(spec_dict[LzTots[j]]):
                        if eigval_inner is not None:
                            if abs(eigval - eigval_inner) < tolerance:
                                spec_dict[LzTots[j]][k] = None
                                break

    for key, value in spec_dict.items():
        filtered_list = [vals for vals in value if vals is not None]
        spec_dict[key] = filtered_list
    return spec_dict


def lists_to_dict(LL, EE):
    spec_dict = {}
    for i, LzTot in enumerate(LL):
        if LzTot not in spec_dict.keys():
            spec_dict[LzTot] = []
        spec_dict[LzTot].append(EE[i])

    return spec_dict

def dict_to_lists(spec_dict):

    spec_list = []
    L_list = []

    for L in spec_dict.keys():
        for spec in spec_dict[L]:
            L_list.append(L)
            spec_list.append(spec)

    return L_list, spec_list

spec_dict = lists_to_dict(LL, EE)

spec_dict = organize_spectrum(spec_dict)

L_list, spec_list = dict_to_lists(spec_dict)
L_list = -np.array(L_list)

directory = f'./Simulation_Data/cutoff_{cutoff}/Q_{Q}'
if not os.path.exists(directory):
    os.makedirs(directory)

int_vs_kin = np.round(intConst / kinConst, 8)
filename = f'spec_{int_vs_kin:.8f}.dat'
full_filepath = os.path.join(directory, filename)
np.savetxt(full_filepath, np.c_[L_list, spec_list])

Plot_Directory = f'./Plots/cutoff_{cutoff}/Q_{Q}'
if not os.path.exists(Plot_Directory):
    os.makedirs(Plot_Directory)
plot_filename = f'spec_{int_vs_kin:.8f}.pdf'
plot_filename = os.path.join(Plot_Directory, plot_filename)


plt.figure(1)
plt.scatter(L_list, spec_list)
plt.xlabel('L')
plt.ylabel('$U_{int}$/$L_{kin}$')
plt.grid(True)
plt.savefig(plot_filename)
plt.show()























