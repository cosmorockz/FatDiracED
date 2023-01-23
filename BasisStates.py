'''This script finds and saves states with the same total
angular momentum, for a given monopole charge and cutoff'''

from header import *

def decimalToBinary(n):
    # Converts a given number in decimal representation to a binary number
    return "{0:b}".format(int(n))

def ContainerArray(n,ArraySize):
    # Given a decimal number, this function outputs a binary number
    # in an array of fixed legth.
    # For example, if n=4 and ArraySize=6,
    # The output will be [0,0,0,1,0,0]
    Container = [0 for i in range(ArraySize)]
    BinaryNumber = decimalToBinary(n)
    LenBinary = len(BinaryNumber)
    for i in range(LenBinary):
        Container[ArraySize-1-i] = int(BinaryNumber[LenBinary-1-i])
    return Container

def FockBasisToDecimal(state):
    BinaryString = ''.join(str(e) for e in state)
    n = int(BinaryString,2)
    return n

def BasisGenerator(Nlevels):
    # Given a Hilbert space size, generates all the basis vectors that spans the space
    MaxN = 2**Nlevels # Highest number in the basis in decimal (By highest we mean the fully filled system)
    basis = []
    for i in range(MaxN):
        basis.append(ContainerArray(i,Nlevels))
    return basis

def Bset(Q,cutoff,m):
    # Set of quantum numbers if magentic field and cutoff is specified for a given m
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

def FockBasisToOperator(BasisElement):
    # Inputs a state in the Fock basis e.g. [0,1,0,1]
    # Outputs it as [2,4]
    basis = []
    for i,element in enumerate(BasisElement):
        if element == 1:
            basis.append(i)
    return basis

def OperatorToFockBasis(state):
    ArraySize = TotalSPStates
    Container = [0 for i in range(ArraySize)]
    for element in state:
        Container[element] = 1
    return Container

def DecimalToFockBasis(n):
    Binary = "{0:b}".format(int(n))
    ArraySize = TotalSPStates
    Container = [0 for i in range(ArraySize)]
    LenBinary = len(Binary)
    for i in range(LenBinary):
        Container[ArraySize-1-i] = int(Binary[LenBinary-1-i])
    return Container

def DecimalToFockBasis(n):
    Binary = "{0:b}".format(int(n))
    ArraySize = TotalSPStates
    Container = [0 for i in range(ArraySize)]
    LenBinary = len(Binary)
    for i in range(LenBinary):
        Container[ArraySize-1-i] = int(Binary[LenBinary-1-i])
    return Container

mm = [-l_max + i for i in range(int(2*l_max+1))]
MM = [] # Angular momentum of states
AllStates = [] # Contains tuples of all levels
for m in mm:
    LevelSet = Bset(Q,cutoff,m)
    NN = len(LevelSet)
    for i in range(NN):
        MM.append(m)
    for element in LevelSet:
        AllStates.append((element,m))

# Finds and saves states with the same total angular momentum
data = {}
AMlist = []
ArraySize = len(MM)
for n in range(0,LargestNumber):
    Number = ContainerArray(n,ArraySize)
    AngularMomentum = np.dot(MM,Number)
    AMlist.append(AngularMomentum)
    key = str(AngularMomentum)
    if key in data:
        data[key].append(n)
    else:
        data[key] = [n]

AMs = list(set(AMlist)) # List of all the total angular momentum

del AMlist
gc.collect()

AMfolder = "./DecimalAMCorrCutoff"+str(cutoff)+"Q"+str(Q)
try:
    shutil.rmtree(AMfolder)
except OSError:
    pass

os.mkdir(AMfolder)

AMsector = [] # Angular momentum sector
LenAMsector = [] # Total number of elements in each angular
# momentum sector
for key, values in data.items():
    np.savetxt(AMfolder+"/AM"+str(round(float(key),1))\
        +".dat", np.c_[values], fmt="%i") #, np.array(values).astype(int))
    AMsector.append(float(key))
    LenAMsector.append(len(values))

del data
gc.collect()

np.savetxt(AMfolder+"/AMinfo.dat",np.c_[AMsector,LenAMsector],fmt="%.1f %i")





