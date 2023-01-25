'''This script finds and saves states with the same total
angular momentum, for a given monopole charge and cutoff'''

import itertools
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

def BasisGeneratorNumberConserving(BasisElement):
    # Given a basis element, this function generates all the other basis elements 
    # of same angular momentum and same particle number
    return list(set(list(itertools.permutations(BasisElement))))

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

ParallelRoutes = pFactor * NCPU # Size of the parallized parts

def ChunkCreator(iPoint,fPoint,pRoutes):
    # Given a starting point, ending point and the number of chunks,
    # this function returns all the chunks
    chunks = [] # Stores the starting and ending point 
    #of all the chunks
    chunkSize = int((fPoint - iPoint) / pRoutes)
    for i in range(pRoutes):
        chunkSP = 0 + i * chunkSize # Starting point of the chunk
        chunkEP = chunkSP + chunkSize - 1 # Ending point of the chunk
        chunks.append((chunkSP,chunkEP))
    return chunks

AllChunks = ChunkCreator(0,LargestNumber,ParallelRoutes)
# Contains all The parallelizable chunks

BS_folder = "./BasisStates"+str(cutoff)+"Q"+str(Q)
# Folder for the basis states

# Checks if there is the Basis states folder and creates if
# it does not exist
if os.path.isdir(BS_folder) == True:
    pass
else:
    os.mkdir(BS_folder)

NP_folder = BS_folder + "/NParticles"+str(NParticles)
# Folder inside BS_folder with a given particle number

if os.path.isdir(NP_folder) == True:
    pass
else:
    os.mkdir(NP_folder)

def AM_Assigner(chunk):
    # Given a Chunk of decimal numbers, this function forms a
    # dictionary of them on the basis of angular momentum.
    # Then it Saves them in different files according
    # to their angular momentum. Notice that
    # this function only saves the decimal whose corresponding
    # binary representation has the required particle number

    chunkSP = chunk[0]
    chunkEP = chunk[1]

    data = {}

    for n in range(chunkSP, chunkEP):
        
        BNumber = ContainerArray(n,TotalSPStates)
        
        if sum(BNumber) != NParticles:
            continue
        else:
            key = str(np.dot(BNumber,MM))
            if key in data:
                data[key].append(n)
            else:
                data[key] = [n]
    
    for key, values in data.items():
        
        filename = NP_folder + "/AM" + str(round(float(key),1)) + ".dat"
        # if os.path.isfile(filename) == True:
        with lock:
            with open(filename, 'a+') as fd:
                for value in values:
                    fd.write("%i \n" % value)
        # else:
        #     np.savetxt(filename,np.c_[values],fmt="%i")
    
    del data
    gc.collect()

    return None

with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(AM_Assigner,AllChunks)
        





