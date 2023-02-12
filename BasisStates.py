'''This script finds and saves states with the same total
angular momentum, for a given monopole charge and cutoff'''

from header import *
from CommonFunctions import *

ParallelRoutes = pFactor * NCPU # Size of the parallized parts

def ChunkCreator(iPoint,fPoint,pRoutes):
    # Given a starting point, ending point and the number of chunks,
    # this function returns all the chunks
    chunks = [] # Stores the starting and ending point 
    #of all the chunks
    chunkSize = int((fPoint - iPoint) / pRoutes)
    for i in range(pRoutes):
        chunkSP = 0 + i * chunkSize # Starting point of the chunk
        chunkEP = chunkSP + chunkSize # Ending point of the chunk
        chunks.append((chunkSP,chunkEP))
    if chunkEP < fPoint:
        chunkEP = fPoint
        chunks[pRoutes-1][1] = chunkEP
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


NP_folder = BS_folder + "/NParticles"+str(NParticles)
# Folder inside BS_folder with a given particle number

if os.path.isdir(NP_folder) == True:
    pass
else:
    os.mkdir(NP_folder)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(AM_Assigner,AllChunks)
        





