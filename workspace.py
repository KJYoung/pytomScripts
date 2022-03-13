from pytom.basic.files import recenterVolume, naivePDBParser, mmCIFParser
from pytom.basic.files import read
from pytom.simulation.EMSimulation import simpleSimulation
from pytom_volume import vol, initSphere
from pytom.basic.structures import WedgeInfo

# VKJY
import wget
import os.path
import re
import requests
from urllib.error import HTTPError
import numpy as np
import mrcfile
import math 
import random
import time
# ------------
# Custom volume generator from pdb, cif.
# Compact volume generaor is added.
#
#

rootDIR = "/cdata"
pdbDataDIR = f"{rootDIR}/pdbData"
singleParticleEMCuboidDIR = f"{rootDIR}/singleParticleEM_cuboid"
singleParticleEMCubeDIR = f"{rootDIR}/singleParticleEM_cube"
singleParticleMRCCuboidDIR = f"{rootDIR}/singleParticleMRC_cuboid"
singleParticleMRCCubeDIR = f"{rootDIR}/singleParticleMRC_cube"
scenarioDIR = f"{rootDIR}/scenario"
metadataFILE = f"{rootDIR}/metadata.txt"

global called
called =  0
def appendMetaDataln(metadata):
    fmeta = open(metadataFILE, "a")
    fmeta.write(metadata + "\n")
    fmeta.close()
##########################################################################################################################################################################
#  Section for PDB ID to volume.
##########################################################################################################################################################################
def em2mrc(filename,newfilename):
    # Not USED : incompatible with mrcfile
    from pytom_volume import read
    from pytom.tools.files import checkFileExists,checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('EM file not found! ',filename)

    emfile = read(filename)
    emfile.write(newfilename,'mrc')

def atomList2emCube(atomList, pixelSize, densityNegative=False, resolutionFactor=None, verbose=False):
    atoms = []
    for atom in atomList:
        if not atom.getAtomType() in atoms:
            atoms.append(atom.getAtomType())
       
    print(atoms)
    """
    atomList2emCube : generate cube em file containing single particles.
    @param atomList:
    @param pixelSize:
    @param resolutionFactor : resolution tuning.. before rotation? after rotation?
    @param cubeSize:
    @param densityNegative:
    @return:    
    """
    from math import floor, sqrt
    from pytom_volume import vol

    if len(atomList) == 0:
        raise RuntimeError('atomList2em : Your atom list is empty!')

    for i in range(len(atomList)):
        atomList[i].setX(round(atomList[i].getX() / pixelSize))
        atomList[i].setY(round(atomList[i].getY() / pixelSize))
        atomList[i].setZ(round(atomList[i].getZ() / pixelSize))

    periodicTableAvailable = True
    try:
        # searching for periodic table library http://pypi.python.org/pypi/periodictable
        from periodictable import elements
    except ImportError:
        periodicTableAvailable = False

    maxValues = [ -10000.0, -10000.0, -10000.0 ]
    minValues = [ 10000.0, 10000.0, 10000.0 ]
    
    for i in range(len(atomList)):
        x = int(atomList[i].getX())
        y = int(atomList[i].getY())
        z = int(atomList[i].getZ())
        currentValues = [x, y, z]

        for i in [0,1,2]:
            if currentValues[i] > maxValues[i]:
                maxValues[i] = currentValues[i]
            if currentValues[i] < minValues[i]:
                minValues[i] = currentValues[i]
    
    if verbose:
        print("---------------------")
        print("maxValues : ", maxValues)
        print("minValues : ", minValues)

    #################### COMPACT CUBE VOLUME ####################
    compactX, compactY, compactZ = maxValues[0]-minValues[0], maxValues[1]-minValues[1], maxValues[2]-minValues[2]
    cubeSize = int(sqrt( compactX**2 + compactY**2 + compactZ**2 ))
    if cubeSize % 2 == 0:
        cubeSize += 1 # Let cubeSize to be odd.
    
    volumeCompact = vol(cubeSize, cubeSize, cubeSize)

    # add 1 is crucial, basically
    volumeCompact.setAll(0.0)
    overlap = 0
    for i in range(len(atomList)):
        x = int(atomList[i].getX() - minValues[0] + 0.5 * ( cubeSize - compactX ))
        y = int(atomList[i].getY() - minValues[1] + 0.5 * ( cubeSize - compactY ))
        z = int(atomList[i].getZ() - minValues[2] + 0.5 * ( cubeSize - compactZ ))

        currentValue = volumeCompact.getV(x, y, z)
        if periodicTableAvailable:
            atomName = atomList[i].getAtomType()[0]
            element = elements.symbol(atomName)
            mass = element.mass
            if currentValue != 0:
                overlap += 1
            volumeCompact.setV(currentValue + mass, x, y, z)
            
        else:
            if atomList[i].getAtomType()[0] == 'H':  ##maybe take this out
                volumeCompact.setV(currentValue + 1.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'C':
                volumeCompact.setV(currentValue + 6.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'N':
                volumeCompact.setV(currentValue + 7.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'O':
                volumeCompact.setV(currentValue + 8.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'P':
                volumeCompact.setV(currentValue + 15.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'S':
                volumeCompact.setV(currentValue + 16.0, x, y, z)

    if densityNegative:
        volumeCompact = volumeCompact * -1

    print(overlap, " : is the overlap counted")
    return volumeCompact, cubeSize/2, cubeSize/2, cubeSize/2

def atomList2emCompact(atomList, pixelSize, densityNegative=False, verbose=False):
    """
    atomList2em:
    @param atomList:
    @param pixelSize:
    @param cubeSize:
    @param densityNegative:
    @return:    
    """
    from math import floor
    from pytom_volume import vol

    if len(atomList) == 0:
        raise RuntimeError('atomList2em : Your atom list is empty!')

    centroidX = 0
    centroidY = 0
    centroidZ = 0

    for i in range(len(atomList)):
        centroidX += atomList[i].getX()
        centroidY += atomList[i].getY()
        centroidZ += atomList[i].getZ()

    # centroidX = centroidX / len(atomList)
    # centroidY = centroidY / len(atomList)
    # centroidZ = centroidZ / len(atomList)

    # centerX = floor(float(cubeSize) / 2.0)
    # centerY = floor(float(cubeSize) / 2.0)
    # centerZ = floor(float(cubeSize) / 2.0)

    # shiftX = centroidX - centerX
    # shiftY = centroidY - centerY
    # shiftZ = centroidZ - centerZ

    for i in range(len(atomList)):
        # atomList[i].setX(round(atomList[i].getX() / pixelSize) + centerX)
        # atomList[i].setY(round(atomList[i].getY() / pixelSize) + centerY)
        # atomList[i].setZ(round(atomList[i].getZ() / pixelSize) + centerZ)
        atomList[i].setX(round(atomList[i].getX() / pixelSize))
        atomList[i].setY(round(atomList[i].getY() / pixelSize))
        atomList[i].setZ(round(atomList[i].getZ() / pixelSize))

    periodicTableAvailable = True
    try:
        # searching for periodic table library http://pypi.python.org/pypi/periodictable
        from periodictable import elements
    except ImportError:
        periodicTableAvailable = False

    maxValues = [ -1000.0, -1000.0, -1000.0 ]
    minValues = [ 1000.0, 1000.0, 1000.0 ]
    
    for i in range(len(atomList)):
        x = int(atomList[i].getX())
        y = int(atomList[i].getY())
        z = int(atomList[i].getZ())
        currentValues = [x, y, z]

        for i in [0,1,2]:
            if currentValues[i] > maxValues[i]:
                maxValues[i] = currentValues[i]
            if currentValues[i] < minValues[i]:
                minValues[i] = currentValues[i]
    
    if verbose:
        print("---------------------")
        print("maxValues : ", maxValues)
        print("minValues : ", minValues)
        print("centroids : ", [centroidX, centroidY, centroidZ])
        print("centers   : ", [centerX, centerY, centerZ])
        print("shifts    : ", [shiftX, shiftY, shiftZ])

    #################### COMPACT VOLUME ####################
    compactX, compactY, compactZ = maxValues[0]-minValues[0], maxValues[1]-minValues[1], maxValues[2]-minValues[2]
    volumeCompact = vol(compactX+1, compactY+1, compactZ+1)
    # add 1 is crucial, basically
    volumeCompact.setAll(0.0)

    for i in range(len(atomList)):
        x = int(atomList[i].getX()) - minValues[0]
        y = int(atomList[i].getY()) - minValues[1]
        z = int(atomList[i].getZ()) - minValues[2]

        currentValue = volumeCompact(x, y, z)
        if periodicTableAvailable:
            atomName = atomList[i].getAtomType()[0]
            element = elements.symbol(atomName)
            mass = element.mass
            volumeCompact.setV(currentValue + mass, x, y, z)
            
        else:
            if atomList[i].getAtomType()[0] == 'H':  ##maybe take this out
                volumeCompact.setV(currentValue + 1.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'C':
                volumeCompact.setV(currentValue + 6.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'N':
                volumeCompact.setV(currentValue + 7.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'O':
                volumeCompact.setV(currentValue + 8.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'P':
                volumeCompact.setV(currentValue + 15.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'S':
                volumeCompact.setV(currentValue + 16.0, x, y, z)

    if densityNegative:
        volumeCompact = volumeCompact * -1

    return volumeCompact, compactX, compactY, compactZ

def atomList2em(atomList, pixelSize, cubeSize, densityNegative=False):
    """
    atomList2em:
    @param atomList:
    @param pixelSize:
    @param cubeSize:
    @param densityNegative:
    @return:    
    """
    from math import floor
    from pytom_volume import vol

    if len(atomList) == 0:
        raise RuntimeError('atomList2em : Your atom list is empty!')

    # get map
    volume = vol(cubeSize, cubeSize, cubeSize)
    volume.setAll(0.0)

    centroidX = 0
    centroidY = 0
    centroidZ = 0

    for i in range(len(atomList)):
        centroidX += atomList[i].getX()
        centroidY += atomList[i].getY()
        centroidZ += atomList[i].getZ()

    centroidX = centroidX / len(atomList)
    centroidY = centroidY / len(atomList)
    centroidZ = centroidZ / len(atomList)

    centerX = floor(float(cubeSize) / 2.0)
    centerY = floor(float(cubeSize) / 2.0)
    centerZ = floor(float(cubeSize) / 2.0)

    shiftX = centroidX - centerX
    shiftY = centroidY - centerY
    shiftZ = centroidZ - centerZ

    for i in range(len(atomList)):
        atomList[i].setX(round(atomList[i].getX() / pixelSize) + centerX)
        atomList[i].setY(round(atomList[i].getY() / pixelSize) + centerY)
        atomList[i].setZ(round(atomList[i].getZ() / pixelSize) + centerZ)

    periodicTableAvailable = True
    try:
        # searching for periodic table library http://pypi.python.org/pypi/periodictable
        from periodictable import elements
    except ImportError:
        periodicTableAvailable = False

    for i in range(len(atomList)):
        x = int(atomList[i].getX())
        y = int(atomList[i].getY())
        z = int(atomList[i].getZ())

        if x >= cubeSize or y >= cubeSize or z >= cubeSize:
            raise RuntimeError('Cube size is too small. Please specify a larger cube for PDB structure!')

        currentValue = volume(x, y, z)

        if periodicTableAvailable:
            atomName = atomList[i].getAtomType()[0]
            element = elements.symbol(atomName)
            mass = element.mass
            volume.setV(currentValue + mass, x, y, z)
        else:
            if atomList[i].getAtomType()[0] == 'H':  ##maybe take this out
                volume.setV(currentValue + 1.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'C':
                volume.setV(currentValue + 6.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'N':
                volume.setV(currentValue + 7.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'O':
                volume.setV(currentValue + 8.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'P':
                volume.setV(currentValue + 15.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'S':
                volume.setV(currentValue + 16.0, x, y, z)


    if densityNegative:
        volume = volume * -1

    return volume

def cifpdb2em(inputPath, pixelSize, cubeSize=0.0, toCompact=False, chain=None, densityNegative=False, fname='', recenter=True):
    """
    cifpdb2em: Creates an volume out of a mmCIF file or pdb file
    @param inputPath: Path to mmCIF file 
    @param pixelSize: The pixel size to convert to 
    @param cubeSize: Resulting cube size
    @param toCompact: Option for compact cuboid
    @return: A volume
    """
    if inputPath.endswith(".pdb"):
        atomList = naivePDBParser(inputPath, chain)
    else:
        atomList = mmCIFParser(inputPath, chain)

    compactX, compactY, compactZ = 0.0, 0.0, 0.0
    if toCompact:
        volm, compactX, compactY, compactZ = atomList2emCube(atomList, pixelSize, densityNegative)
    else:
        volm = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if (not toCompact) and recenter:
        volm = recenterVolume(volm, densityNegative)
    
    if fname:
        volm.write(fname)
        print(f"MRC file is written in {fname}")
    
    return volm, compactX, compactY, compactZ

def getResolution(filePath):
    from pytom.tools.files import checkFileExists

    if not checkFileExists(filePath):
        raise RuntimeError('resolutionResize : input File not found! ', filePath)
    
    if filePath.endswith(".pdb"):
        resPatternPDB = re.compile("RESOLUTION\..*([0-9]+\.[0-9]+).*ANGSTROMS")
        f = open(filePath, 'r')
        pdbContent = f.read()
        f.close()
        return float(re.findall(resPatternPDB, pdbContent)[0])
    elif filePath.endswith(".cif"):
        resPatternCIF = re.compile("_em_3d_reconstruction.resolution +([0-9]+\.[0-9]+)")
        f = open(filePath, 'r')
        cifContent = f.read()
        f.close()
        return float(re.findall(resPatternCIF, cifContent)[0])
        print(f"ends with cif resolution is {resolution}")
    else:
        print("Unsupported File Extension")
        raise RuntimeError('Unsupported file extenstion : ', filePath)

def volume2MRC(volPath, mrcPath, floatMRC=False, overwrite=False, verbose=False):
    inputVolume = read(volPath)
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    if verbose:
        print(f"Volume dimension is initially... {x}x{y}x{z}")
    
    if floatMRC:
        volumeData = np.empty([x, y, z], dtype = np.float32)
    else:
        volumeData = np.empty([x, y, z], dtype = np.int8)
    
    for i in range(inputVolume.sizeX()):
        for j in range(inputVolume.sizeY()):
            for k in range(inputVolume.sizeZ()):
                #print(type(inputVolume.getV(i,j,k)))
                volumeData[i,j,k] = inputVolume.getV(i,j,k)
    
    with mrcfile.new(mrcPath, overwrite=overwrite) as mrc:
        mrc.set_data(volumeData)
        print(f"mrc data dimension is converted to... {mrc.data.shape}")
    return

def wgetByPDBID(pdbID, pdbDir):
    pdbPath = f"{pdbDir}/{pdbID}.pdb"
    cifPath = f"{pdbDir}/{pdbID}.cif"
    pdbURL = f"https://files.rcsb.org/view/{pdbID}.pdb"
    cifURL = f"https://files.rcsb.org/view/{pdbID}.cif"

    if os.path.isfile(pdbPath):
        return pdbPath
    if os.path.isfile(cifPath):
        return cifPath

    isPDB = True
    URL = ""
    Path = ""
    response = requests.get(pdbURL)
    if not response.status_code == 200:
        response = requests.get(cifURL)
        if not response.status_code == 200:
            raise RuntimeError("Invalid pdb ID maybe, ", pdbID)
        isPDB = False
        URL = cifURL
        Path = cifPath
    else:
        URL = pdbURL
        Path = pdbPath
    wget.download(URL, out=Path)
    return Path

def wgetPDB2Volume(pdbID, pdbDir, volumeDir, mrcDir, cubeSize=0.0, toCompact=False, generateMRC=False, overwrite=False, verbose=False):
    """
    wgetPDB2Volume : Creates an PDB(CIF) file, EM file, MRC file from a PDB ID.
    @param overwrite : is for overwrite mrcfile(Volume2MRC).
    """
    # pdbDir should not include dangling /
    # densityNegative for default
    volumePath = f"{volumeDir}/{pdbID}.em"
    mrcPath = f"{mrcDir}/{pdbID}.mrc"
    
    if verbose:
        print(f"wgetPDB2Volume is working with PDBID : {pdbID}")
    
    Path = wgetByPDBID(pdbID, pdbDir)
    resolution = getResolution(Path)
    _vol, compactX, compactY, compactZ = cifpdb2em(Path, pixelSize=resolution, cubeSize=cubeSize, toCompact=toCompact, chain=None, fname=volumePath, densityNegative=False, recenter=True)
    if generateMRC:
        volume2MRC(volumePath, mrcPath, overwrite=overwrite)
    return _vol, compactX, compactY, compactZ, resolution

def prepareCubeVolumes(pdbIDList, pdbDir, volumeDir, mrcDir, cubeSize=0.0, toCompact=False, generateMRC=False, overwrite=False, verbose=False):
    createdVolumes = []
    for pdbID in pdbIDList:
        vol, x, y, z, _r = wgetPDB2Volume(pdbID, pdbDir, volumeDir, mrcDir, cubeSize=cubeSize, toCompact=toCompact, generateMRC=generateMRC, overwrite=overwrite, verbose=verbose)
        createdVolumes.append(vol)
    return createdVolumes

def makeCompact(inputVolume):
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    min = [x+1, y+1, z+1]
    max = [-1, -1, -1]
    for i in range(x):
        for j in range(y):
            for k in range(z):
                if(inputVolume.getV(i,j,k) != 0):
                    if i > max[0]:
                        max[0] = i
                    elif i < min[0]:
                        min[0] = i
                    if j > max[1]:
                        max[1] = j
                    elif j < min[1]:
                        min[1] = j
                    if k > max[2]:
                        max[2] = k
                    elif k < min[2]:
                        min[2] = k
    dif = [max[0]-min[0]+1, max[1]-min[1]+1, max[2]-min[2]+1]
    compactVol = vol(dif[0], dif[1], dif[2])
    compactVol.setAll(0.0)
    for i in range(dif[0]):
        for j in range(dif[1]):
            for k in range(dif[2]):
                compactVol.setV( inputVolume.getV(i+min[0], j+min[1], k+min[2]) , i, j, k)
    return compactVol
##########################################################################################################################################################################
#  Section for Multi particle scenario.
##########################################################################################################################################################################
def makeScenarioByPDBIDs(pdbIDList, volumeDir, scenarioDir, scenarioIdentifier="noname", cubeSize=128, pfailedAttempts=9000, pparticleNum=1600, overwrite=False, isRotation=False, verbose=False):
    # cuboidalOccupancyList = [['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]]]
    startTime = time.time()
    scenarioMetaDataFile = f"{scenarioDir}/{scenarioIdentifier}.txt"
    scenarioVolumeFile = f"{scenarioDir}/{scenarioIdentifier}.em"

    f = open(scenarioMetaDataFile, 'w')
    fulltomX, fulltomY, fulltomZ = 2*cubeSize, 2*cubeSize, 2*cubeSize
    volume = vol(fulltomX, fulltomY, fulltomZ)
    volume.setAll(0.0)
    
    failedAttempts = 0
    particleNum = 0
    scenario = []

    classNum = np.random.randint(low=0, high=len(pdbIDList)) #TODO NUMPY.
    currentTemplate = f"{volumeDir}/{pdbIDList[classNum]}.em"
    currentVol = read(currentTemplate)
    if isRotation == True:
        rotatedVol, phi, theta, psi = compactRandomRotation(currentVol)
    else:
        rotatedVol, phi, theta, psi = currentVol, 0, 0, 0

    sizeX, sizeY, sizeZ = rotatedVol.sizeX(), rotatedVol.sizeY(), rotatedVol.sizeZ()
    x, y, z = np.random.randint(low=0, high=fulltomX-sizeX), np.random.randint(low=0, high=fulltomY-sizeY), np.random.randint(low=0, high=fulltomZ-sizeZ)
    centerX, centerY, centerZ = x+ sizeX/2, y+ sizeY/2, z+ sizeZ/2
    
    scenario.append([[x,y,z], [x+sizeX-1, y+sizeY-1, z+sizeZ-1]])
    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                volume.setV( rotatedVol.getV(i,j,k), x+i, y+j, z+k)
    # INCORPORTATE
    f.write(f"{pdbIDList[classNum]},{centerX},{centerY},{centerZ},{phi},{theta},{psi}\n")
    particleNum+=1
    rotatedVol = None
    rotFailNum = 0
    while failedAttempts < pfailedAttempts and particleNum < pparticleNum:
        if verbose and failedAttempts%300 == 0 and failedAttempts != 0:
            print(f"... Now failed Attemps are {failedAttempts}")   
        if isRotation == True:
            if rotFailNum == 20:
                rotatedVol = None
            if rotatedVol == None:
                rotFailNum = 0
                classNum = np.random.randint(low=0, high=len(pdbIDList)) #TODO NUMPY.
                currentTemplate = f"{volumeDir}/{pdbIDList[classNum]}.em"
                currentVol = read(currentTemplate)
                rotatedVol, phi, theta, psi = compactRandomRotation(currentVol)
            else:
                rotFailNum+=1
        else:
            classNum = np.random.randint(low=0, high=len(pdbIDList)) #TODO NUMPY.
            currentTemplate = f"{volumeDir}/{pdbIDList[classNum]}.em"
            currentVol = read(currentTemplate)
            rotatedVol, phi, theta, psi = currentVol, 0, 0, 0

        sizeX, sizeY, sizeZ = rotatedVol.sizeX(), rotatedVol.sizeY(), rotatedVol.sizeZ()
        # x, y, z = random.randint(0, fulltomX-1-sizeX), random.randint(0, fulltomY-1-sizeY), random.randint(0, fulltomZ-1-sizeZ)
        x, y, z = np.random.randint(low=0, high=fulltomX-sizeX), np.random.randint(low=0, high=fulltomY-sizeY), np.random.randint(low=0, high=fulltomZ-sizeZ)
        centerX, centerY, centerZ = x+ sizeX/2, y+ sizeY/2, z+ sizeZ/2
        
        isOccupied = False
        for existingItem in scenario:
            if (x <= existingItem[1][0] and x+sizeX >= existingItem[0][0]) and \
                (y <= existingItem[1][1] and y+sizeY >= existingItem[0][1]) and \
                (z <= existingItem[1][2] and z+sizeZ >= existingItem[0][2]):
                failedAttempts+=1
                isOccupied = True
                break; 
        if isOccupied == False:
            scenario.append([[x, y, z], [x+sizeX-1, y+sizeY-1, z+sizeZ-1]])
            for i in range(sizeX):
                for j in range(sizeY):
                    for k in range(sizeZ):
                        curVal = rotatedVol.getV(i,j,k)
                        if curVal != 0:
                            volume.setV( curVal, x+i, y+j, z+k)
            particleNum+=1
            rotatedVol = None
            f.write(f"{pdbIDList[classNum]},{centerX},{centerY},{centerZ},{phi},{theta},{psi}\n")
            if particleNum % 10 == 0:
                print(f"... Particle Num : {particleNum}")  
    f.close()

    # --- META DATA ---
    appendMetaDataln(f"makeScenarioByPDBIDs {scenarioIdentifier} done - time elapsed : {time.time() - startTime}s")
    appendMetaDataln(f"-output file : {scenarioVolumeFile}, cubeSize : {fulltomX}x{fulltomY}x{fulltomZ}")
    appendMetaDataln(f"-with pdbIDList : {pdbIDList}")
    appendMetaDataln(f"-failedAttempts : {pfailedAttempts}, resultParticleNum : {particleNum}")

    print(f"----------- Scenario generation is done... with Particle Number {particleNum}----------")
    # -----------------
    volume.write(scenarioVolumeFile)

##########################################################################################################################################################################
#  Section for Simulation.
def compactRandomRotation(inputVolume):
    global called
    called+=1
    phi, theta, psi = np.random.randint(low=0, high=360, size=(3,)) # High exclusive
    from pytom_volume import rotate
    rotatedVolume = vol(inputVolume.sizeX(),inputVolume.sizeY(),inputVolume.sizeZ())
    rotatedVolume.setAll(0.0)
    rotate(inputVolume, rotatedVolume, int(phi), int(theta), int(psi))
    comp = makeCompact(rotatedVolume)
    return comp, phi, theta, psi

def customSimulation(volumePath, simulatedPath=None, snrValue=0.1, rotation=None, wedgeAngle=None, shift=None):
    # Rotation : [ x axis , z axis , y axis ]
    wedge = 0.0
    shiftList = [0, 0, 0]

    v = read(volumePath)
    if rotation == None:
        rotation = [0, 0, 0]
    if not wedgeAngle == None:
        wedge = wedgeAngle
    if not shift == None:
        shiftList = shift
    
    wi = WedgeInfo(wedgeAngle=wedge, cutoffRadius=0.0)
    s = simpleSimulation( volume=v, rotation=rotation, shiftV=shiftList, wedgeInfo=wi, SNR=snrValue)
    if simulatedPath:
        s.write(simulatedPath)
    
    appendMetaDataln(f"customSimulation is done : inputfile is {volumePath} / outputfile is {simulatedPath}")
    appendMetaDataln(f"-snr : {snrValue}, rotation : {rotation}, wedgeAngle : {wedgeAngle}, shift : {shift}")
    return s

##########################################################################################################################################################################
#  Section for Utility.

## Write 2D matrix into the volume 
def matrixToVolLayerZ(vol, matrix, zheight):
    for i in range(vol.sizeX()):
        for j in range(vol.sizeY()):
            vol.setV(matrix[i][j], i, j, zheight)
    return vol
def matrixToVolLayerY(vol, matrix, yheight):
    for i in range(vol.sizeX()):
        for j in range(vol.sizeZ()):
            vol.setV(matrix[i][j], i, j, yheight)
    return vol
def matrixToVolLayerX(vol, matrix, xheight):
    for i in range(vol.sizeZ()):
        for j in range(vol.sizeY()):
            vol.setV(matrix[i][j], i, j, xheight)
    return vol

def compactCuboid2rotateCube(cuboidPath, cubePath):
    cuboid = read(cuboidPath)
    x, y, z = cuboid.sizeX(), cuboid.sizeY(), cuboid.sizeZ()
    size = int(math.sqrt( x**2 + y**2 + z**2 )) + 1 # 0.5 * 2 = 1.
    cube = vol(size, size, size)
    
    sx, sy, sz = int((size-x)/2), int((size-y)/2), int((size-z)/2) 
    for i in range(x):
        for j in range(y):
            for k in range(z):
                cube.setV( cuboid.getV(i, j, k), i+sx, j+sy, k+sz )
    cube.write(cubePath)

def volumeOutliner(inputVolumes, outlineValue=10):
    for inputVolume in inputVolumes:
        x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
        for i in range(x):
            inputVolume.setV(outlineValue, i, 0,0)
            inputVolume.setV(outlineValue, i, y-1,0)
            inputVolume.setV(outlineValue, i, 0,z-1)
            inputVolume.setV(outlineValue, i, y-1,z-1)
        for i in range(y):
            inputVolume.setV(outlineValue, 0, i,0)
            inputVolume.setV(outlineValue, x-1, i,0)
            inputVolume.setV(outlineValue, 0, i,z-1)
            inputVolume.setV(outlineValue, x-1, i,z-1)
        for i in range(z):
            inputVolume.setV(outlineValue, 0, 0,i)
            inputVolume.setV(outlineValue, x-1, 0,i)
            inputVolume.setV(outlineValue, 0, y-1,i)
            inputVolume.setV(outlineValue, x-1, y-1,i)

def volumeListWriter(inputVolumes, outputDir, Description):
    index = 0
    for inputVolume in inputVolumes:
        inputVolume.write(f"{outputDir}/{Description}_{index}.em")
        index += 1
##########################################################################################################################################################################
SHREC2021_FULL = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2", "5mrc" ]
SHREC2021_FULLex5mrc = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2" ]

SHREC2021_FULL_OCCLIST = [['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]
SHREC2021_FULLexc2L2S = [ "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr" ]
SHREC2021_1bxn = [ "1bxn" ]

# FROM PDB ID -> Resolution corrected Compact Cuboid!
def resolutionResizeUpper(identifier, pdbDir, volumeDir, outputDir, toResolution):
    inputPDBPath = f"{pdbDir}/{identifier}.pdb"
    from pytom.tools.files import checkFileExists

    if not checkFileExists(inputPDBPath):
        inputPDBPath = f"{pdbDir}/{identifier}.cif"
        if not checkFileExists(inputPDBPath):
            raise RuntimeError('resolutionResize : input File not found! ', filePath)

    inputVolumePath = f"{volumeDir}/{identifier}.em"
    outputVolumePath = f"{outputDir}/{identifier}.em"
    # Assume input is cube form!!
    resolution = getResolution( inputPDBPath )
    inputVolume = read( inputVolumePath )
    X, Y, Z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    resizedSize = math.ceil( (inputVolume.sizeX() - 1)*resolution / toResolution ) 
    outputVolume = vol(resizedSize, Y, Z)
    outputVolume.setAll(0.0)
    print(f" BEFORE SIZE : {inputVolume.sizeX()} x {inputVolume.sizeY()} x {inputVolume.sizeZ()}")
    print(f" AFTER SIZE : {outputVolume.sizeX()} x {outputVolume.sizeX()} x {outputVolume.sizeX()}")
    # interpolate.
    for i in range(resizedSize):
        realcoord_im1 = (i-1) * toResolution
        realcoord_i   =   i   * toResolution
        realcoord_ip1 = (i+1) * toResolution

        lowerIdx = math.ceil(realcoord_im1 / resolution) if i-1 > 0 else 0
        upperIdx = math.floor(realcoord_ip1 / resolution)
        for j in range(Y):
            for k in range(Z):
                curVal = 0.0
                idx = lowerIdx
                #print(lowerIdx, "~", upperIdx, "-----------------------------------------------------------------------")
                while idx <= upperIdx:
                    realcoord_org = round( idx * resolution, 4 )
                    try:
                        modfactor = round( ( toResolution - abs( realcoord_i - realcoord_org ) ) / toResolution, 4 )
                        #print(toResolution, realcoord_i, realcoord_org, "at ", idx)
                        #print(modfactor, "at ", idx)
                        curVal += inputVolume.getV(idx, j, k) * modfactor
                    except:
                        pass
                    idx += 1
                outputVolume.setV(curVal, i, j, k)

    inputVolume = outputVolume
    #print("inputVolume", inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ())
    outputVolume = vol(resizedSize, resizedSize, Z)
    outputVolume.setAll(0.0)
    for j in range(resizedSize):
        realcoord_jm1 = (j-1) * toResolution
        realcoord_j   =   j   * toResolution
        realcoord_jp1 = (j+1) * toResolution

        lowerIdx = math.ceil(realcoord_jm1 / resolution) if j-1 > 0 else 0
        upperIdx = math.floor(realcoord_jp1 / resolution)
        for i in range(resizedSize):
            for k in range(Z):
                curVal = 0.0
                idx = lowerIdx
                while idx <= upperIdx:
                    realcoord_org = idx * resolution
                    try:
                        curVal += inputVolume.getV(i, idx, k) * ( toResolution - abs( realcoord_i - realcoord_org ) ) / toResolution
                    except:
                        pass
                        #print(idx, j, k, "//", X, Y, Z, "//", lowerIdx, "~", upperIdx)
                    idx += 1
                outputVolume.setV(curVal, i, j, k)

    inputVolume = outputVolume
    outputVolume = vol(resizedSize, resizedSize, resizedSize)
    outputVolume.setAll(0.0)
    for k in range(resizedSize):
        realcoord_km1 = (k-1) * toResolution
        realcoord_k   =   k   * toResolution
        realcoord_kp1 = (k+1) * toResolution
        lowerIdx = math.ceil(realcoord_km1 / resolution) if k-1 > 0 else 0
        upperIdx = math.floor(realcoord_kp1 / resolution)
        for i in range(resizedSize):
            for j in range(resizedSize):
                curVal = 0.0
                idx = lowerIdx
                while idx <= upperIdx:
                    realcoord_org = idx * resolution
                    try:
                        curVal += inputVolume.getV(i, j, idx) * ( toResolution - abs( realcoord_i - realcoord_org ) ) / toResolution
                    except:
                        pass
                    idx += 1
                outputVolume.setV(curVal, i, j, k)
    outputVolume.write(outputVolumePath)

def resolutionResizeUnity(identifier, pdbDir, volumeDir, outputDir, toResolution):
    inputPDBPath = f"{pdbDir}/{identifier}.pdb"
    from pytom.tools.files import checkFileExists

    if not checkFileExists(inputPDBPath):
        inputPDBPath = f"{pdbDir}/{identifier}.cif"
        if not checkFileExists(inputPDBPath):
            raise RuntimeError('resolutionResize : input File not found! ', filePath)

    inputVolumePath = f"{volumeDir}/{identifier}.em"
    outputVolumePath = f"{outputDir}/{identifier}.em"
    # Assume input is cube form!!
    resolution = getResolution( inputPDBPath )
    inputVolume = read( inputVolumePath )
    X, Y, Z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    resizedSize = math.ceil( (inputVolume.sizeX() - 1)*resolution / toResolution ) 
    outputVolume = vol(resizedSize, resizedSize, resizedSize)
    outputVolume.setAll(0.0)
    print(f" BEFORE SIZE : {inputVolume.sizeX()} x {inputVolume.sizeY()} x {inputVolume.sizeZ()}")
    print(f" AFTER SIZE : {outputVolume.sizeX()} x {outputVolume.sizeX()} x {outputVolume.sizeX()}")
    # interpolate.
    for i in range(resizedSize):
        for j in range(resizedSize):
            for k in range(resizedSize):
                curVal = 0.0

                realcoord_i   =   i   * toResolution
                realcoord_j   =   j   * toResolution
                realcoord_k   =   k   * toResolution
                
                lowerIdxX = math.ceil((i-1) * toResolution / resolution) if i-1 > 0 else 0
                upperIdxX = math.floor((i+1) * toResolution / resolution)
                lowerIdxY = math.ceil((j-1) * toResolution / resolution) if j-1 > 0 else 0
                upperIdxY = math.floor((j+1) * toResolution / resolution)
                lowerIdxZ = math.ceil((k-1) * toResolution / resolution) if k-1 > 0 else 0
                upperIdxZ = math.floor((k+1) * toResolution / resolution)

                idxX = lowerIdxX
                while idxX <= upperIdxX:
                    realcordX = round( idxX * resolution, 4)
                    idxY = lowerIdxY
                    while idxY <= upperIdxY:
                        realcordY = round( idxY * resolution, 4)
                        idxZ = lowerIdxZ
                        while idxZ <= upperIdxZ:
                            realcordZ = round( idxZ * resolution, 4 )
                            modfactorX = round( ( toResolution - abs( realcoord_i - realcordX ) ) / toResolution, 4 )
                            modfactorY = round( ( toResolution - abs( realcoord_j - realcordY ) ) / toResolution, 4 )
                            modfactorZ = round( ( toResolution - abs( realcoord_k - realcordZ ) ) / toResolution, 4 )

                            modfactor = modfactorX * modfactorY * modfactorZ
                            try:
                                val = inputVolume.getV(idxX, idxY, idxZ) * modfactor
                                curVal += val
                            except:
                                pass
                            idxZ += 1
                        idxY += 1
                    idxX += 1
                outputVolume.setV(curVal, i, j, k)
    outputVolume.write(outputVolumePath)
# FROM PDB ID -> Resolution corrected Compact Cuboid!
def wgetPDB2VolumeR(pdbID, pdbDir, volumeDir, mrcDir, overwrite=False, verbose=False):
    # toCompact = True.
    # Always generate MRC. (generateMRC = True)
    _vol, _x, _y, _z, resolution = wgetPDB2Volume(pdbID, pdbDir, volumeDir, mrcDir, toCompact=True, generateMRC=True, overwrite=overwrite, verbose=verbose)
    
    
    #volume -> 

if __name__ == "__main__":
    executionStart = time.time()
    #################### Workspace ##################    
    now = time.localtime()
    programTime = f"{now.tm_year}/{now.tm_mon}/{now.tm_mday} {now.tm_hour}:{now.tm_min}:{now.tm_sec}"
    appendMetaDataln(f"===> Scripts running : {programTime}")
    # Put some description.
    DESCRIPTION = "5_0_outline test"
    appendMetaDataln(f"===> {DESCRIPTION}")

    # 01 Work : Just test tomogram without rotation. #######################################################################
    # makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="1_withoutrotation", cubeSize=256, pfailedAttempts=4000, isRotation=False, verbose=True)
    # makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="2_withoutrotation", cubeSize=256, pfailedAttempts=4000, isRotation=False, verbose=True)
    #volume2MRC("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1_withoutrotation_vol.mrc")
    #volume2MRC("/cdata/scenario/2_withoutrotation.em", "/cdata/scenario/2_withoutrotation_vol.mrc")
    # customSimulation("/cdata/scenario/1_withoutrotation.em", simulatedPath="/cdata/simulated/1_withoutrotation.em", snrValue=2.0)
    # customSimulation("/cdata/scenario/2_withoutrotation.em", simulatedPath="/cdata/simulated/2_withoutrotation.em", snrValue=2.0)
    # volume2MRC("/cdata/simulated/1_withoutrotation.em", "/cdata/simulated/1_withoutrotation.mrc")
    # volume2MRC("/cdata/simulated/2_withoutrotation.em", "/cdata/simulated/2_withoutrotation.mrc")
    ########################################################################################################################
    # 02 Get resolution information.
    # prepareCubeVolumes(SHREC2021_FULL, pdbDir=pdbDataDIR, volumeDir=singleParticleEMCubeDIR, mrcDir=None, toCompact=True, overwrite=True, verbose=True)
    # 03 Improve : Execute Rotation only None. 
    # BEFORE :: makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="3_before improve", cubeSize=256, pfailedAttempts=8000, isRotation=True, verbose=True)
    #makeScenarioByPDBIDs(SHREC2021_FULLex5mrc, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="4_6_rotFailMax", cubeSize=256, pfailedAttempts=400000, isRotation=True, verbose=True)
    #makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="4_7_rotFailMaxFull", cubeSize=256, pfailedAttempts=400000, isRotation=True, verbose=True)
    #makeScenarioByPDBIDs(SHREC2021_FULLex5mrc, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="4_6_rotFailMax2", cubeSize=256, pfailedAttempts=400000, isRotation=True, verbose=True)
    #makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="4_7_rotFailMaxFull2", cubeSize=256, pfailedAttempts=400000, isRotation=True, verbose=True)
    #makeScenarioByPDBIDs(SHREC2021_FULLex5mrc, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="4_6_rotFailMax", cubeSize=256, pfailedAttempts=400000, isRotation=True, verbose=True)
    #makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="4_8_after improve(720)", cubeSize=256, pfailedAttempts=400000, pparticleNum=720, isRotation=True, verbose=True)
    #time1 = time.time()
    #volume2MRC("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1worot_timetestMine.mrc")
    #time2 = time.time()
    #em2mrc("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1worot_timetestPYTOM.mrc")
    #time3 = time.time()
    #volume2MRC("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1_withoutrotation_float16.mrc")

    #volList = prepareCubeVolumes(SHREC2021_FULL, pdbDataDIR, singleParticleEMCubeDIR, singleParticleMRCCubeDIR, cubeSize=0.0, toCompact=True, generateMRC=False, overwrite=False, verbose=False)
    #volumeOutliner(volList, outlineValue = 40)
    #volumeListWriter(volList, "/cdata/outlined", "0313_beforeRes")

    #resolutionResizeUnity("1bxn", "/cdata/pdbData", "/cdata/singleParticleEM_cube", "/cdata/resolution2", 10.0)
    # volList = [ read("/cdata/resolution/1bxn.em") ]
    # volumeOutliner(volList, outlineValue = 30)
    # volumeListWriter(volList, "/cdata/outlined", "0313_RETEST")


    for pdbID in SHREC2021_FULL:
        resolutionResizeUnity(pdbID, "/cdata/pdbData", "/cdata/singleParticleEM_cube", "/cdata/resolution2", 10.0)

    volList = []
    for pdbID in SHREC2021_FULL:
        volList.append(read(f"/cdata/resolution2/{pdbID}.em"))
        
    volumeOutliner(volList, outlineValue = 300)
    volumeListWriter(volList, "/cdata/outlined", "0313_UNITY")


    #resolutionResize("/cdata/pdbData/1bxn.pdb", "/cdata/singleParticleEM_cube/1bxn.em", 10.0, "/cdata/resolution/1bxn.em")
    #resolutionResize("/cdata/pdbData/5mrc.cif", "/cdata/singleParticleEM_cube/5mrc.em", 10.0, "/cdata/resolution/5mrc.em")

    #print(f" 1 : {time2 - time1}")
    #print(f" 2 : {time3 - time2}")
    #appendMetaDataln(f"called num of compact rotation {called}")
    ################### Workspace Ended #############
    print(f" All of the Jobs completed with elapsed time : {time.time()-executionStart}")
    #################### Program Ended ##############