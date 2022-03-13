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
def atomList2emCube(atomList, pixelSize, densityNegative=False, resolutionFactor=None, verbose=False):
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

    for i in range(len(atomList)):
        x = int(atomList[i].getX() - minValues[0] + 0.5 * ( cubeSize - compactX ))
        y = int(atomList[i].getY() - minValues[1] + 0.5 * ( cubeSize - compactY ))
        z = int(atomList[i].getZ() - minValues[2] + 0.5 * ( cubeSize - compactZ ))

        currentValue = volumeCompact.getV(x, y, z)
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

def cifpdb2em(inputPath, isPDB, pixelSize, cubeSize=0.0, toCompact=False, chain=None, densityNegative=False, fname='', recenter=True):
    """
    mmCIF2em: Creates an volume out of a mmCIF file
    @param mmCIFPath: Path to mmCIF file 
    @param pixelSize: The pixel size to convert to 
    @param cubeSize: Resulting cube size
    @param toCompact: Option for compact cuboid
    @return: A volume
    """
    from math import floor
    from pytom_volume import vol

    if isPDB:
        atomList = mmCIFParser(inputPath, chain)
    else:
        atomList = naivePDBParser(inputPath, chain)
    vol = None
    compactX, compactY, compactZ = 0.0, 0.0, 0.0
    if toCompact:
        vol, compactX, compactY, compactZ = atomList2emCube(atomList, pixelSize, densityNegative)
    else:
        vol = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if (not toCompact) and recenter:
        vol = recenterVolume(vol, densityNegative)
    
    if fname:
        vol.write(fname)
        print(f"MRC file is written in {fname}")
    
    return vol, compactX, compactY, compactZ

def getResolution(filePath):
    from pytom.tools.files import checkFileExists

    if not checkFileExists(filePath):
        raise RuntimeError('resolutionResize : input File not found! ', filePath)
    
    if filePath.endswith(".pdb"):
        resPatternPDB = re.compile("RESOLUTION\..*([0-9]+\.[0-9]+).*ANGSTROMS")
        f = open(filePath, 'r')
        pdbContent = f.read()
        f.close()
        return re.findall(resPatternPDB, pdbContent)[0]
    elif filePath.endswith(".cif"):
        resPatternCIF = re.compile("_em_3d_reconstruction.resolution +([0-9]+\.[0-9]+)")
        f = open(filePath, 'r')
        cifContent = f.read()
        f.close()
        return re.findall(resPatternCIF, cifContent)[0]
        print(f"ends with cif resolution is {resolution}")
    else:
        print("Unsupported File Extension")
        raise RuntimeError('Unsupported file extenstion : ', filePath)

def volume2MRCconverter(volPath, mrcPath, floatMRC=False, overwrite=False, verbose=False):
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

def wgetPDB2Volume(pdbID, pdbDir, volumeDir, mrcDir, cubeSize=0.0, toCompact=False, generateMRC=False, overwrite=False, verbose=False):
    # pdbDir should not include dangling /
    # densityNegative for default
    if verbose:
        print(f"wgetPDB2Volume is working with PDBID : {pdbID}")
    
    pdbPath = f"{pdbDir}/{pdbID}.pdb"
    cifPath = f"{pdbDir}/{pdbID}.cif"
    pdbURL = f"https://files.rcsb.org/view/{pdbID}.pdb"
    cifURL = f"https://files.rcsb.org/view/{pdbID}.cif"
    volumePath = f"{volumeDir}/{pdbID}.em"
    mrcPath = f"{mrcDir}/{pdbID}.mrc"

    templateExists = False
    if os.path.isfile(pdbPath) or os.path.isfile(cifPath):
        templateExists = True

    isPDB = True
    URL = "", Path = ""
    response = requests.get(pdbURL)
    if not response.status_code == 200:
        response = requests.get(cifURL)
        if not response.status_code == 200:
            print("Invalid pdb ID maybe.")
            raise RuntimeError("Invalid pdb ID maybe, " pdbID)
        
        isPDB = False
        URL = cifURL
        Path = cifPath
    else:
        URL = pdbURL
        Path = pdbPath
    if not templateExists:
        wget.download(URL, out=Path)

    resolution = getResolution(Path)
    print("Resolution of {} is {}".format(pdbID, resolution))
    _vol, compactX, compactY, compactZ = cifpdb2em(cifPath, isPDB=isPDB, pixelSize=float(resolution), cubeSize=cubeSize, toCompact=toCompact, chain=None, fname=volumePath, densityNegative=False, recenter=True)
    if generateMRC:
        volume2MRCconverter(volumePath, mrcPath, overwrite=overwrite)
    return _vol, compactX, compactY, compactZ

def prepareCubeVolumes(pdbIDList, pdbDir, volumeDir, mrcDir, cubeSize=0.0, toCompact=False, generateMRC=False, overwrite=False, verbose=False):
    createdVolumes = []
    for pdbID in pdbIDList:
        vol, x, y, z = wgetPDB2Volume(pdbID, pdbDir, volumeDir, mrcDir, cubeSize=cubeSize, toCompact=toCompact, generateMRC=generateMRC, overwrite=overwrite, verbose=verbose)
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

def volumeOutliner(inputVolume, outlineValue=10):
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
    return inputVolume
##########################################################################################################################################################################
SHREC2021_FULL = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2", "5mrc" ]
SHREC2021_FULLex5mrc = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2" ]

SHREC2021_FULL_OCCLIST = [['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]
SHREC2021_FULLexc2L2S = [ "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr" ]
SHREC2021_1bxn = [ "1bxn" ]

def em2mrc(filename,newfilename):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists,checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('EM file not found! ',filename)

    emfile = read(filename)
    emfile.write(newfilename,'mrc')

# FROM PDB ID -> Resolution corrected Compact Cuboid!
def resolutionResize(inputPDBPath, inputVolumePath, toResolution, outputPath):
    resolution = getResolution( inputPDBPath )
    
    # MRC convert.
    # Interpolate.
    # Save?
    

if __name__ == "__main__":
    executionStart = time.time()
    #################### Workspace ##################    
    now = time.localtime()
    programTime = f"{now.tm_year}/{now.tm_mon}/{now.tm_mday} {now.tm_hour}:{now.tm_min}:{now.tm_sec}"
    appendMetaDataln(f"===> Scripts running : {programTime}")
    # Put some description.
    DESCRIPTION = "4_8_after improve(720)"
    appendMetaDataln(f"===> {DESCRIPTION}")

    # 01 Work : Just test tomogram without rotation. #######################################################################
    # makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="1_withoutrotation", cubeSize=256, pfailedAttempts=4000, isRotation=False, verbose=True)
    # makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="2_withoutrotation", cubeSize=256, pfailedAttempts=4000, isRotation=False, verbose=True)
    #volume2MRCconverter("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1_withoutrotation_vol.mrc")
    #volume2MRCconverter("/cdata/scenario/2_withoutrotation.em", "/cdata/scenario/2_withoutrotation_vol.mrc")
    # customSimulation("/cdata/scenario/1_withoutrotation.em", simulatedPath="/cdata/simulated/1_withoutrotation.em", snrValue=2.0)
    # customSimulation("/cdata/scenario/2_withoutrotation.em", simulatedPath="/cdata/simulated/2_withoutrotation.em", snrValue=2.0)
    # volume2MRCconverter("/cdata/simulated/1_withoutrotation.em", "/cdata/simulated/1_withoutrotation.mrc")
    # volume2MRCconverter("/cdata/simulated/2_withoutrotation.em", "/cdata/simulated/2_withoutrotation.mrc")
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
    #volume2MRCconverter("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1worot_timetestMine.mrc")
    #time2 = time.time()
    #em2mrc("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1worot_timetestPYTOM.mrc")
    #time3 = time.time()
    volume2MRCconverter("/cdata/scenario/1_withoutrotation.em", "/cdata/scenario/1_withoutrotation_float16.mrc")

    #resolutionResize("/cdata/pdbData/1bxn.pdb", "/cdata/singleParticleEM_cube/1bxn.em", 10.0, "/cdata/resolution/1bxn.em")
    #resolutionResize("/cdata/pdbData/5mrc.cif", "/cdata/singleParticleEM_cube/5mrc.em", 10.0, "/cdata/resolution/5mrc.em")

    #print(f" 1 : {time2 - time1}")
    #print(f" 2 : {time3 - time2}")
    #appendMetaDataln(f"called num of compact rotation {called}")
    ################### Workspace Ended #############
    print(f" All of the Jobs completed with elapsed time : {time.time()-executionStart}")
    #################### Program Ended ##############