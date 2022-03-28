from pytom_volume import vol, initSphere
from pytom.basic.structures import WedgeInfo

# VKJY
from urllib.error import HTTPError
import numpy as np
import mrcfile
import os.path, wget, re, requests, math, random, time, json

# My Other files
from utils import volumeOutliner, volumeListWriter, volObj2Numpy, newNumpyByXYZ, volume2MRC
from pytomLib import recenterVolume, naivePDBParser, mmCIFParser, read, noiseApplier
# from pytomLib import recenterVolume, naivePDBParser, mmCIFParser, vol, initSphere, WedgeInfo

rootDIR = "/cdata"
pdbDataDIR = f"{rootDIR}/pdbData"
singleParticleEMCuboidDIR = f"{rootDIR}/singleParticleEM_cuboid"
singleParticleEMCubeDIR = f"{rootDIR}/singleParticleEM_cube"
singleParticleMRCCuboidDIR = f"{rootDIR}/singleParticleMRC_cuboid"
singleParticleMRCCubeDIR = f"{rootDIR}/singleParticleMRC_cube"
scenarioDIR = f"{rootDIR}/scenario"
metadataFILE = f"{rootDIR}/metadata.txt"

def appendMetaDataln(metadata):
    fmeta = open(metadataFILE, "a")
    fmeta.write(metadata + "\n")
    fmeta.close()

def jsonWriterV5(Dict, Write):
    pass
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
    # atoms = []
    # for atom in atomList:
    #     if not atom.getAtomType() in atoms:
    #         atoms.append(atom.getAtomType())
       
    # print(atoms)
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
    
    # if verbose:
    #     print("---------------------")
    #     print("maxValues : ", maxValues)
    #     print("minValues : ", minValues)

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

    # print(overlap, " : is the overlap counted")
    # return volumeCompact, cubeSize/2, cubeSize/2, cubeSize/2
    return volumeCompact

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
        volm = atomList2emCube(atomList, pixelSize, densityNegative)
    else:
        volm = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if (not toCompact) and recenter:
        volm = recenterVolume(volm, densityNegative)
    
    if fname:
        volm.write(fname)
        print(f"MRC file is written in {fname}")
    
    return volm

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

# TODO
def getResolutionsByID(pdbIDList):
    pass

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

def wgetPDB2Volume(pdbID, pdbDir, volumeDir, toSave=True, cubeSize=0.0, toCompact=False, overwrite=False, verbose=False):
    """
    wgetPDB2Volume : Creates an PDB(CIF) file, EM file, MRC file from a PDB ID.
    @param overwrite : is for overwrite mrcfile(Volume2MRC).
    """
    # pdbDir should not include dangling /
    # densityNegative for default
    volumePath = f"{volumeDir}/{pdbID}.em"
    
    if verbose:
        print(f"wgetPDB2Volume is working with PDBID : {pdbID}")
    
    Path = wgetByPDBID(pdbID, pdbDir)
    resolution = getResolution(Path)
    _vol = cifpdb2em(Path, pixelSize=resolution, cubeSize=cubeSize, toCompact=toCompact, chain=None, fname=None, densityNegative=False, recenter=True)
    
    return _vol, resolution

def prepareCubeVolumes(pdbIDList, pdbDir, volumeDir, cubeSize=0.0, toCompact=False, overwrite=False, verbose=False):
    createdVolumes = []
    resolutionList = []
    for pdbID in pdbIDList:
        vol, resolution = wgetPDB2Volume(pdbID, pdbDir, volumeDir, toSave=False, cubeSize=cubeSize, toCompact=toCompact, overwrite=overwrite, verbose=verbose)
        createdVolumes.append(vol)
        resolutionList.append(resolution)
    return createdVolumes, resolutionList

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
def getMedadataJsonTemplate():
    ''' JSON format
    "v5": {
        "header" : "STRING",
        "pdbIDs" : "STRING LIST",
        "resolutions" : "FLOAT LIST",
        "particles" : [ 
            {
                "classNum" : "INT",
                "min"      : [ "x", "y", "z" ],
                "max"      : [ "x", "y", "z" ],
                "coord"    : [ [ "x", "y", "z" ] ]
            }
        ]
    }
    '''
    return {
        "header" : None,
        "pdbIDs": None,
        "resolutions": None,
        "particles": None,
    }
def minmaxUpdate(minList, maxList, coordList):
    for i in range(3):
        if minList[i] > coordList[i]:
            minList[i] = coordList[i]
        if maxList[i] < coordList[i]:
            maxList[i] = coordList[i]

def checkOverlap(minList, maxList, curList):
    if  minList[0] <= curList[0] and curList[0] <= maxList[0] and \
        minList[1] <= curList[1] and curList[1] <= maxList[1] and \
        minList[2] <= curList[2] and curList[2] <= maxList[2] :
        return True
    
    return False

def findMinMaxByList(list):
    min = [ list[0][0], list[0][1], list[0][2] ]
    max = [ list[0][0], list[0][1], list[0][2] ]
    for i in list:
        for j in range(3):
            if min[j] > i[j]:
                min[j] = i[j]
            if max[j] < i[j]:
                max[j] = i[j]
    return min, max

def makeScenarioByPDBIDs(pdbIDList, volumeDir, scenarioDir, scenarioIdentifier="noname", toSave=True, withClassMask=False, tomoSize=128, pfailedAttempts=9000, pparticleNum=1600, rotationStep=0, JSONCOMPACT=True, verbose=False):
    # cuboidalOccupancyList = [['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]]]
    startTime = time.time()
    scenarioMetaDataFile = f"{scenarioDir}/{scenarioIdentifier}.txt"
    scenarioJsonFile = f"{scenarioDir}/{scenarioIdentifier}.json"
    scenarioVolumeFile = f"{scenarioDir}/{scenarioIdentifier}.em"
    classmaskFile = f"{scenarioDir}/{scenarioIdentifier}_class_mask.em"

    jsonMetadataObject = getMedadataJsonTemplate()
    jsonMetadataObject["header"] = f"resolution : {10.0}, without Noise"
    jsonMetadataObject["pdbIDs"] = pdbIDList
    jsonMetadataObject["resolutions"] = [3.0]
    jsonMetadataObject["particles"] = []

    f = open(scenarioMetaDataFile, 'w')
    f.write("PDBID,centerX,centerY,centerZ,phi,theta,psi\n")
    fulltomX, fulltomY, fulltomZ = 2*tomoSize, 2*tomoSize, 2*tomoSize
    volume = vol(fulltomX, fulltomY, fulltomZ)
    volume.setAll(0.0)

    volumeTemplateList = []
    for pdbid in pdbIDList:
        currentTemplate = f"{volumeDir}/{pdbid}.em"
        volumeTemplateList.append(read(currentTemplate))

    if withClassMask:
        class_mask = vol(fulltomX, fulltomY, fulltomZ)
        class_mask.setAll(0)
    else:
        class_mask = None
    
    failedAttempts = 0
    particleNum = 0
    scenario = []

    classNum = np.random.randint(low=0, high=len(pdbIDList))
    currentVol = volumeTemplateList[classNum]

    if rotationStep != 0:
        rotatedVol, phi, theta, psi = compactRandomRotation(currentVol, rotationStep=rotationStep)
    else:
        rotatedVol, phi, theta, psi = currentVol, 0, 0, 0

    sizeX, sizeY, sizeZ = rotatedVol.sizeX(), rotatedVol.sizeY(), rotatedVol.sizeZ()
    x, y, z = np.random.randint(low=0, high=fulltomX-sizeX), np.random.randint(low=0, high=fulltomY-sizeY), np.random.randint(low=0, high=fulltomZ-sizeZ)
    centerX, centerY, centerZ = x+ sizeX/2, y+ sizeY/2, z+ sizeZ/2
    
    scenario.append([[x,y,z], [x+sizeX-1, y+sizeY-1, z+sizeZ-1]])

    occupyVoxels = []
    minCoord = [ 4*tomoSize, 4*tomoSize, 4*tomoSize ]
    maxCoord = [ -1, -1, -1]
    for i in range(sizeX):
        for j in range(sizeY):
            for k in range(sizeZ):
                curVal = rotatedVol.getV(i,j,k)
                if curVal != 0.0:
                    volume.setV( curVal , x+i, y+j, z+k)
                    minmaxUpdate(minCoord, maxCoord, [ x+i, y+j, z+k ])

                    if withClassMask:
                        class_mask.setV( classNum, x+i, y+j, z+k)
                    
                    occupyVoxels.append( [x+i, y+j, z+k] )
    
    if JSONCOMPACT:
        jsonMetadataObject["particles"].append( { "classNum" : classNum, "min" : minCoord, "max" : maxCoord } )
    else:
        jsonMetadataObject["particles"].append( { "classNum" : classNum, "min" : minCoord, "max" : maxCoord, "coord" : occupyVoxels } )
    # INCORPORTATE
    f.write(f"{pdbIDList[classNum]},{centerX},{centerY},{centerZ},{phi},{theta},{psi}\n")
    particleNum+=1
    rotatedVol = None
    rotFailNum = 0
    while failedAttempts < pfailedAttempts and particleNum < pparticleNum:
        occupyVoxels = []
        minCoord = [ 4*tomoSize, 4*tomoSize, 4*tomoSize ]
        maxCoord = [ -1, -1, -1]

        if verbose and failedAttempts%500 == 0 and failedAttempts != 0:
            print(f"... Now failed Attemps are {failedAttempts}")   
        if rotationStep != 0:
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
            classNum = np.random.randint(low=0, high=len(pdbIDList))
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
                        if curVal != 0: # TODO
                            volume.setV( curVal, x+i, y+j, z+k)
                            minmaxUpdate(minCoord, maxCoord, [ x+i, y+j, z+k ])
                            if withClassMask:
                                class_mask.setV( classNum, x+i, y+j, z+k)

                            occupyVoxels.append( [x+i, y+j, z+k] )
    
            if JSONCOMPACT:
                jsonMetadataObject["particles"].append( { "classNum" : classNum, "min" : minCoord, "max" : maxCoord } )
            else:
                jsonMetadataObject["particles"].append( { "classNum" : classNum, "min" : minCoord, "max" : maxCoord, "coord" : occupyVoxels } )
            particleNum+=1
            rotatedVol = None
            f.write(f"{pdbIDList[classNum]},{centerX},{centerY},{centerZ},{phi},{theta},{psi}\n")
            if verbose and particleNum % 50 == 0:
                print(f"... Particle Num : {particleNum}")  
    f.close()

    with open(scenarioJsonFile, "w") as json_file:
        json.dump(jsonMetadataObject, json_file)
    
    # --- META DATA ---
    appendMetaDataln(f"makeScenarioByPDBIDs {scenarioIdentifier} done - time elapsed : {time.time() - startTime}s")
    appendMetaDataln(f"-output file : {scenarioVolumeFile}, cubeSize : {fulltomX}x{fulltomY}x{fulltomZ}")
    appendMetaDataln(f"-with pdbIDList : {pdbIDList}")
    appendMetaDataln(f"-failedAttempts : {pfailedAttempts}, resultParticleNum : {particleNum}")

    print(f"----------- Scenario generation is done... with Particle Number {particleNum}----------")
    # -----------------
    if toSave:
        volume.write(scenarioVolumeFile)
        if withClassMask:
            class_mask.write(classmaskFile)
    else:
        return volume, class_mask

def makeGrandModelByPDBIDs(pdbIDList, pdbDir, volumeDir, scenarioDir, scenarioIdentifier="noname", withClassMask=True, toResolution=10.0, tomoSize=128, pfailedAttempts=8000, pparticleNum=1500, rotationStep=0, JSONCOMPACT=True, verbose=False):
    targetPath = f"{scenarioDir}/{scenarioIdentifier}.em"
    targetVoxelOccupyPath = f"{scenarioDir}/{scenarioIdentifier}_voxelOccupy.txt"
    maskPath = f"{scenarioDir}/{scenarioIdentifier}_class_mask.em"
    # First, PDB IDs -> PDB files -> Volume(.em) List
    print("makeGrandModelByPDBIDs : 1. prepare volume object from the Internet. -----------")
    volumes, resolutions = prepareCubeVolumes(pdbIDList, pdbDir=pdbDir, volumeDir=volumeDir, toCompact=True, overwrite=True, verbose=True)
    print("-- resolution list : ", resolutions)

    # Second, Resolution adjustment.
    print("makeGrandModelByPDBIDs : 2. resize volume object with respect to resolution. ---")
    f = open(targetVoxelOccupyPath, 'w')
    f.write("pdbID,occupyVoxelNum,resolution,toResolution\n")
    for pdbID, volume, resolution in zip(pdbIDList, volumes, resolutions):
        occupyVoxelNum = resolutionResizeUnity(volume, pdbID, 5.0, volumeDir, toResolution, verbose=verbose)
        f.write(f"{pdbID},{occupyVoxelNum},{resolution},{toResolution}\n")
    f.close()
    # Now, volume file is ready.
    print("makeGrandModelByPDBIDs : 3. make grandmodel. -----------------------------------")
    volume, class_mask = makeScenarioByPDBIDs(pdbIDList, volumeDir, toSave=False, withClassMask=withClassMask, scenarioDir=scenarioDir, scenarioIdentifier=scenarioIdentifier, tomoSize=tomoSize, pfailedAttempts=pfailedAttempts, pparticleNum=pparticleNum, rotationStep=rotationStep, JSONCOMPACT=JSONCOMPACT, verbose=verbose)

    volume.write(targetPath)
    class_mask.write(maskPath)
##########################################################################################################################################################################
#  Section for Simulation.
def compactRandomRotation(inputVolume, rotationStep = 1, toSave = False):
    #phi, theta, psi = np.random.randint(low=0, high=360, size=(3,)) # High exclusive
    phi, theta, psi = random.randrange(0, 359, rotationStep),random.randrange(0, 360, rotationStep),random.randrange(0, 360, rotationStep) # High inclusive. For step.
    from pytom_volume import rotate
    rotatedVolume = vol(inputVolume.sizeX(),inputVolume.sizeY(),inputVolume.sizeZ())
    rotatedVolume.setAll(0.0)
    rotate(inputVolume, rotatedVolume, int(phi), int(theta), int(psi))
    comp = makeCompact(rotatedVolume)
    return comp, phi, theta, psi

# FROM PDB ID -> Resolution corrected Compact Cuboid!
def resolutionResizeUnity(volume, identifier, resolution, outputDir, toResolution, verbose=False):
    if resolution > toResolution:
        raise RuntimeError(f"Target resolution {toResolution} is smaller than Resolution {resolution}. PDBID is {identifier}")
    
    outputVolumePath = f"{outputDir}/{identifier}.em"
    occupyVoxelNum = 0
    if resolution == toResolution:
        volume.write(outputVolumePath)
        return
    
    X, Y, Z = volume.sizeX(), volume.sizeY(), volume.sizeZ()
    resizedSizeX = math.ceil( (volume.sizeX() - 1) * resolution / toResolution ) 
    resizedSizeY = math.ceil( (volume.sizeY() - 1) * resolution / toResolution ) 
    resizedSizeZ = math.ceil( (volume.sizeZ() - 1) * resolution / toResolution ) 
    outputVolume = vol(resizedSizeX, resizedSizeY, resizedSizeZ)
    outputVolume.setAll(0.0)

    if verbose:
        print(f" BEFORE SIZE : {volume.sizeX()} x {volume.sizeY()} x {volume.sizeZ()}")
        print(f" AFTER SIZE : {outputVolume.sizeX()} x {outputVolume.sizeY()} x {outputVolume.sizeZ()}")
    
    for i in range(resizedSizeX):
        for j in range(resizedSizeY):
            for k in range(resizedSizeZ):
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
                                val = volume.getV(idxX, idxY, idxZ) * modfactor
                                curVal += val
                            except:
                                pass
                            idxZ += 1
                        idxY += 1
                    idxX += 1
                
                modfactor = resolution / toResolution
                voxelVal = curVal * modfactor
                if voxelVal != 0.0:
                    occupyVoxelNum += 1
                    outputVolume.setV(voxelVal , i, j, k)
    outputVolume.write(outputVolumePath)
    return occupyVoxelNum

def subtomoSampleSaver(tomoIdentifier, scenarioDir, subtomoIdentifier, subtomoDir, crowdLevel, SNR=1.0, generateNum=1, subtomoSizeX=50, subtomoSizeY=0, subtomoSizeZ=0):
    # Missing size is filled with X axis.
    if subtomoSizeY == 0:
        subtomoSizeY = subtomoSizeX
    if subtomoSizeZ == 0:
        subtomoSizeZ = subtomoSizeX
    
    scenarioVolume = read(f'{scenarioDir}/{tomoIdentifier}.em')
    
    particleCenters = []
    # Load txt file
    with open(f'{scenarioDir}/{tomoIdentifier}.txt') as scenarioParticleTxt:
        txt_contents = scenarioParticleTxt.readlines()
        particleTxtPattern = re.compile("(.*),(.*),(.*),(.*),(.*),(.*),(.*)")
        particleID = 0
        for line in txt_contents:
            parsedInfo = re.findall(particleTxtPattern, line)[0]
            particleCenter = [ parsedInfo[0], parsedInfo[1], parsedInfo[2], parsedInfo[3], particleID ]
            particleID += 1
            particleCenters.append(particleCenter)
    
    particleJsonList = []
    # Load json file
    with open(f'{scenarioDir}/{tomoIdentifier}.json') as scenarioJsonFile:
        json_object = json.load(scenarioJsonFile)
        particleJsonList = json_object['particles']
        pdbIDList = json_object['pdbIDs']

    for i in range(generateNum):
        particleList = []
        metadataParticleListTXT = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}_particles.txt"
        subtomoMRCFile = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}.mrc"
        subtomoJSONFile = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}.json"

        jsonMetadataObject = getMedadataJsonTemplate()
        jsonMetadataObject["header"] = f"{10.0}A/vx with white noise"
        jsonMetadataObject["pdbIDs"] = pdbIDList
        jsonMetadataObject["resolutions"] = [5.0]
        jsonMetadataObject["particles"] = []

        subtomo = vol(subtomoSizeX, subtomoSizeY, subtomoSizeZ)
        subtomo.setAll(0.0)

        if True: # type1 : just pick the center of the particle.
            # np.random.randint :: low exclusive high inclusive
            particleInfo = particleCenters[ np.random.randint(low = 0, high = len(particleCenters)) ]
            cX, cY, cZ = int(float(particleInfo[1])), int(float(particleInfo[2])), int(float(particleInfo[3]))
        
        if cX <= subtomoSizeX/2:
            lrX = 0
        elif cX >= scenarioVolume.sizeX() - subtomoSizeX/2:
            lrX = scenarioVolume.sizeX() - subtomoSizeX
        else:
            lrX = cX - subtomoSizeX//2
        if cY <= subtomoSizeY/2:
            lrY = 0
        elif cY >= scenarioVolume.sizeY() - subtomoSizeY/2:
            lrY = scenarioVolume.sizeY() - subtomoSizeY
        else:
            lrY = cY - subtomoSizeY//2
        if cZ <= subtomoSizeZ/2:
            lrZ = 0
        elif cZ >= scenarioVolume.sizeZ() - subtomoSizeZ/2:
            lrZ = scenarioVolume.sizeZ() - subtomoSizeZ
        else:
            lrZ = cZ - subtomoSizeZ//2
        
        for x in range(subtomoSizeX):
            for y in range(subtomoSizeY):
                for z in range(subtomoSizeZ):
                    getValue = scenarioVolume.getV( lrX + x , lrY + y, lrZ + z )
                    
                    if getValue != 0.0 :
                        subtomo.setV( getValue , x, y, z)
                        debug1 = []
                        for particle in particleJsonList:
                            
                            minList = particle['min']
                            maxList = particle['max']
                            curCord = [ lrX + x, lrY + y, lrZ + z]
                
                            debug1.append( [ minList, maxList, curCord ])
                            if checkOverlap(minList, maxList, curCord):
                                particleList.append( [  particleJsonList.index(particle), [x,y,z] ] )
                                break
                        else:
                            raise RuntimeError("SHOULD not happen.")
        
        subtomo = noiseApplier(subtomo, SNR=SNR)
        
        particleKeys = []
        particleDicts = []
        particleIDList = []

        for p in particleList:
            if p[0] in particleKeys:
                for pp in particleDicts:
                    if pp['particleID'] == p[0]:
                        pp['coord'].append(p[1])
                        break
            else:
                particleKeys.append( p[0] )
                particleDicts.append( { "particleID" : p[0], "min" : [] , "max" : [] ,  "coord": [ p[1] ] } )
        
        index = 1
        mrcSubtomo = volObj2Numpy(subtomo, floatMRC = True)
        f = open(metadataParticleListTXT, 'w')
        f.write("index,classNum,particleID\n")
        for p in particleDicts:
            p['min'], p['max'] = findMinMaxByList(p['coord'])
            p['classNum'] = particleJsonList[  p['particleID'] ]['classNum']
            f.write( f"{index},{p['classNum']},{p['particleID']}\n")
            # del p['particleID']
            mrcFileName = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}_subtomoParticle_{index}.mrc"
            subtomoP = newNumpyByXYZ(subtomoSizeX, subtomoSizeY, subtomoSizeZ, floatMRC = False) # Minimizing space -> integer.
            for particleCoord in p['coord']:
                subtomoP[particleCoord[0], particleCoord[1], particleCoord[2]] = 10 
            # Write each particle MRC metadata.
            with mrcfile.new(mrcFileName, overwrite=True) as mrc:
                mrc.set_data(subtomoP)
            index+=1
        f.close()

        jsonMetadataObject["particles"] = particleDicts

        print(f"--------- subtomogram generated : {i + 1}")
        
        # Write Subtomo MRC.
        with mrcfile.new(subtomoMRCFile, overwrite=True) as mrc:
            mrc.set_data(mrcSubtomo)
        # Write Subtomo Json metadata.
        with open(subtomoJSONFile, "w") as json_file:
            json.dump(jsonMetadataObject, json_file)

def subtomoDirectSaver(pdbIDList, volumesDir, subtomoIdentifier, subtomoDir, crowdLevel, SNR=1.0, generateNum=1, subtomoSizeX=50, subtomoSizeY=0, subtomoSizeZ=0):
    # Missing size is filled with X axis.
    if subtomoSizeY == 0:
        subtomoSizeY = subtomoSizeX
    if subtomoSizeZ == 0:
        subtomoSizeZ = subtomoSizeX
    
    volumeList = []
    for pdbID in pdbIDList:
        # TODO : also update the make grand!
    
    scenarioVolume = read(f'{scenarioDir}/{identifier}.em')
    
    particleCenters = []
    # Load txt file
    with open(f'{scenarioDir}/{identifier}.txt') as scenarioParticleTxt:
        txt_contents = scenarioParticleTxt.readlines()
        particleTxtPattern = re.compile("(.*),(.*),(.*),(.*),(.*),(.*),(.*)")
        particleID = 0
        for line in txt_contents:
            parsedInfo = re.findall(particleTxtPattern, line)[0]
            particleCenter = [ parsedInfo[0], parsedInfo[1], parsedInfo[2], parsedInfo[3], particleID ]
            particleID += 1
            particleCenters.append(particleCenter)
    
    particleJsonList = []
    # Load json file
    with open(f'{scenarioDir}/{identifier}.json') as scenarioJsonFile:
        json_object = json.load(scenarioJsonFile)
        particleJsonList = json_object['particles']
        pdbIDList = json_object['pdbIDs']

    for i in range(generateNum):
        particleList = []
        metadataParticleListTXT = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}_particles.txt"
        subtomoMRCFile = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}.mrc"
        subtomoJSONFile = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}.json"

        jsonMetadataObject = getMedadataJsonTemplate()
        jsonMetadataObject["header"] = f"{10.0}A/vx with white noise"
        jsonMetadataObject["pdbIDs"] = pdbIDList
        jsonMetadataObject["resolutions"] = [5.0]
        jsonMetadataObject["particles"] = []

        subtomo = vol(subtomoSizeX, subtomoSizeY, subtomoSizeZ)
        subtomo.setAll(0.0)

        if True: # type1 : just pick the center of the particle.
            # np.random.randint :: low exclusive high inclusive
            particleInfo = particleCenters[ np.random.randint(low = 0, high = len(particleCenters)) ]
            cX, cY, cZ = int(float(particleInfo[1])), int(float(particleInfo[2])), int(float(particleInfo[3]))
        
        if cX <= subtomoSizeX/2:
            lrX = 0
        elif cX >= scenarioVolume.sizeX() - subtomoSizeX/2:
            lrX = scenarioVolume.sizeX() - subtomoSizeX
        else:
            lrX = cX - subtomoSizeX//2
        if cY <= subtomoSizeY/2:
            lrY = 0
        elif cY >= scenarioVolume.sizeY() - subtomoSizeY/2:
            lrY = scenarioVolume.sizeY() - subtomoSizeY
        else:
            lrY = cY - subtomoSizeY//2
        if cZ <= subtomoSizeZ/2:
            lrZ = 0
        elif cZ >= scenarioVolume.sizeZ() - subtomoSizeZ/2:
            lrZ = scenarioVolume.sizeZ() - subtomoSizeZ
        else:
            lrZ = cZ - subtomoSizeZ//2
        
        for x in range(subtomoSizeX):
            for y in range(subtomoSizeY):
                for z in range(subtomoSizeZ):
                    getValue = scenarioVolume.getV( lrX + x , lrY + y, lrZ + z )
                    
                    if getValue != 0.0 :
                        subtomo.setV( getValue , x, y, z)
                        debug1 = []
                        for particle in particleJsonList:
                            
                            minList = particle['min']
                            maxList = particle['max']
                            curCord = [ lrX + x, lrY + y, lrZ + z]
                
                            debug1.append( [ minList, maxList, curCord ])
                            if checkOverlap(minList, maxList, curCord):
                                particleList.append( [  particleJsonList.index(particle), [x,y,z] ] )
                                break
                        else:
                            raise RuntimeError("SHOULD not happen.")
        
        subtomo = noiseApplier(subtomo, SNR=SNR)
        
        particleKeys = []
        particleDicts = []
        particleIDList = []

        for p in particleList:
            if p[0] in particleKeys:
                for pp in particleDicts:
                    if pp['particleID'] == p[0]:
                        pp['coord'].append(p[1])
                        break
            else:
                particleKeys.append( p[0] )
                particleDicts.append( { "particleID" : p[0], "min" : [] , "max" : [] ,  "coord": [ p[1] ] } )
        
        index = 1
        mrcSubtomo = volObj2Numpy(subtomo, floatMRC = True)
        f = open(metadataParticleListTXT, 'w')
        f.write("index,classNum,particleID\n")
        for p in particleDicts:
            p['min'], p['max'] = findMinMaxByList(p['coord'])
            p['classNum'] = particleJsonList[  p['particleID'] ]['classNum']
            f.write( f"{index},{p['classNum']},{p['particleID']}\n")
            # del p['particleID']
            mrcFileName = f"{subtomoDir}/{subtomoIdentifier}_{i + 1}_subtomoParticle_{index}.mrc"
            subtomoP = newNumpyByXYZ(subtomoSizeX, subtomoSizeY, subtomoSizeZ, floatMRC = False) # Minimizing space -> integer.
            for particleCoord in p['coord']:
                subtomoP[particleCoord[0], particleCoord[1], particleCoord[2]] = 10 
            # Write each particle MRC metadata.
            with mrcfile.new(mrcFileName, overwrite=True) as mrc:
                mrc.set_data(subtomoP)
            index+=1
        f.close()

        jsonMetadataObject["particles"] = particleDicts

        print(f"--------- subtomogram generated : {i + 1}")
        
        # Write Subtomo MRC.
        with mrcfile.new(subtomoMRCFile, overwrite=True) as mrc:
            mrc.set_data(mrcSubtomo)
        # Write Subtomo Json metadata.
        with open(subtomoJSONFile, "w") as json_file:
            json.dump(jsonMetadataObject, json_file)

##########################################################################################################################################################################
#  Main Code Workspace
##########################################################################################################################################################################
SHREC2021_FULL = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2", "5mrc" ]
SHREC2021_FULLex5mrc = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2" ]

SHREC2021_FULL_OCCLIST = [['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]
SHREC2021_FULLexc2L2S = [ "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr" ]
SHREC2021_1bxn = [ "1bxn" ]

if __name__ == "__main__":
    executionStart = time.time()
    #################### Workspace ##################    
    now = time.localtime()
    programTime = f"{now.tm_year}/{now.tm_mon}/{now.tm_mday} {now.tm_hour}:{now.tm_min}:{now.tm_sec}"
    appendMetaDataln(f"===> Scripts running : {programTime}")
    # Put some description.
    DESCRIPTION = "6__subtomo test data"
    appendMetaDataln(f"===> {DESCRIPTION}")

    #######################################################################################################################
    #ov = volumeOutliner("/cdata/scenario/0315_merge2.em", isFile=True, outlineValue = 500)
    #ov.write("/cdata/outlined/0315_merge2_utilstest.em")
    # for test dataset 2.
    # for - grandmodel.
    # makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0321_gmwoN_5.0compact", 10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=2200, rotationStep=2, JSONCOMPACT=True, verbose=True)
    # makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0321_gmwoN_5.0verbose", 10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=2200, rotationStep=2, JSONCOMPACT=False, verbose=True)
    # for - subtomogram.
    # previous approach.
    # sl, js = subtomoSampler("0321_gmwoN_5.0compact", "/cdata/scenario", 2, generateNum=3, subtomoSizeX=50)
    # volumeListWriter(sl, "/cdata/scenario", "0318_2_6gmwoN_5.0compact_subtomo", JSON=js)
    # unified approach.
    #subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0321_gmwoN_5.0compact_1", "/cdata/subtomo", 0, SNR=1.0, generateNum=10, subtomoSizeX=50)
    #subtomoSampleSaver("0321_gmwoN_5.0verbose", "/cdata/scenario", "0321_gmwoN_5.0verbose_1", "/cdata/subtomo", 0, SNR=1.0, generateNum=10, subtomoSizeX=50)
    
    # for test dataset 3. : mrc output / mrc metadata
    #subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0328_gmwoN_5.0wonoise_4", "/cdata/subtomo", 0, SNR=-1, generateNum=15, subtomoSizeX=50)
    #subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0328_gmwoN_5.0noise2.0_4", "/cdata/subtomo", 0, SNR=2.0, generateNum=15, subtomoSizeX=50)
    
    # for test dataset 4. : subtomo occupy percentage threshold.
    makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0328_5.0compact", 10.0, tomoSize=256, pfailedAttempts=11000, pparticleNum=2300, rotationStep=3, JSONCOMPACT=True, verbose=True)
    ################### Workspace Ended #############
    print(f" All of the Jobs completed with elapsed time : {time.time()-executionStart}")
    #################### Program Ended ##############
