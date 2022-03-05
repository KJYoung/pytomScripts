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
# ------------
# Custom volume generator from pdb, cif.
# Compact volume generaor is added.
#
#

rootDIR = "/cdata"
pdbDataDIR = f"{rootDIR}/pdbData"
singleParticleEMCuboidDIR = f"{rootDIR}/singleParticleEM_cuboid"
singleParticleMRCCuboidDIR = f"{rootDIR}/singleParticleMRC_cuboid"
##########################################################################################################################################################################
#  Section for PDB ID to volume.
##########################################################################################################################################################################
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

def pdb2em(pdbPath, pixelSize, cubeSize=0.0, toCompact=False, chain=None, densityNegative=False, fname='', recenter=True):
    """
    pdb2em: Creates an volume out of a PDB file
    @param pdbPath: Path to PDB file or PDB id for online download
    @param pixelSize: The pixel size to convert to 
    @param cubeSize: Resulting cube size
    @param toCompact: Option for compact cuboid
    @return: A volume
    @author: Thomas Hrabe & Luis Kuhn
    """
    from math import floor
    from pytom_volume import vol

    atomList = naivePDBParser(pdbPath, chain)

    vol = None
    compactX, compactY, compactZ = 0.0, 0.0, 0.0
    if toCompact:
        vol, compactX, compactY, compactZ = atomList2emCompact(atomList, pixelSize, densityNegative)
    else:
        vol = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if (not toCompact) and recenter:
        vol = recenterVolume(vol, densityNegative)

    if fname:
        vol.write(fname)
        print("MRC file is written in {}".format(fname))

    return vol, compactX, compactY, compactZ

def mmCIF2em(mmCIFPath, pixelSize, cubeSize=0.0, toCompact=False, chain=None, densityNegative=False, fname='', recenter=True):
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

    atomList = mmCIFParser(mmCIFPath, chain)

    vol = None
    compactX, compactY, compactZ = 0.0, 0.0, 0.0
    if toCompact:
        vol, compactX, compactY, compactZ = atomList2emCompact(atomList, pixelSize, densityNegative)
    else:
        vol = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if (not toCompact) and recenter:
        vol = recenterVolume(vol, densityNegative)
    
    if fname:
        vol.write(fname)
        print(f"MRC file is written in {fname}")
    
    return vol, compactX, compactY, compactZ

def volume2MRCconverter(volPath, mrcPath, overwrite=False, verbose=False):
    inputVolume = read(volPath)
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    if verbose:
        print(f"Volume dimension is initially... {x}x{y}x{z}")
    volumeData = np.empty([x,y,z], dtype=np.int8)
    for i in range(inputVolume.sizeX()):
        for j in range(inputVolume.sizeY()):
            for k in range(inputVolume.sizeZ()):
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

    response = requests.get(pdbURL)
    if not response.status_code == 200:
        response = requests.get(cifURL)
        if not response.status_code == 200:
            print("Invalid pdb ID maybe.")
            return
        if not templateExists:
            wget.download(cifURL, out=cifPath)
        resPatternCIF = re.compile("_em_3d_reconstruction.resolution +([0-9]+\.[0-9]+)")

        f = open(cifPath, 'r')
        cifContent = f.read()
        f.close()

        cifPixResolution = re.findall(resPatternCIF, cifContent)[0]
        print("cif resolution of {} is {}".format(pdbID, cifPixResolution))
        _vol, compactX, compactY, compactZ = mmCIF2em(cifPath, pixelSize=float(cifPixResolution), cubeSize=cubeSize, toCompact=toCompact, chain=None, fname=volumePath, densityNegative=False, recenter=True)
        if generateMRC:
            volume2MRCconverter(volumePath, mrcPath, overwrite=overwrite)
        return compactX, compactY, compactZ
    if not templateExists:
        wget.download(pdbURL, out=pdbPath)
    resPatternPDB = re.compile("RESOLUTION\..*([0-9]+\.[0-9]+).*ANGSTROMS")
    
    f = open(pdbPath, 'r')
    pdbContent = f.read()
    f.close()

    pdbPixResolution = re.findall(resPatternPDB, pdbContent)[0]
    print("pdb resolution of {} is {}".format(pdbID, pdbPixResolution))
    _vol, compactX, compactY, compactZ = pdb2em(pdbPath, pixelSize=float(pdbPixResolution), cubeSize=cubeSize, toCompact=toCompact, chain=None, fname=volumePath, densityNegative=False, recenter=True)
    if generateMRC:
        volume2MRCconverter(volumePath, mrcPath, overwrite=overwrite)
    return compactX, compactY, compactZ

##########################################################################################################################################################################
#  Section for Multi particle scenario.
##########################################################################################################################################################################
import random
def occupancyCalculator(pdbIDList, pdbDir, volumeDir, mrcDir, overwrite=False, verbose=False):
    cuboidalOccupancyList = []
    for pdbID in pdbIDList:
        x, y, z = wgetPDB2Volume(pdbID, pdbDir=pdbDir, volumeDir=volumeDir, mrcDir=mrcDir, toCompact=True, generateMRC=False, overwrite=overwrite, verbose=verbose)
        cuboidalOccupancyList.append([pdbID, [x, y, z]])
        if verbose:
            print(f"Volume building with {pdbID} is done...")
    return cuboidalOccupancyList
[['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]


def makeScenarioByPDBIDs(pdbIDList, scenarioPATH, occList=None, cubeSize=256, overwrite=False, verbose=False):
    # Prepare Volumes.
    # cuboidalOccupancyList = [['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]]]
    if occList:
        cuboidalOccupancyList = occList
    else:
        cuboidalOccupancyList = occupancyCalculator(pdbIDList, overwrite, verbose)
    if verbose:
        print("----------- occupancy List is done... -----------")
        #print(cuboidalOccupancyList)
    failedAttempts = 0
    particleNum = 0
    scenario = []
    while failedAttempts < 1600000 and particleNum < 1600:
        if verbose and failedAttempts%400000 == 0:
            print(f"... Now failed Attemps are {failedAttempts}")   
        classNum = random.randint(0, len(cuboidalOccupancyList)-1)
        classOccpy = cuboidalOccupancyList[classNum]
        x, y, z = random.randint(0, cubeSize-1-classOccpy[1][0]), random.randint(0, cubeSize-1-classOccpy[1][1]), random.randint(0, cubeSize-1-classOccpy[1][2])
        
        if len(scenario) == 0:
            scenario.append([classNum, [x,y,z], classOccpy[1]])
            particleNum+=1
        else:
            for existingItem in scenario:
                if (x <= existingItem[1][0]+existingItem[2][0] and x+classOccpy[1][0] >= existingItem[1][0]) and \
                    (y <= existingItem[1][1]+existingItem[2][1] and y+classOccpy[1][1] >= existingItem[1][1]) and \
                    (z <= existingItem[1][2]+existingItem[2][2] and z+classOccpy[1][2] >= existingItem[1][2]):
                    failedAttempts+=1
                    break; 
                if existingItem == scenario[-1]:
                    scenario.append([classNum, [x, y, z], classOccpy[1]])
                    particleNum+=1
    
    print(f"----------- Scenario generation is done... with Particle Number {particleNum}----------")
    # print(scenario)
    # [[0, [305, 116, 185], [46, 32, 38]], [3, [404, 398, 234], [34, 69, 67]], [3, [158, 219, 313], [34, 69, 67]], [4, [31, 471, 392], [31, 36, 44]], [7, [215, 168, 1], [45, 41, 54]], [1, [445, 347, 63], [39, 32, 37]], [2, [75, 46, 145], [41, 34, 27]], [2, [361, 299, 102], [41, 34, 27]], [5, [267, 6, 446], [25, 36, 21]], [1, [118, 81, 469], [39, 32, 37]]]

    volume = vol(cubeSize, cubeSize, cubeSize)
    volume.setAll(0.0)

    for items in scenario:
        particleVol = read(f"/cdata/volumesCompact/{cuboidalOccupancyList[items[0]][0]}.em")
        x, y, z = particleVol.sizeX(), particleVol.sizeY(), particleVol.sizeZ()
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    if volume.getV(items[1][0]+i, items[1][1]+j, items[1][2]+k) != 0:
                        print("ERROR : PARTICLE OVERAPPED!")
                        print(f"debug information : {items}")
                    volume.setV( particleVol.getV(i,j,k), items[1][0]+i, items[1][1]+j, items[1][2]+k)
        
        # with mrcfile.new(mrcPath, overwrite=overwrite) as mrc:
        #     mrc.set_data(volumeData)
        #     print(f"mrc data dimension is converted to... {mrc.data.shape}")
        # return
        # COPY
    
    volume.write(scenarioPATH)
    #END

##########################################################################################################################################################################
#  Section for Simulation.

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
    return s

##########################################################################################################################################################################
#  Section for Toy Test Generation.


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


def syntheticGen1(outputPath):
    volume = vol(64, 64, 64)
    volume.setAll(0.0)
    for i in range(64):
        volume.setV(20, 32, 32, i)

        volume.setV(20, i, 32, 32)
        volume.setV(20, i, 31, 31)
        volume.setV(20, i, 33, 33)

        volume.setV(20, 32, i, 32)
        volume.setV(20, 29, i, 29)
        volume.setV(20, 35, i, 35)
    
    mat = np.zeros([64, 64])
    for i in range(64):
        mat[i][i] = 20
    matrixToVolLayerX(volume, mat, 40)
    volume.write(outputPath)

def getWeightedSum(inputVolume):
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    sumOfPoints = 0
    for i in range(x):
        for j in range(y):
            for k in range(z):
                sumOfPoints += inputVolume.getV(i,j,k)    
    return sumOfPoints

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

def rotationIntegrityTest2(volumePath):
    em = read(volumePath)
    angles = range(0, 360, 20)
    f = open("/cdata/0305/rotateIntegrity_cube.txt", 'w')
    for angleX in angles:
        print(f"angleX : {angleX} is started...")
        for angleY in angles:
            print(f"--angleY : {angleY} is started...")
            for angleZ in angles:
                from pytom_volume import rotate
                rotatedCopy = vol(em.sizeX(),em.sizeY(),em.sizeZ())
                rotatedCopy.setAll(0.0)
                rotate(em, rotatedCopy, angleX, angleY, angleZ)
                compactRotated = makeCompact(rotatedCopy)
                result = getWeightedSum(compactRotated)
                f.write(f"{angleX},{angleY},{angleZ},{compactRotated.sizeX()},{compactRotated.sizeY()},{compactRotated.sizeZ()},{result}\n")
    f.close()

def rotationIntegrityTest1(volumePath): #1bxn 45 45 37
    em = read(volumePath)
    #print(em.sizeX(),em.sizeY(),em.sizeZ())
    angles = range(0, 360, 20)
    f = open("/cdata/0305/rotateIntegrity_cuboid.txt", 'w')
    for angleX in angles:
        print(f"angleX : {angleX} is started...")
        for angleY in angles:
            print(f"--angleY : {angleY} is started...")
            for angleZ in angles:
                from pytom_volume import rotate
                rotatedCopy = vol(em.sizeX(),em.sizeY(),em.sizeZ())
                rotatedCopy.setAll(0.0)
                rotate(em, rotatedCopy, angleX, angleY, angleZ)
                result = getWeightedSum(rotatedCopy)
                f.write(f"{angleX},{angleY},{angleZ},{result}\n")
    f.close()
##########################################################################################################################################################################
SHREC2021_FULL = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2", "5mrc" ]
SHREC2021_FULL_OCCLIST = [['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]
SHREC2021_FULLexc2L2S = [ "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr" ]
SHREC2021_1bxn = [ "1bxn" ]

#for pdbID in SHREC2021_1bxn:
#    wgetPDB2Volume(pdbID)
if __name__ == "__main__":
    import time
    #################### Workspace ##################
    executionStart = time.time()

    #volume2MRCconverter('/cdata/volumesV/1bxn_200.mrc', '/cdata/mrcs/1bxn_200.mrc', verbose=True)
    #wgetPDB2Volume('1bxn', toCompact=True, generateMRC=True, verbose=True)
    dateDir = "/cdata/0305"
    scenarioFile = f'{dateDir}/scenario256_1.em'
    simulatedFile = f'{dateDir}/simulated(256_1).em'

    # Rotation integrity Test
    #wgetPDB2Volume('1bxn', pdbDir=pdbDataDIR, volumeDir=singleParticleEMCuboidDIR, mrcDir=None, toCompact=True, overwrite=True, verbose=True)
    #compactCuboid2rotateCube("/cdata/singleParticleEM_cuboid/1bxn.em", "/cdata/singleParticleEM_cube/1bxn.em")
    wgetPDB2Volume('5mrc', pdbDir=pdbDataDIR, volumeDir=singleParticleEMCuboidDIR, mrcDir=None, toCompact=True, overwrite=True, verbose=True)
    #rotationIntegrityTest2("/cdata/singleParticleEM_cube/1bxn.em")

    #syntheticGen1(f'{dateDir}/rotest4X.em')
    #customSimulation(f'{dateDir}/rotest4X.em', f'{dateDir}/rotest4X_0_0_30.em', 100, rotation=[0,0,30])
    #customSimulation(f'{dateDir}/rotest4X.em', f'{dateDir}/rotest4X_0_30_0.em', 100, rotation=[0,30,0])
    #customSimulation(f'{dateDir}/rotest4X.em', f'{dateDir}/rotest4X_30_0_0.em', 100, rotation=[30,0,0])


    # makeScenarioByPDBIDs(SHREC2021_FULL, scenarioPATH=scenarioFile, occList=SHREC2021_FULL_OCCLIST, cubeSize=256,overwrite=True, verbose=True)
    # customSimulation(scenarioFile, simulatedFile, 0.1, wedgeAngle=70)
    print(f" All of the Jobs completed with elapsed time : {time.time()-executionStart}")
    #################### Program Ended ##############