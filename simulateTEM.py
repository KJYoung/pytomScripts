from pytom_volume import vol, initSphere
from pytom.basic.structures import WedgeInfo

# VKJY
from urllib.error import HTTPError
import numpy as np
import mrcfile
import os.path, os, wget, re, requests, math, random, time, json
import sys 
# My Other files
from type_convert import em2mrc, atomList2emCube, cifpdb2em, mrc2em, volume2MRC
from utils import volumeOutliner, volumeListWriter, volObj2Numpy, newNumpyByXYZ, volumeResizer, makeCompact
from pytomLib import recenterVolume, naivePDBParser, mmCIFParser, read, noiseApplier

PYTOM = 1 # PDB + CIF OK!
EMAN2 = 2 # ONLY PDB file!
SITUS = 3 # TOBEIMPLEMENTED(PDB only).

rootDIR = "/cdata"
pdbDataDIR = f"{rootDIR}/pdbData"
singleParticleEMCuboidDIR = f"{rootDIR}/singleParticleEM_cuboid"
singleParticleEMCubeDIR = f"{rootDIR}/singleParticleEM_cube"
singleParticleMRCCuboidDIR = f"{rootDIR}/singleParticleMRC_cuboid"
singleParticleMRCCubeDIR = f"{rootDIR}/singleParticleMRC_cube"
scenarioDIR = f"{rootDIR}/scenario"
metadataFILE = f"{rootDIR}/metadata.txt"

##########################################################################################################################################################################
#  Section for Utils.
##########################################################################################################################################################################
def appendMetaDataln(metadata):
    fmeta = open(metadataFILE, "a")
    fmeta.write(metadata + "\n")
    fmeta.close()
def getMedadataJsonTemplate():
    ''' JSON format
    "v6": {
        "header" : "STRING",
        "pdbIDs" : "STRING LIST",
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
        "particles": None,
    }
def minmaxUpdate(minList, maxList, coordList):
    for i in range(2):
        if minList[i] > coordList[i]:
            minList[i] = coordList[i]
        if maxList[i] < coordList[i]:
            maxList[i] = coordList[i]

##########################################################################################################################################################################
#  Section for wget.
##########################################################################################################################################################################
def getResolution(filePath):
    from pytom.tools.files import checkFileExists

    if not checkFileExists(filePath):
        raise RuntimeError('resolutionResize : input File not found! ', filePath)
    
    if filePath.endswith(".pdb"):
        resPatternPDB = re.compile("RESOLUTION\..*([0-9]+\.[0-9]+).*ANGSTROMS")
        f = open(filePath, 'r')
        pdbContent = f.read()
        f.close()
        regexList = re.findall(resPatternPDB, pdbContent)
        if len(regexList) == 0:
            raise RuntimeError('resolution is not applicable from the ', filePath)
        return float(regexList[0])
    elif filePath.endswith(".cif"):
        resPatternCIF = re.compile("_em_3d_reconstruction.resolution +([0-9]+\.[0-9]+)")
        f = open(filePath, 'r')
        cifContent = f.read()
        f.close()
        regexList = re.findall(resPatternCIF, cifContent)
        if len(regexList) == 0:
            raise RuntimeError('resolution is not applicable from the ', filePath)
        return float(regexList[0])
    else:
        raise RuntimeError('Unsupported file extenstion : ', filePath)

# Download pdbID.pdb or pdbID.cif(if .pdb not available) file into pdbDir
def wgetByPDBID(pdbID, pdbDir):
    pdbPath = f"{pdbDir}/{pdbID}.pdb"
    cifPath = f"{pdbDir}/{pdbID}.cif"
    pdbURL = f"https://files.rcsb.org/view/{pdbID}.pdb"
    cifURL = f"https://files.rcsb.org/view/{pdbID}.cif"

    if os.path.isfile(pdbPath):
        return pdbPath
    if os.path.isfile(cifPath):
        return cifPath

    if not requests.get(pdbURL).status_code == 200:
        if not requests.get(cifURL).status_code == 200:
            raise RuntimeError("Invalid pdb ID(Cannot download)! ", pdbID)
        isPDB = False
        URL = cifURL
        Path = cifPath
    else:
        isPDB = True
        URL = pdbURL
        Path = pdbPath
    wget.download(URL, out=Path)
    return Path

def em2mrc2norm2em(emPath, outPath, mrcPath, floatMRC=True, overwrite=True):
    inputVolume = read(emPath)
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    print(f"Volume dimension is initially... {x}x{y}x{z}")

    if floatMRC:
        volumeData = np.zeros([x, y, z], dtype = np.float32)
    else:
        volumeData = np.zeros([x, y, z], dtype = np.int8)
    
    for i in range(inputVolume.sizeX()):
        for j in range(inputVolume.sizeY()):
            for k in range(inputVolume.sizeZ()):
                try:
                    volumeData[i,j,k] = inputVolume.getV(i,j,k)
                except:
                    pass
    
    average = np.average(volumeData)
    std = np.std(volumeData)
    print(f" average : {average}, std : {std}")
    
    volumeData = (volumeData - average ) / 100 * std + 1
    average = np.average(volumeData)
    std = np.std(volumeData)
    print(f" average : {average}, std : {std}")

    # outputVolume = vol(x, y, z)
    # outputVolume.setAll(0.0)
    # for i in range(x):
    #     for j in range(y):
    #         for k in range(z):
    #             outputVolume.setV( volumeData[i,j,k] ,i,j,k)
    # outputVolume.write(outPath)
    with mrcfile.new(mrcPath, overwrite=overwrite) as mrc:
        mrc.set_data(volumeData)
        print(f"mrc data dimension is converted to... {mrc.data.shape}")
    mrc2em(mrcPath, outPath)
    return

"""
wgetPDB2Volume : 
    1. Creates an PDB(CIF) file by wgetByPDBID.
    2.  EM file, MRC file from a PDB ID.
@param pdbDir : should not include the dangling '/'
@param overwrite : is for overwrite mrcfile(Volume2MRC).
"""
def wgetPDB2Volume(pdbID, pdbDir, volumeDir, pixelSize=1, pdb2em=PYTOM, overwrite=False):
    # if type(pixelSize) != type(1):
    #     raise RuntimeError("wgetPDB2Volume : pixelSize should be Integer! ", pixelSize)
    print(f"wgetPDB2Volume is working with PDBID : {pdbID} -----------------------------------")

    volumePath = f"{volumeDir}/{pdbID}.em"
    mrcPath = f"{volumeDir}/{pdbID}.mrc"
    Path = wgetByPDBID(pdbID, pdbDir)
    resolution = getResolution(Path)

    if pdb2em == PYTOM:
        vol = cifpdb2em(Path, pixelSize=pixelSize, chain=None, fname=None)
    elif pdb2em == EMAN2:
        ## --omit OMIT : Randomly omit this percentage of atoms in the output map!
        if not Path.endswith(".pdb"):
            print(f"Warning : EMAN2 pdb2mrc only supports pdb type. But, {pdbID} doesn't have pdb format file. PYTOM mode conversion will be executed.")
            vol = cifpdb2em(Path, pixelSize=pixelSize, chain=None, fname=None)
        else:
            if resolution:
                eman2Command = f"e2pdb2mrc.py {Path} {mrcPath} --apix 1 --res {resolution} --center"
            else:
                eman2Command = f"e2pdb2mrc.py {Path} {mrcPath} --apix 1 --center"
            print(f"Executing ... {eman2Command}")
            os.system(eman2Command) # executing EMAN2
            mrc2em(mrcPath, volumePath) # mrc2em
            vol = read(volumePath) # em file to em object.
            vol = volumeResizer(vol, pixelSize) # average resizer.
            vol = makeCompact(vol) # compact!
    return vol, resolution

def prepareCubeVolumes(pdbIDList, pdbDir, volumeDir, pixelSize=1, pdb2em=PYTOM, overwrite=False):
    createdVolumes = []
    resolutionList = []
    for pdbID in pdbIDList:
        vol, resolution = wgetPDB2Volume(pdbID, pdbDir, volumeDir, pixelSize=pixelSize, pdb2em=pdb2em, overwrite=overwrite)
        createdVolumes.append(vol)
        resolutionList.append(resolution)
    return createdVolumes, resolutionList

##########################################################################################################################################################################
#  Section for Multi particle scenario.
##########################################################################################################################################################################
def makeImageByPDBIDs(pdbIDList, volumeDir, scenarioDir, scenarioIdentifier="noname", toSave=True, withClassMask=False, imageSize=128, pfailedAttempts=9000, pparticleNum=1600, rotationStep=0, verbose=False):
    # cuboidalOccupancyList = [['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]]]
    startTime = time.time()
    scenarioMetaDataFile = f"{scenarioDir}/{scenarioIdentifier}.txt"
    scenarioJsonFile = f"{scenarioDir}/{scenarioIdentifier}.json"
    scenarioVolumeFile = f"{scenarioDir}/{scenarioIdentifier}.em"
    classmaskFile = f"{scenarioDir}/{scenarioIdentifier}_class_mask.em"

    jsonMetadataObject = getMedadataJsonTemplate()
    jsonMetadataObject["header"] = f"{10.0} A/voxel, without Noise"
    jsonMetadataObject["pdbIDs"] = pdbIDList
    jsonMetadataObject["particles"] = []

    f = open(scenarioMetaDataFile, 'w')
    f.write("PDBID,centerX,centerY,phi,theta,psi\n")
    volume = vol(1, imageSize, imageSize)
    volume.setAll(0.0)

    volumeTemplateList = []
    for pdbid in pdbIDList:
        currentTemplate = f"{volumeDir}/{pdbid}.em"
        volumeTemplateList.append(read(currentTemplate))

    if withClassMask:
        class_mask = vol(1, imageSize, imageSize)
        class_mask.setAll(0)
    else:
        class_mask = None
    
    scenario = []

    classNum = np.random.randint(low=0, high=len(pdbIDList))
    currentVol = volumeTemplateList[classNum]

    # Rotate the particle.
    if rotationStep != 0:
        rotatedVol, phi, theta, psi = compactRandomRotation(currentVol, rotationStep=rotationStep)
    else:
        rotatedVol, phi, theta, psi = currentVol, 0, 0, 0
    # rotatedVol : rotated particle volume.

    sizeX, sizeY, sizeZ = rotatedVol.sizeX(), rotatedVol.sizeY(), rotatedVol.sizeZ()
    xProjectionHeight   = np.random.randint(low=0, high=sizeX) # [low, high)

    y, z = np.random.randint(low=0, high=imageSize-sizeY), np.random.randint(low=0, high=imageSize-sizeZ)
    centerY, centerZ = y+ sizeY/2, z+ sizeZ/2
    
    scenario.append([[y,z], [y+sizeY-1, z+sizeZ-1]])

    occupyVoxels = []
    minCoord = [ 4*imageSize, 4*imageSize ]
    maxCoord = [ -1, -1 ]
    for i in range(sizeY):
        for j in range(sizeZ):
                curVal = rotatedVol.getV(xProjectionHeight, i, j)
                if curVal != 0.0:
                    volume.setV( curVal , 0, y+i, z+j)
                    minmaxUpdate(minCoord, maxCoord, [ y+i, z+j ])

                    if withClassMask:
                        class_mask.setV( classNum + 1, 0, y+i, z+j)
                    occupyVoxels.append( [ y+i, z+j ] )
    
    jsonMetadataObject["particles"].append( { "classNum" : classNum, "min" : minCoord, "max" : maxCoord } )

    # INCORPORTATE
    f.write(f"{pdbIDList[classNum]},{centerY},{centerZ},{phi},{theta},{psi}\n")
    rotatedVol = None
    rotFailNum = 0

    failedAttempts = 0
    particleNum = 1
    while failedAttempts < pfailedAttempts and particleNum < pparticleNum:
        occupyVoxels = []
        minCoord = [ 4*imageSize, 4*imageSize ]
        maxCoord = [ -1, -1 ]

        if verbose and failedAttempts%2000 == 0 and failedAttempts != 0:
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
        xProjectionHeight   = np.random.randint(low=0, high=sizeX) # [low, high)
        y, z = np.random.randint(low=0, high=imageSize-sizeY), np.random.randint(low=0, high=imageSize-sizeZ)
        centerY, centerZ = y+ sizeY/2, z + sizeZ/2
        
        isOccupied = False
        for existingItem in scenario:
            if (y <= existingItem[1][0] and y+sizeY >= existingItem[0][0]) and (z <= existingItem[1][1] and z+sizeZ >= existingItem[0][1]):
                # print(f"Failed {failedAttempts} : {y}<={existingItem[1][0]} and {y}+{sizeY} >= {existingItem[0][0]}) and ({z} <= {existingItem[1][1]} and {z}+{sizeZ} >= {existingItem[0][1]})")
                failedAttempts+=1
                isOccupied = True
                break; 
        if isOccupied == False:
            scenario.append([[y, z], [y+sizeY-1, z+sizeZ-1]])
            for i in range(sizeY):
                for j in range(sizeZ):
                        curVal = rotatedVol.getV(xProjectionHeight, i, j)
                        if curVal != 0: # TODO
                            volume.setV( curVal, 0, y+i, z+j)
                            minmaxUpdate(minCoord, maxCoord, [ y+i, z+j ])
                            if withClassMask:
                                class_mask.setV( classNum + 1, 0, y+i, z+j)
                            occupyVoxels.append( [ y+i, z+j ] )
    
            jsonMetadataObject["particles"].append( { "classNum" : classNum, "min" : minCoord, "max" : maxCoord } )
            particleNum+=1
            rotatedVol = None
            f.write(f"{pdbIDList[classNum]},{centerY},{centerZ},{phi},{theta},{psi}\n")
            if verbose and particleNum % 500 == 0:
                print(f"... Particle Num : {particleNum}")  
    f.close()

    with open(scenarioJsonFile, "w") as json_file:
        json.dump(jsonMetadataObject, json_file)
    
    # --- META DATA ---
    appendMetaDataln(f"makeScenarioByPDBIDs {scenarioIdentifier} done - time elapsed : {time.time() - startTime}s")
    appendMetaDataln(f"-output file : {scenarioVolumeFile}, imageSize : {imageSize}x{imageSize}")
    appendMetaDataln(f"-with pdbIDList : {pdbIDList}")
    appendMetaDataln(f"-Parameter - pfailedAttempts : {pfailedAttempts}, pParticleNum : {pparticleNum}")
    appendMetaDataln(f"-Result - failedAttempts : {failedAttempts}, resultParticleNum : {particleNum}")

    print(f"----------- Scenario generation is done... with Particle Number {particleNum}----------")
    # -----------------
    if toSave:
        volume.write(scenarioVolumeFile)
        if withClassMask:
            class_mask.write(classmaskFile)
    return volume, class_mask, particleNum
    
def makeVolumeByPDBIDs(pdbIDList, pdbDir, volumeDir, pixelSize=10):
   # PDB IDs -> PDB files -> Volume(.em) List
    print("makeVolumeByPDBIDs : prepare volume object from the Internet. -----------")
    volumes, _resolutions = prepareCubeVolumes(pdbIDList, pdbDir=pdbDir, pixelSize=pixelSize, volumeDir=volumeDir, overwrite=True)
    for pdbID, volume in zip(pdbIDList, volumes):
        volume.write(f"{volumeDir}/{pdbID}.em")

def simulateTEM( pdbIDList, pdbDir, volumeDir, scenarioDir, scenarioIdentifier="noname", withClassMask=True, 
                            newVolume=True, pixelSize=10, imageSize=128, pfailedAttempts=8000, pparticleNum=1500, rotationStep=0, 
                            pdb2em=PYTOM, JSONCOMPACT=True, verbose=False):
    targetPath = f"{scenarioDir}/{scenarioIdentifier}.em"
    targetVoxelOccupyPath = f"{scenarioDir}/{scenarioIdentifier}_voxelOccupy.txt"
    maskPath = f"{scenarioDir}/{scenarioIdentifier}_class_mask.em"
    # PDB IDs -> PDB files -> Volume(.em) List
    print("makeGrandModelByPDBIDs : 1. prepare volume object from the Internet. -----------")
    if newVolume:
        volumes, _resolutions = prepareCubeVolumes(pdbIDList, pdbDir=pdbDir, pixelSize=pixelSize, volumeDir=volumeDir, pdb2em=pdb2em, overwrite=True)
        for pdbID, volume in zip(pdbIDList, volumes):
            volume.write(f"{volumeDir}/{pdbID}.em")

    # Now, volume file is ready.
    print("makeGrandModelByPDBIDs : 2. make grandmodel. -----------------------------------")
    volume, class_mask, particleNum = makeImageByPDBIDs(pdbIDList, volumeDir, toSave=False, withClassMask=withClassMask, scenarioDir=scenarioDir, scenarioIdentifier=scenarioIdentifier, imageSize=imageSize, pfailedAttempts=pfailedAttempts, pparticleNum=pparticleNum, rotationStep=rotationStep, verbose=verbose)

    volume.write(targetPath)
    class_mask.write(maskPath)
    return particleNum
##########################################################################################################################################################################
#  Section for Simulation.
def compactRandomRotation(inputVolume, rotationStep = 1, toSave = False):
    #phi, theta, psi = np.random.randint(low=0, high=360, size=(3,)) # High exclusive
    phi, theta, psi = random.randrange(0, 359, rotationStep),random.randrange(0, 359, rotationStep),random.randrange(0, 359, rotationStep) # High inclusive. For step.
    from pytom_volume import rotate
    rotatedVolume = vol(inputVolume.sizeX(),inputVolume.sizeY(),inputVolume.sizeZ())
    rotatedVolume.setAll(0.0)
    rotate(inputVolume, rotatedVolume, int(phi), int(theta), int(psi))
    comp = makeCompact(rotatedVolume)
    return comp, phi, theta, psi

##########################################################################################################################################################################
#  Main Code Workspace
##########################################################################################################################################################################
if __name__ == "__main__":
    executionStart = time.time()
    #################### Workspace ##################    
    now = time.localtime()
    programTime = f"{now.tm_year}/{now.tm_mon}/{now.tm_mday} {now.tm_hour}:{now.tm_min}:{now.tm_sec}"
    appendMetaDataln(f"===> Scripts running : {programTime}")
    # Put some description.
    DESCRIPTION = "8_InsilicoTEM replace?"
    appendMetaDataln(f"===> {DESCRIPTION}")

    # generate emByEMAN2/2wrj_norm2.em ,mrc.
    # em2mrc2norm2em("/cdata/tomsimDIR/emByEMAN2/2wrj.em", "/cdata/tomsimDIR/emByEMAN2/2wrj_norm2.em", "/cdata/tomsimDIR/emByEMAN2/2wrj_norm2.mrc", floatMRC=True, overwrite=True)
    
    # simulateTEM(["2wrj"], "/cdata/tomsimDIR/pdbData", "/cdata/tomsimDIR/emByEMAN2", "/cdata/tomsimDIR/scenario", "null", pixelSize=2, imageSize=4096, pfailedAttempts=0, pparticleNum=0, rotationStep=0, pdb2em=EMAN2, verbose=True)
    
    print("makeGrandModelByPDBIDs : 2. make grandmodel. -----------------------------------")
    for num in range(250):
        identifier = f"0714_EMAN2_2wrj_{num}"
        pNum = np.random.randint(low=60, high=120)
        makeImageByPDBIDs(["2wrj"], "/cdata/tomsimDIR/emByEMAN2", toSave=True, withClassMask=False, scenarioDir="/cdata/tomsimDIR/scenario", scenarioIdentifier=identifier, imageSize=4096, pfailedAttempts=50, pparticleNum=pNum, rotationStep=5, verbose=True)
        volume2MRC(f"/cdata/tomsimDIR/scenario/{identifier}.em", f"/cdata/tomsimDIR/scenario/{identifier}.mrc", floatMRC=True, overwrite=False, verbose=True)
        num2 = num+1
        print("{}/250".format(num2), file=sys.stderr, end='\r')
    # print(particleNum)
    # for test dataset 7. : EMAN2 + various crowd level.
    # test dataset 7 scenario generation. -------------
    # particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523(PDB)_crowd6(Max)_EMAN2", pixelSize=10, tomoSize=256, pfailedAttempts=2000000, pparticleNum=9999999, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
    # print(particleNum)
    # particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523(PDB)_crowd5_EMAN2", pixelSize=10, tomoSize=256, pfailedAttempts=300000, pparticleNum=9999999, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
    # particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd4", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=10000, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
    # particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd3", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=4000, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
    # particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd2", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=1500, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
    # particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd1", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=600, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
    # volume2MRC("/cdata/scenario/0523(PDB)_crowd6(Max)_EMAN2.em", "/cdata/scenario/0523(PDB)_crowd6(Max)_EMAN2.mrc", floatMRC=True, overwrite=False, verbose=True)
    # volume2MRC("/cdata/scenario/0523(PDB)_crowd5(Max)_EMAN2.em", "/cdata/scenario/0523(PDB)_crowd5(Max)_EMAN2.mrc", floatMRC=True, overwrite=False, verbose=True)
    # volume2MRC("/cdata/scenario/0523_crowd4.em", "/cdata/scenario/0523_crowd4.mrc", floatMRC=True, overwrite=False, verbose=True)
    # volume2MRC("/cdata/scenario/0523_crowd3.em", "/cdata/scenario/0523_crowd3.mrc", floatMRC=True, overwrite=False, verbose=True)
    # volume2MRC("/cdata/scenario/0523_crowd2.em", "/cdata/scenario/0523_crowd2.mrc", floatMRC=True, overwrite=False, verbose=True)
    # volume2MRC("/cdata/scenario/0523_crowd1.em", "/cdata/scenario/0523_crowd1.mrc", floatMRC=True, overwrite=False, verbose=True)
    # Generate subtomo from those template scenario.
    print("Generating subtomograms..")

    ################### Workspace Ended #############
    print(f" All of the Jobs completed with elapsed time : {time.time()-executionStart}")
    #################### Program Ended ##############
