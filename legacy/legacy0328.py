# Dimension concept already exist in the PyTom..
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

    # for i in range(len(atomList)):
    #     centroidX += atomList[i].getX()
    #     centroidY += atomList[i].getY()
    #     centroidZ += atomList[i].getZ()

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
        #print("centroids : ", [centroidX, centroidY, centroidZ])
        #print("centers   : ", [centerX, centerY, centerZ])
        #print("shifts    : ", [shiftX, shiftY, shiftZ])

    #################### COMPACT VOLUME ####################
    compactX, compactY, compactZ = maxValues[0]-minValues[0], maxValues[1]-minValues[1], maxValues[2]-minValues[2]
    volumeCompact = vol(compactX+1, compactY+1, compactZ+1)
    # add 1 is crucial, basically
    volumeCompact.setAll(0.0)

    for i in range(len(atomList)):
        x = int(atomList[i].getX()) - minValues[0]
        y = int(atomList[i].getY()) - minValues[1]
        z = int(atomList[i].getZ()) - minValues[2]

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

    return volumeCompact, compactX, compactY, compactZ
