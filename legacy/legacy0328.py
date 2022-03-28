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