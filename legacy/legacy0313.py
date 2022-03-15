# Rotate each axis : value explode..
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
