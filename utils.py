################################
# matrixToVolLayerX, Y, Z
# volumeOutliner
# volumeListWriter
################################
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