from pytom.basic.files import read
from pytom_volume import vol
import json
import mrcfile
import numpy as np
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

# isFile = True -> filePath
# isFile = False -> list of volume
def volumeOutliner(inputVolumes, isFile=False, outlineValue=10):
    if isFile:
        inputVolumes = [ read(inputVolumes) ]
     
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
    
    if isFile:
        return inputVolumes[0]

def volumeListWriter(inputVolumes, outputDir, Description, JSON=None):
    index = 0
    if JSON:
        for inputVolume, jsonOb in zip(inputVolumes, JSON):
            inputVolume.write(f"{outputDir}/{Description}_{index}.em")
            with open(f"{outputDir}/{Description}_{index}.json", "w") as json_file:
                json.dump(jsonOb, json_file)

            index += 1
    else:
        for inputVolume in inputVolumes:
            inputVolume.write(f"{outputDir}/{Description}_{index}.em")
            
            index += 1

def volume2MRC(volPath, mrcPath, floatMRC=False, overwrite=False, verbose=False):
    inputVolume = read(volPath)
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    if verbose:
        print(f"Volume dimension is initially... {x}x{y}x{z}")
    
    if floatMRC:
        volumeData = np.zeros([x, y, z], dtype = np.float32)
    else:
        volumeData = np.zeros([x, y, z], dtype = np.int8)
    
    for i in range(inputVolume.sizeX()):
        for j in range(inputVolume.sizeY()):
            for k in range(inputVolume.sizeZ()):
                #print(type(inputVolume.getV(i,j,k)))
                volumeData[i,j,k] = inputVolume.getV(i,j,k)
    
    with mrcfile.new(mrcPath, overwrite=overwrite) as mrc:
        mrc.set_data(volumeData)
        print(f"mrc data dimension is converted to... {mrc.data.shape}")
    return

def volObj2Numpy(inputVolume, floatMRC=False):
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    
    if floatMRC:
        volumeData = np.zeros([x, y, z], dtype = np.float32)
    else:
        volumeData = np.zeros([x, y, z], dtype = np.int8)
    
    for i in range(inputVolume.sizeX()):
        for j in range(inputVolume.sizeY()):
            for k in range(inputVolume.sizeZ()):
                volumeData[i,j,k] = inputVolume.getV(i,j,k)
    
    return volumeData

def newNumpyByXYZ(x, y, z, floatMRC=False):
    if floatMRC:
        volumeData = np.zeros([x, y, z], dtype = np.float32)
    else:
        volumeData = np.zeros([x, y, z], dtype = np.int8)
    
    return volumeData

