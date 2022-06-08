
# em2mrc
# atomList2emCube
# cifpdb2em
# mrc2em
# volume2mrc
from pytom_volume import vol
from pytomLib import recenterVolume, naivePDBParser, mmCIFParser, read
import numpy as np
import mrcfile

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
    #     print("maxValues : ", maxValues)
    #     print("minValues : ", minValues)

    #################### COMPACT CUBE VOLUME ####################
    compactX, compactY, compactZ = maxValues[0]-minValues[0], maxValues[1]-minValues[1], maxValues[2]-minValues[2]
    cubeSize = int(sqrt( compactX**2 + compactY**2 + compactZ**2 ))
    # if cubeSize % 2 == 0:
    #     cubeSize += 1 # Let cubeSize to be odd.
    # Let cubeSize to be even!
    if cubeSize % 2 != 0:
        cubeSize += 1 # Let cubeSize to be even.
    
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
            volumeCompact.setV(currentValue + mass / pixelSize / pixelSize / pixelSize, x, y, z)
            
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

def cifpdb2em(inputPath, pixelSize, chain=None, fname=''):
    """
    cifpdb2em: Creates an volume out of a mmCIF file or pdb file
    @param inputPath: Path to mmCIF file 
    @param pixelSize: The pixel size to convert to 
    @param cubeSize: Resulting cube size
    @param toCompact: Option for compact cuboid
    @param recenter: center the volume - default and only True supported.
    @param densityNegartive: minus all the value - default and only False supported.
    @return: A volume
    """
    if inputPath.endswith(".pdb"):
        atomList = naivePDBParser(inputPath, chain)
    elif inputPath.endswith(".cif"):
        atomList = mmCIFParser(inputPath, chain)
    else:
        raise RuntimeError("cifpdb2em : Unsupported file format! ", inputPath)

    # Default and only support toCompact : True. Automatically recenter.
    volm = atomList2emCube(atomList, pixelSize, False) # densityNegartive : False.
    
    if fname: # if fname is given, write to disk.
        volm.write(fname)
        print(f"cifpdb2em : MRC file {fname} is saved.")
    return volm

def mrc2em(filename,destname):
    from pytom.basic.files import read
    from pytom.tools.files import checkFileExists,checkDirExists
    if not checkFileExists(filename):
        raise RuntimeError('MRC file not found! ',filename)
    emfile = read(filename)
    emfile.write(destname,'em')

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
                try:
                    volumeData[i,j,k] = inputVolume.getV(i,j,k)
                except:
                    pass
    
    with mrcfile.new(mrcPath, overwrite=overwrite) as mrc:
        mrc.set_data(volumeData)
        print(f"mrc data dimension is converted to... {mrc.data.shape}")
    return