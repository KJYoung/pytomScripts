from pytom.basic.files import recenterVolume, naivePDBParser, mmCIFParser

def atomList2emMulti(atomList, pixelSize, cubeSize, densityNegative=False):
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


def pdb2emMulti(pdbPath, pixelSize, cubeSize, chain=None, densityNegative=False, fname='', recenter=True):
    """
    pdb2em: Creates an volume out of a PDB file
    @param pdbPath: Path to PDB file or PDB id for online download
    @param pixelSize: The pixel size to convert to 
    @param cubeSize: Resulting cube size
    @return: A volume
    @author: Thomas Hrabe & Luis Kuhn
    """
    from math import floor
    from pytom_volume import vol

    atomList = naivePDBParser(pdbPath, chain)

    vol = atomList2emMulti(atomList, pixelSize, cubeSize, densityNegative)

    if recenter:
        vol = recenterVolume(vol, densityNegative)

    if fname:
        vol.write(fname)

    else:
        return vol

def mmCIF2emMulti(mmCIFPath, pixelSize, cubeSize, chain=None, densityNegative=False, fname='', recenter=True):
    """
    pdb2em: Creates an volume out of a mmCIF file
    @param mmCIFPath: Path to mmCIF file 
    @param pixelSize: The pixel size to convert to 
    @param cubeSize: Resulting cube size
    @return: A volume
    @author: Thomas Hrabe 
    """
    from math import floor
    from pytom_volume import vol

    atomList = mmCIFParser(mmCIFPath, chain)

    vol = atomList2emMulti(atomList, pixelSize, cubeSize, densityNegative)

    if recenter:
        vol = recenterVolume(vol, densityNegative)

    if fname:
        vol.write(fname)
        
    return vol

# VKJY
import wget
import os.path
import re
import requests
from urllib.error import HTTPError

MIN_CUBE_SIZE_1bxn = 125

def wgetPDB2Volume(pdbID, cubeSize):
    # densityNegative for default
    print("wgetPDB2Volume is working with PDBID : {} and cubeSize : {}".format(pdbID, cubeSize))
    pdbPath = "/cdata/pdbData/{}.pdb".format(pdbID)
    cifPath = "/cdata/pdbData/{}.cif".format(pdbID)
    
    pdbURL = "https://files.rcsb.org/view/{}.pdb".format(pdbID)
    cifURL = "https://files.rcsb.org/view/{}.cif".format(pdbID)
    
    volumePath = "/cdata/volumesV/{}_{}.mrc".format(pdbID, cubeSize)
    
    templateExists = False

    if os.path.isfile(pdbPath) or os.path.isfile(cifPath):
        # print("pdb File with ID {} is already exists! Returned.".format(pdbID));
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
        volume = mmCIF2emMulti(cifPath, chain=None, pixelSize=float(cifPixResolution), cubeSize=cubeSize, densityNegative=False)
        volume.write(volumePath)
        return
    if not templateExists:
        wget.download(pdbURL, out=pdbPath)
    resPatternPDB = re.compile("RESOLUTION\..*([0-9]+\.[0-9]+).*ANGSTROMS")
    
    f = open(pdbPath, 'r')
    pdbContent = f.read()
    f.close()

    pdbPixResolution = re.findall(resPatternPDB, pdbContent)[0]
    print("pdb resolution of {} is {}".format(pdbID, pdbPixResolution))
    pdb2emMulti(pdbPath, chain=None, pixelSize=float(pdbPixResolution), cubeSize=cubeSize, fname=volumePath, densityNegative=False)
    return
        
SHREC2021_FULL = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2", "5mrc" ]
SHREC2021_1bxn = ["1bxn"]

for pdbID in SHREC2021_1bxn:
    wgetPDB2Volume(pdbID, 240)

# Auto incremental system.
# for pdbID in SHREC2021:
#     cubeSize = 100
#     while True:
#         try:
#             wgetPDB2Volume(pdbID, cubeSize)
#             break
#         except RuntimeError:
#             cubeSize+=10

