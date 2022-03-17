# TODO
    from scipy.ndimage import center_of_mass
    from pytom.tompy.io import read, write
    from pytom.tompy.tools import paste_in_center
    from pytom.gpu.initialize import xp
    from pytom_numpy import vol2npy
    from pytom.tools.files import checkFileExists
    from pytom_volume import read







def read(file, subregion=[0, 0, 0, 0, 0, 0], sampling=[0, 0, 0], binning=[0, 0, 0]):
    """
    read: Reads a file
    @param file: Path to file. Supports EM, MRC and CCP4 files 
    @type file: str
    @param subregion: Specifies a subregion to be read. The first three 
    values specify the upper left coordinates within the large volume, 
    the last three the length of the subregion along each dimension. 
    @type subregion: List of 6 integers  
    @param sampling: Change read sampling. Read every second (2), 
    third (3) pixel along each dimension.
    @type sampling:  List of 3 integers
    @param binning: Bin volume along each dimension. Note, 1 will do nothing,
    2 will bin with a kernelsize of 2 pixels along each dimension, 3 will bin
    with a kernelsize of 3 pixels along each dimension and so forth. 
    @type binning:  List of 3 integers
    @return: A volume object. L{pytom_volume.vol}
    @author: Thomas Hrabe
    """
    from pytom.tools.files import checkFileExists
    from pytom_volume import read

    if not file.__class__ == str:
        raise TypeError('File parameter must be a string!')

    if not checkFileExists(file):
        raise IOError('File not found or path is wrong: ' + file)

    try:
        f = read(file, subregion[0], subregion[1], subregion[2], subregion[3],
                 subregion[4], subregion[5], sampling[0], sampling[1], sampling[2],
                 binning[0], binning[1], binning[2])
        return f
    except (RuntimeError, errorNumber, errorString):
        # redundant to code above, but just in case it goes through
        if "Wrong file format or file doesn't exist!" in errorString:
            raise IOError('File not found or path is wrong: ' + file)
        else:
            raise

class NaiveAtom:

    def __init__(self, atomSeq, atomType, x, y, z, resSeq, resType):

        self._atomSeq = atomSeq
        self._atomType = atomType
        self._x = x
        self._y = y
        self._z = z
        self._resSeq = resSeq
        self._resType = resType

    def getAtomType(self):
        return self._atomType

    def getAtomSeq(self):
        return self._atomSeq

    def getX(self):
        return self._x

    def getY(self):
        return self._y

    def getZ(self):
        return self._z

    def setX(self, value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setX()')

        self._x = value

    def setY(self, value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setY()')

        self._y = value

    def setZ(self, value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setZ()')

        self._z = value

def naivePDBParser(pdbPath, chainName=None):
    atomList = []

    pdbFile = open(pdbPath, 'r')
    try:
        for line in pdbFile:

            name = line[:6]

            if name == "ATOM  ":
                chain = line[21]
                if chainName != 'all' and chainName != 'All':
                    if chainName is not None and chain != chainName:
                        continue
                atomdata = line.split()
                if len(atomdata) > 17:
                    line = atomdata
                    atomSeq  = line[1]
                    atomType = line[3]
                    resType  = line[5]
                    resSeq   = line[8]
                    x        = float(line[10])
                    y        = float(line[11])
                    z        = float(line[12])

                else:
                    atomSeq  = line[5:11]
                    atomType = line[11:17].strip()
                    resType  = line[17:20]
                    resSeq   = line[22:26]
                    x        = float(line[29:38])
                    y        = float(line[38:46])
                    z        = float(line[46:54])

                atomList.append(NaiveAtom(atomSeq, atomType, x, y, z, resSeq, resType))
    finally:
        pdbFile.close()

    return atomList

def mmCIFParser(filePath, chainName=None):
    """
    mmCIFParser: Parses mmCIF files from PDB and returns a list of Atom coordinates
    @param filePath: Path to the file
    @param chainName: Focus on one specific chain. Optional, if not specified, the whole file will be used. (This is the oposite to L{pytom.files.naivePDBParser}).
    @return: List of L{pytom.files.NaiveAtom}s. 
    """
    import re

    mmCIFFile = open(filePath, 'r')
    lines = mmCIFFile.readlines()
    mmCIFFile.close()

    atoms = []

    for line in lines:
        try:
            if line[:4] == 'ATOM':
                parts = re.sub('\s+', ' ', line).split(' ')
                chain = parts[6]

                if chainName is not None and chain != chainName:
                    continue

                atomType = parts[2].strip()
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                atom = NaiveAtom('', atomType, x, y, z, '', '')
                atoms.append(atom)

        except:
            continue
        finally:
            pass

    return atoms

def recenterVolume(volume, densityNegative=False):
    from scipy.ndimage import center_of_mass
    from pytom.tompy.io import read, write
    from pytom.tompy.tools import paste_in_center
    from pytom.gpu.initialize import xp
    from pytom_numpy import vol2npy
    import os

    try:
        a = vol2npy(volume).copy()
        vol =True
    except:
        a = volume
        vol = False

    if densityNegative:
        a *= -1

    x, y, z = list(map(int, center_of_mass(a)))
    cx, cy, cz = a.shape[0] // 2, a.shape[1] // 2, a.shape[2] // 2

    sx = min(x, a.shape[0] - x)
    sy = min(y, a.shape[0] - y)
    sz = min(z, a.shape[0] - z)

    ac = a[x - sx:x + sx, y - sy:y + sy, z - sz:z + sz]
    b = xp.zeros_like(a)

    b = paste_in_center(ac, b)

    if densityNegative: b *= -1

    if vol:
        write('recenteredDBV21.em', b)
        from pytom.basic.files import read
        vol = read('recenteredDBV21.em')
        os.system('rm recenteredDBV21.em')
        return vol
    else:
        return b



