from pytom.basic.files import pdb2em, mmCIF2em

# VKJY
import wget
import os.path
import re
import requests
from urllib.error import HTTPError

MIN_CUBE_SIZE_1bxn = 125

def wgetPDB2Volume(pdbID, cubeSize):
    # densityNegative for default
    pdbPath = "/cdata/pdbData/{}.pdb".format(pdbID)
    cifPath = "/cdata/pdbData/{}.cif".format(pdbID)
    
    pdbURL = "https://files.rcsb.org/view/{}.pdb".format(pdbID)
    cifURL = "https://files.rcsb.org/view/{}.cif".format(pdbID)
    
    volumePath = "/cdata/volumes/{}_{}.mrc".format(pdbID, cubeSize)
    
    if False: #os.path.isfile(pdbPath) or os.path.isfile(cifPath):
        print("pdb File with ID {} is already exists! Returned.".format(pdbID));
        return
    else:
        try:
            response = requests.get(pdbURL)
            wget.download(pdbURL, out=pdbPath)
            resPatternPDB = re.compile("RESOLUTION\..*([0-9]+\.[0-9]+).*ANGSTROMS")
            
            f = open(pdbPath, 'r')
            pdbContent = f.read()
            f.close()

            pdbPixResolution = re.findall(resPatternPDB, pdbContent)[0]
            print("pdb resolution of {} is {}".format(pdbID, pdbPixResolution))
            pdb2em(pdbPath, chain=None, pixelSize=float(pdbPixResolution), cubeSize=cubeSize, fname=volumePath, densityNegative=False)
            return
        except HTTPError as exception:
            try:
                response = requests.get(cifURL)
                wget.download(cifURL, out=cifPath)
                resPatternCIF = re.compile("_em_3d_reconstruction.resolution +([0-9]+\.[0-9]+)")

                f = open(cifPath, 'r')
                cifContent = f.read()
                f.close()

                cifPixResolution = re.findall(resPatternCIF, cifContent)[0]
                print("cif resolution of {} is {}".format(pdbID, cifPixResolution))
                volume = mmCIF2em(cifPath, chain=None, pixelSize=float(cifPixResolution), cubeSize=cubeSize, densityNegative=False)
                volume.write(volumePath)
            except HTTPError as exception:
                print("Invalid pdb ID maybe.")
                return

        
SHREC2021 = [ "1s3x", "3qm1", "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr", "4cr2", "5mrc" ]
# SHREC2021 = ["1s3x"]
for pdbID in SHREC2021:
    wgetPDB2Volume(pdbID, 200)

