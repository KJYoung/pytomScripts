from pytom.basic.files import pdb2em, read
from pytom.simulation.EMSimulation import simpleSimulation
from pytom_volume import vol, initSphere
from pytom.basic.structures import WedgeInfo

min_Cube_Size = 125

modeVol = False
modeSim = True
# PDB to MRC
# for i in [100, 150, 1000, 1900]: for first settings.

if modeVol:
    for i in [125]:
        pdbPath = "/cdata/pdbData/1bxn.pdb"
        pdbPixelSize = 2.7
        volumePath = "/cdata/volumes/1bxnfalse_{}.mrc".format(i)
        pdb2em(pdbPath, chain=None, pixelSize=pdbPixelSize, cubeSize=i, fname=volumePath, densityNegative=False)
        volumePath = "/cdata/volumes/1bxntrue_{}.mrc".format(i)
        pdb2em(pdbPath, chain=None, pixelSize=pdbPixelSize, cubeSize=i, fname=volumePath, densityNegative=True)

def simulation(volumePath, simulatedPath):
    v = read(volumePath)
    wedge = 0.
    shift = [-1, 2, 3]
    rotation = [0, 0, 0]
    # there is a slight inconsistency when smoothing > 0 -
    # cleaner implementation would be multipliction with sqrt(mask) in corr function
    wi = WedgeInfo(wedgeAngle=wedge, rotation=[10.0,20.0,30.0], cutoffRadius=0.0)
    s = simpleSimulation( volume=v, rotation=rotation, shiftV=shift, wedgeInfo=wi, SNR=10.)
    s.write(simulatedPath)

def simulationSNR_Rotation_Wedge(volumePath, simulatedPath, snrValue, rotation=None, wedgeAngle=None, shift=None):
    wedge = 0.0
    shiftList = [0, 0, 0]

    v = read(volumePath)
    if rotation == None:
        rotation = [0, 0, 0]
    if not wedgeAngle == None:
        wedge = wedgeAngle
    if not shift == None:
        shiftList = shift
    
    wi = WedgeInfo(wedgeAngle=wedge, cutoffRadius=0.0)
    s = simpleSimulation( volume=v, rotation=rotation, shiftV=shiftList, wedgeInfo=wi, SNR=snrValue)
    s.write(simulatedPath)

def call_simulationSNR(ident):
    for snr in [1000.0, 100.0, 10.0, 1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
        volumePath = "/cdata/volumes/1bxntrue_{}.mrc".format(ident)
        simulatedPath = "/cdata/simulated/1bxn_true_{}_SNR_{}.mrc".format(ident,snr)
        simulationSNR(volumePath, simulatedPath, snr)
        volumePath = "/cdata/volumes/1bxnfalse_{}.mrc".format(ident)
        simulatedPath = "/cdata/simulated/1bxn_false_{}_SNR_{}.mrc".format(ident, snr)
        simulationSNR_Rotation_Wedge(volumePath, simulatedPath, snr)

def call_simulationRotation(ident):
    for snr in [100.0]:
        for axis in [0, 1, 2]:
            for angle in [0, 30, 45, 60, 90]:
                rotation = []
                if axis == 0:
                    rotation = [angle, 0, 0]
                elif axis == 1:
                    rotation = [0, angle, 0]
                else:
                    rotation = [0, 0, angle]
                
                volumePath = "/cdata/volumes/1bxnfalse_{}.mrc".format(ident)
                simulatedPath = "/cdata/simulated/1bxn_false_{}_SNR_{}_{}axis{}.mrc".format(ident, snr, axis, angle)
                simulationSNR_Rotation_Wedge(volumePath, simulatedPath, snr, rotation)

def call_simulationWedge(ident):
    for snr in [0.5]:
        for angle in [0, 50, 60, 70, 80]:
            rotation = [0,0,0]
            
            volumePath = "/cdata/volumes/1bxnfalse_{}.mrc".format(ident)
            simulatedPath = "/cdata/simulated/1bxn_false_{}_SNR_{}_wedge{}.mrc".format(ident, snr, angle)
            simulationSNR_Rotation_Wedge(volumePath, simulatedPath, snr, rotation, angle)

def call_simulationShift(ident):
    for snr in [0.7]:
        for axis in [0, 1, 2]:
            for shift in [0, 10, 30, 50, 90]:
                shiftA = []
                if axis == 0:
                    shiftA = [shift, 0, 0]
                elif axis == 1:
                    shiftA = [0, shift, 0]
                else:
                    shiftA = [0, 0, shift]
                
                volumePath = "/cdata/volumes/1bxnfalse_{}.mrc".format(ident)
                simulatedPath = "/cdata/simulated/1bxn_false_{}_SNR_{}_{}shift{}.mrc".format(ident, snr, axis, shift)
                simulationSNR_Rotation_Wedge(volumePath, simulatedPath, snr, shift=shiftA)

if modeSim:
    call_simulationShift(125)

