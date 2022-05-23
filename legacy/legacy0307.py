# Mainly, this file is with test functions for rotation scheme.
# First Version of naive Simulator : without rotation.
def makeScenarioByPDBIDs(pdbIDList, scenarioPATH, occList=None, cubeSize=256, overwrite=False, verbose=False):
    # Prepare Volumes.
    # cuboidalOccupancyList = [['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]]]
    if occList:
        cuboidalOccupancyList = occList
    else:
        cuboidalOccupancyList = occupancyCalculator(pdbIDList, overwrite, verbose)
    if verbose:
        print("----------- occupancy List is done... -----------")
        #print(cuboidalOccupancyList)
    failedAttempts = 0
    particleNum = 0
    scenario = []
    while failedAttempts < 1600000 and particleNum < 1600:
        if verbose and failedAttempts%400000 == 0:
            print(f"... Now failed Attemps are {failedAttempts}")   
        classNum = random.randint(0, len(cuboidalOccupancyList)-1)
        classOccpy = cuboidalOccupancyList[classNum]
        x, y, z = random.randint(0, cubeSize-1-classOccpy[1][0]), random.randint(0, cubeSize-1-classOccpy[1][1]), random.randint(0, cubeSize-1-classOccpy[1][2])
        
        if len(scenario) == 0:
            scenario.append([classNum, [x,y,z], classOccpy[1]])
            particleNum+=1
        else:
            for existingItem in scenario:
                if (x <= existingItem[1][0]+existingItem[2][0] and x+classOccpy[1][0] >= existingItem[1][0]) and \
                    (y <= existingItem[1][1]+existingItem[2][1] and y+classOccpy[1][1] >= existingItem[1][1]) and \
                    (z <= existingItem[1][2]+existingItem[2][2] and z+classOccpy[1][2] >= existingItem[1][2]):
                    failedAttempts+=1
                    break; 
                if existingItem == scenario[-1]:
                    scenario.append([classNum, [x, y, z], classOccpy[1]])
                    particleNum+=1
    
    print(f"----------- Scenario generation is done... with Particle Number {particleNum}----------")
    # print(scenario)
    # [[0, [305, 116, 185], [46, 32, 38]], [3, [404, 398, 234], [34, 69, 67]], [3, [158, 219, 313], [34, 69, 67]], [4, [31, 471, 392], [31, 36, 44]], [7, [215, 168, 1], [45, 41, 54]], [1, [445, 347, 63], [39, 32, 37]], [2, [75, 46, 145], [41, 34, 27]], [2, [361, 299, 102], [41, 34, 27]], [5, [267, 6, 446], [25, 36, 21]], [1, [118, 81, 469], [39, 32, 37]]]

    volume = vol(cubeSize, cubeSize, cubeSize)
    volume.setAll(0.0)

    for items in scenario:
        particleVol = read(f"/cdata/volumesCompact/{cuboidalOccupancyList[items[0]][0]}.em")
        x, y, z = particleVol.sizeX(), particleVol.sizeY(), particleVol.sizeZ()
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    if volume.getV(items[1][0]+i, items[1][1]+j, items[1][2]+k) != 0:
                        print("ERROR : PARTICLE OVERAPPED!")
                        print(f"debug information : {items}")
                    volume.setV( particleVol.getV(i,j,k), items[1][0]+i, items[1][1]+j, items[1][2]+k)
        
        # with mrcfile.new(mrcPath, overwrite=overwrite) as mrc:
        #     mrc.set_data(volumeData)
        #     print(f"mrc data dimension is converted to... {mrc.data.shape}")
        # return
        # COPY
    
    volume.write(scenarioPATH)
    #END

def occupancyCalculator(pdbIDList, pdbDir, volumeDir, mrcDir, overwrite=False, verbose=False):
    cuboidalOccupancyList = []
    for pdbID in pdbIDList:
        _vol, x, y, z = wgetPDB2Volume(pdbID, pdbDir=pdbDir, volumeDir=volumeDir, mrcDir=mrcDir, toCompact=True, generateMRC=False, overwrite=overwrite, verbose=verbose)
        cuboidalOccupancyList.append([pdbID, [x, y, z]])
        if verbose:
            print(f"Volume building with {pdbID} is done...")
    return cuboidalOccupancyList
# [['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]

def syntheticGen1(outputPath):
    volume = vol(64, 64, 64)
    volume.setAll(0.0)
    for i in range(64):
        volume.setV(20, 32, 32, i)

        volume.setV(20, i, 32, 32)
        volume.setV(20, i, 31, 31)
        volume.setV(20, i, 33, 33)

        volume.setV(20, 32, i, 32)
        volume.setV(20, 29, i, 29)
        volume.setV(20, 35, i, 35)
    
    mat = np.zeros([64, 64])
    for i in range(64):
        mat[i][i] = 20
    matrixToVolLayerX(volume, mat, 40)
    volume.write(outputPath)

def getWeightedSum(inputVolume):
    x, y, z = inputVolume.sizeX(), inputVolume.sizeY(), inputVolume.sizeZ()
    sumOfPoints = 0
    for i in range(x):
        for j in range(y):
            for k in range(z):
                sumOfPoints += inputVolume.getV(i,j,k)    
    return sumOfPoints

def rotationIntegrityTest2(volumePath):
    em = read(volumePath)
    angles = range(0, 360, 20)
    f = open("/cdata/0305/rotateIntegrity_cube.txt", 'w')
    for angleX in angles:
        print(f"angleX : {angleX} is started...")
        for angleY in angles:
            print(f"--angleY : {angleY} is started...")
            for angleZ in angles:
                from pytom_volume import rotate
                rotatedCopy = vol(em.sizeX(),em.sizeY(),em.sizeZ())
                rotatedCopy.setAll(0.0)
                rotate(em, rotatedCopy, angleX, angleY, angleZ)
                compactRotated = makeCompact(rotatedCopy)
                result = getWeightedSum(compactRotated)
                f.write(f"{angleX},{angleY},{angleZ},{compactRotated.sizeX()},{compactRotated.sizeY()},{compactRotated.sizeZ()},{result}\n")
    f.close()

def rotationIntegrityTest1(volumePath): #1bxn 45 45 37
    em = read(volumePath)
    #print(em.sizeX(),em.sizeY(),em.sizeZ())
    angles = range(0, 360, 20)
    f = open("/cdata/0305/rotateIntegrity_cuboid.txt", 'w')
    for angleX in angles:
        print(f"angleX : {angleX} is started...")
        for angleY in angles:
            print(f"--angleY : {angleY} is started...")
            for angleZ in angles:
                from pytom_volume import rotate
                rotatedCopy = vol(em.sizeX(),em.sizeY(),em.sizeZ())
                rotatedCopy.setAll(0.0)
                rotate(em, rotatedCopy, angleX, angleY, angleZ)
                result = getWeightedSum(rotatedCopy)
                f.write(f"{angleX},{angleY},{angleZ},{result}\n")
    f.close()

def randomRotationAndCompact(pdbID, volumeDir, rotationDir, overwrite=False, iter=1):
    rotList = []
    for i in range(iter):
        phi, theta, psi = np.random.randint(low=0, high=360, size=(3,))
        print(type(phi), phi, theta, psi)
        volumePath = f"{volumeDir}/{pdbID}.em"
        rotPath = f"{rotationDir}/{pdbID}_{phi}_{theta}_{psi}.em"
        comPath = f"{rotationDir}/compact/{pdbID}_{phi}_{theta}_{psi}.em"
        em = read(volumePath)
        from pytom_volume import rotate
        rotatedCopy = vol(em.sizeX(),em.sizeY(),em.sizeZ())
        rotatedCopy.setAll(0.0)
        rotate(em, rotatedCopy, int(phi), int(theta), int(psi))
        rotatedCopy.write(rotPath)
        comp = makeCompact(rotatedCopy)
        comp.write(comPath)
        rotList.append(comp)
    return rotList



    # 01 # volList = prepareCubeVolumes(SHREC2021_FULL, pdbDir=pdbDataDIR, volumeDir=singleParticleEMCubeDIR, mrcDir=None, toCompact=True, overwrite=True, verbose=True)
    ###### get Cube volume(.em) files by PDB ID List.
    #def makeScenarioByPDBIDs(pdbIDList, volumeDir, scenarioDir, scenarioIdentifier="noname", cubeSize=128, overwrite=False, verbose=False):
    
    #volume2MRCconverter('/cdata/volumesV/1bxn_200.mrc', '/cdata/mrcs/1bxn_200.mrc', verbose=True)
    #wgetPDB2Volume('1bxn', toCompact=True, generateMRC=True, verbose=True)
    # Rotation integrity Test
    #wgetPDB2Volume('1bxn', pdbDir=pdbDataDIR, volumeDir=singleParticleEMCuboidDIR, mrcDir=None, toCompact=True, overwrite=True, verbose=True)
    #compactCuboid2rotateCube("/cdata/singleParticleEM_cuboid/1bxn.em", "/cdata/singleParticleEM_cube/1bxn.em")
    

    #f.write(f"first job start : {executionStart}\n")
    #makeScenarioByPDBIDs(["1bxn"], volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="only1bxn2", cubeSize=128)
    #f.write(f"first job end & second start : {firstEnd - executionStart}\n")
    #makeScenarioByPDBIDs(SHREC2021_FULL, volumeDir=singleParticleEMCubeDIR, scenarioDir=scenarioDIR, scenarioIdentifier="fulltest3_reducedfail", cubeSize=256, verbose=True)
    ##v = read("/cdata/scenario/fulltest3_reducedfail.em")
    #volumeOutliner(v).write("/cdata/scenario/outlined_fulltest3_reducedfail.em")
    #f.write(f"second finished : {time.time() - firstEnd}\n")
    #volList = randomRotationAndCompact("1bxn", volumeDir=singleParticleEMCubeDIR, rotationDir="/cdata/0307/temp_rotation_test", iter=10)
    #f.close()
    # vol, _x, _y, _z = wgetPDB2Volume('1bxn', pdbDir=pdbDataDIR, volumeDir=singleParticleEMCubeDIR, mrcDir=None, toCompact=True, overwrite=True, verbose=True)
    # i = 0
    # for v in volList:
    #     volumeOutliner(v).write(f"/cdata/outlined/rotation_outlined_{i}.em")
    #     i+=1
    #rotationIntegrityTest2("/cdata/singleParticleEM_cube/1bxn.em")

    #syntheticGen1(f'{dateDir}/rotest4X.em')
    #customSimulation(f'{dateDir}/rotest4X.em', f'{dateDir}/rotest4X_0_0_30.em', 100, rotation=[0,0,30])
    #customSimulation(f'{dateDir}/rotest4X.em', f'{dateDir}/rotest4X_0_30_0.em', 100, rotation=[0,30,0])
    #customSimulation(f'{dateDir}/rotest4X.em', f'{dateDir}/rotest4X_30_0_0.em', 100, rotation=[30,0,0])


    # makeScenarioByPDBIDs(SHREC2021_FULL, scenarioPATH=scenarioFile, occList=SHREC2021_FULL_OCCLIST, cubeSize=256,overwrite=True, verbose=True)
    # customSimulation(scenarioFile, simulatedFile, 0.1, wedgeAngle=70)