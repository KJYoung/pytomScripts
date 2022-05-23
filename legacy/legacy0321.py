def subtomoSampler(identifier, scenarioDir, crowdLevel, generateNum=1, subtomoSizeX=50, subtomoSizeY=0, subtomoSizeZ=0):
    # Missing size is filled with X axis.
    if subtomoSizeY == 0:
        subtomoSizeY = subtomoSizeX
    if subtomoSizeZ == 0:
        subtomoSizeZ = subtomoSizeX
    
    scenarioVolume = read(f'{scenarioDir}/{identifier}.em')
    
    particleCenters = []
    # Load txt file
    with open(f'{scenarioDir}/{identifier}.txt') as scenarioParticleTxt:
        txt_contents = scenarioParticleTxt.readlines()
        particleTxtPattern = re.compile("(.*),(.*),(.*),(.*),(.*),(.*),(.*)")
        particleID = 0
        for line in txt_contents:
            parsedInfo = re.findall(particleTxtPattern, line)[0]
            particleCenter = [ parsedInfo[0], parsedInfo[1], parsedInfo[2], parsedInfo[3], particleID ]
            particleID += 1
            particleCenters.append(particleCenter)
    
    particleJsonList = []
    # Load json file
    with open(f'{scenarioDir}/{identifier}.json') as scenarioJsonFile:
        json_object = json.load(scenarioJsonFile)
        particleJsonList = json_object['particles']
        pdbIDList = json_object['pdbIDs']
        
    subtomos = []
    jsons = []
    particleList = []

    index = 0
    for i in range(generateNum):
        jsonMetadataObject = getMedadataJsonTemplate()
        jsonMetadataObject["header"] = f"resolution : {10.0} with white noise"
        jsonMetadataObject["pdbIDs"] = pdbIDList
        jsonMetadataObject["resolutions"] = [5.0]
        jsonMetadataObject["particles"] = []

        subtomo = vol(subtomoSizeX, subtomoSizeY, subtomoSizeZ)
        subtomo.setAll(0.0)

        if True: # type1 : just pick the center of the particle.
            # np.random.randint :: low exclusive high inclusive
            particleInfo = particleCenters[ np.random.randint(low = 0, high = len(particleCenters)) ]
            cX, cY, cZ = int(float(particleInfo[1])), int(float(particleInfo[2])), int(float(particleInfo[3]))
        
        if cX <= subtomoSizeX/2:
            lrX = 0
        elif cX >= scenarioVolume.sizeX() - subtomoSizeX/2:
            lrX = scenarioVolume.sizeX() - subtomoSizeX
        else:
            lrX = cX - subtomoSizeX//2
        if cY <= subtomoSizeY/2:
            lrY = 0
        elif cY >= scenarioVolume.sizeY() - subtomoSizeY/2:
            lrY = scenarioVolume.sizeY() - subtomoSizeY
        else:
            lrY = cY - subtomoSizeY//2
        if cZ <= subtomoSizeZ/2:
            lrZ = 0
        elif cZ >= scenarioVolume.sizeZ() - subtomoSizeZ/2:
            lrZ = scenarioVolume.sizeZ() - subtomoSizeZ
        else:
            lrZ = cZ - subtomoSizeZ//2
        
        for x in range(subtomoSizeX):
            for y in range(subtomoSizeY):
                for z in range(subtomoSizeZ):
                    getValue = scenarioVolume.getV( lrX + x , lrY + y, lrZ + z )
                    subtomo.setV( getValue , x, y, z)

                    if getValue != 0:
                        for particle in particleJsonList:
                            minList = particle['min']
                            maxList = particle['max']
                            curCord = [ lrX + x, lrY + y, lrZ + z]

                            if checkOverlap(minList, maxList, curCord):
                                # This Particle!
                                #print("index(particle)", particleJsonList.index(particle))
                                #print("min --------- \n", minList)
                                #print("max --------- \n", maxList)
                                #print("cur --------- \n", curCord)
                                particleList.append( [  particleJsonList.index(particle), [x,y,z] ] )
                                break
                            if particle == particleJsonList[-1]:
                                #print("Noise!")
                                pass
        particleKeys = []
        particleDicts = []
        particleIDList = []

        for p in particleList:
            if p[0] in particleKeys:
                for pp in particleDicts:
                    if pp['particleID'] == p[0]:
                        pp['coord'].append(p[1])
                        break
            else:
                particleKeys.append( p[0] )
                particleDicts.append( { "particleID" : p[0], "min" : [] , "max" : [] ,  "coord": [ p[1] ] } )
        
        for p in particleDicts:
            # particleIDList.append(p['classNum'])
            p['min'], p['max'] = findMinMaxByList(p['coord'])
            p['classNum'] = particleJsonList[  p['particleID'] ]['classNum']
            # del p['particleID']

        #print(particleDicts)
        jsonMetadataObject["particles"] = particleDicts
        subtomos.append(subtomo)
        jsons.append(jsonMetadataObject)
        print(f"--------- subtomogram generated : {index}")
        index += 1
    return subtomos, jsons