def subtomoSampleSaverCSV(tomoIdentifier, scenarioDir, subtomoIdentifier, subtomoDir, crowdLevel, SNR=1.0, generateNum=1, subtomoSizeX=50, subtomoSizeY=0, subtomoSizeZ=0):
    metadataCSV = f"{subtomoDir}/{subtomoIdentifier}_files.csv"
    csvFile = open(metadataCSV, "w")
    # Format : [subtomogram mrc file name],[subtomogram particle mask file name]
    for i in range(generateNum):
        mrcFileName = f"{subtomoIdentifier}_{i + 1}_particleMask.mrc"
        subtomoMRCFile = f"{subtomoIdentifier}_{i + 1}.mrc"
        csvFile.write(f"{subtomoMRCFile},{mrcFileName}\n")
    csvFile.close()