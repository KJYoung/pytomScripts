SHREC2021_FULL_OCCLIST = [['1s3x', [34, 31, 28]], ['3qm1', [23, 32, 22]], ['3gl1', [46, 32, 38]], ['3h84', [39, 32, 37]], ['2cg9', [41, 34, 27]], ['3d2f', [34, 69, 67]], ['1u6g', [31, 36, 44]], ['3cf3', [25, 36, 21]], ['1bxn', [44, 44, 36]], ['1qvr', [45, 41, 54]], ['4cr2', [36, 27, 27]], ['5mrc', [95, 86, 74]]]
SHREC2021_FULLexc2L2S = [ "3gl1", "3h84", "2cg9", "3d2f", "1u6g", "3cf3", "1bxn", "1qvr" ]
SHREC2021_1bxn = [ "1bxn" ]

ov = volumeOutliner("/cdata/scenario/0315_merge2.em", isFile=True, outlineValue = 500)
ov.write("/cdata/outlined/0315_merge2_utilstest.em")
# for test dataset 2.
# for - grandmodel.
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0321_gmwoN_5.0compact", 10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=2200, rotationStep=2, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0321_gmwoN_5.0verbose", 10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=2200, rotationStep=2, JSONCOMPACT=False, verbose=True)
# for - subtomogram.
# previous approach.
sl, js = subtomoSampler("0321_gmwoN_5.0compact", "/cdata/scenario", 2, generateNum=3, subtomoSizeX=50)
volumeListWriter(sl, "/cdata/scenario", "0318_2_6gmwoN_5.0compact_subtomo", JSON=js)
# unified approach.
subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0321_gmwoN_5.0compact_1", "/cdata/subtomo", 0, SNR=1.0, generateNum=10, subtomoSizeX=50)
subtomoSampleSaver("0321_gmwoN_5.0verbose", "/cdata/scenario", "0321_gmwoN_5.0verbose_1", "/cdata/subtomo", 0, SNR=1.0, generateNum=10, subtomoSizeX=50)

# for test dataset 3. : mrc output / mrc metadata
subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0328_gmwoN_5.0wonoise_4", "/cdata/subtomo", 0, SNR=-1, generateNum=15, subtomoSizeX=50)
subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0328_gmwoN_5.0noise2.0_4", "/cdata/subtomo", 0, SNR=2.0, generateNum=15, subtomoSizeX=50)
# for test dataset 3.1 : merged mrc metadata.
subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0328_5.0wonoise_5", "/cdata/subtomo", 0, SNR=-1, generateNum=15, subtomoSizeX=50)
subtomoSampleSaver("0321_gmwoN_5.0compact", "/cdata/scenario", "0328_5.0noise2.0_5", "/cdata/subtomo", 0, SNR=2.0, generateNum=15, subtomoSizeX=50)

# for test dataset 4. : subtomo occupy percentage threshold.
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0328_compact10pV", pixelSize=10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=3000, rotationStep=5, JSONCOMPACT=True, verbose=True)
subtomoSampleSaver("0328_compact10pV", "/cdata/scenario", "0329_wonoise_5", "/cdata/subtomo", 0, SNR=-1, generateNum=15, subtomoSizeX=50)
subtomoSampleSaver("0328_compact10pV", "/cdata/scenario", "0329_noise2.0_5", "/cdata/subtomo", 0, SNR=2.0, generateNum=15, subtomoSizeX=50)

