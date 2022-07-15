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

# for test dataset 4.2 : same as v4. : 05/23.
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0523_compact10pV4", pixelSize=10, tomoSize=256, pfailedAttempts=10000, pparticleNum=3000, rotationStep=5, JSONCOMPACT=True, verbose=True)
subtomoSampleSaver("0523_compact10pV4", "/cdata/scenario", "0523_wonoise", "/cdata/subtomo", 0, SNR=-1, generateNum=100, subtomoSizeX=50)
subtomoSampleSaver("0523_compact10pV4", "/cdata/scenario", "0523_SNR_2.0", "/cdata/subtomo", 0, SNR=2.0, generateNum=100, subtomoSizeX=50)

# for test dataset 5. : same.
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd7(Max)", pixelSize=10.0, tomoSize=256, pfailedAttempts=200000, pparticleNum=15000, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd6", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=200000, pparticleNum=10000, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd5", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=100000, pparticleNum=5000, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd4", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=50000, pparticleNum=3000, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd3", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=30000, pparticleNum=2000, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd2", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=20000, pparticleNum=1000, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_crowd1(Min)", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=20000, pparticleNum=800, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_hypocrowd4", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=15000, pparticleNum=600, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_hypocrowd3", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=400, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_hypocrowd2", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=200, rotationStep=5, JSONCOMPACT=True, verbose=True)
makeGrandModelByPDBIDs(SHREC2021_FULL, "/cdata/pdbData", "/cdata/resolution4", "/cdata/scenario", "0405_hypocrowd1", newVolume=False, pixelSize=10.0, tomoSize=256, pfailedAttempts=10000, pparticleNum=100, rotationStep=5, JSONCOMPACT=True, verbose=True)
testSet5 = [ "0405_crowd7(Max)", "0405_crowd6", "0405_crowd5", "0405_crowd4", "0405_crowd3", "0405_crowd2", "0405_crowd1(Min)", "0405_hypocrowd4", "0405_hypocrowd3", "0405_hypocrowd2", "0405_hypocrowd1" ]
for test in testSet5:
    subtomoSampleSaverCSV(test, "/cdata/scenario/", f"{test}_noise2.0", "/cdata/subtomo0405", 0, SNR=2.0, generateNum=20, subtomoSizeX=50)


# for test dataset 7. : EMAN2 + various crowd level.
test dataset 7 scenario generation. -------------
particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523(PDB)_crowd6(Max)_EMAN2", pixelSize=10, tomoSize=256, pfailedAttempts=2000000, pparticleNum=9999999, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
print(particleNum)
particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523(PDB)_crowd5_EMAN2", pixelSize=10, tomoSize=256, pfailedAttempts=300000, pparticleNum=9999999, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd4", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=10000, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd3", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=4000, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd2", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=1500, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
particleNum = makeGrandModelByPDBIDs(SHREC2021_PDB, "/cdata/pdbData", "/cdata/emByEMAN2", "/cdata/scenario", "0523_crowd1", pixelSize=10, tomoSize=256, pfailedAttempts=250000, pparticleNum=600, rotationStep=5, pdb2em=EMAN2, JSONCOMPACT=True, verbose=True)
volume2MRC("/cdata/scenario/0523(PDB)_crowd6(Max)_EMAN2.em", "/cdata/scenario/0523(PDB)_crowd6(Max)_EMAN2.mrc", floatMRC=True, overwrite=False, verbose=True)
volume2MRC("/cdata/scenario/0523(PDB)_crowd5(Max)_EMAN2.em", "/cdata/scenario/0523(PDB)_crowd5(Max)_EMAN2.mrc", floatMRC=True, overwrite=False, verbose=True)
volume2MRC("/cdata/scenario/0523_crowd4.em", "/cdata/scenario/0523_crowd4.mrc", floatMRC=True, overwrite=False, verbose=True)
volume2MRC("/cdata/scenario/0523_crowd3.em", "/cdata/scenario/0523_crowd3.mrc", floatMRC=True, overwrite=False, verbose=True)
volume2MRC("/cdata/scenario/0523_crowd2.em", "/cdata/scenario/0523_crowd2.mrc", floatMRC=True, overwrite=False, verbose=True)
volume2MRC("/cdata/scenario/0523_crowd1.em", "/cdata/scenario/0523_crowd1.mrc", floatMRC=True, overwrite=False, verbose=True)