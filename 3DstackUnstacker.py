import sys,os,mrcfile

srcDir    = sys.argv[1]  # micrographs
dstDir    = sys.argv[2]  # raw micrographs

input_list= os.listdir(srcDir)

for i, srcName in enumerate(input_list):
    if srcName.endswith(".mrc"):
        with mrcfile.open(srcDir + srcName) as mrc:
            npdata = mrc.data
            assert len(npdata.shape) == 3
            if npdata.shape[0] < npdata.shape[2]:
                # 0 unstack
                unstackNum  = npdata.shape[0]
                unstackSize = (npdata.shape[1], npdata.shape[2])
                for j in range(unstackNum):
                    dstName = srcName[:-4] + "_{}".format(j) + ".mrc"
                    with mrcfile.new(dstDir + dstName) as mrcUnstack:
                        mrcUnstack.set_data(npdata[j, :, :])
            elif npdata.shape[0] > npdata.shape[0]:
                # 2 unstack
                unstackNum = npdata.shape[2]
                unstackSize = (npdata.shape[0], npdata.shape[1]) 
                for j in range(unstackNum):
                    dstName = srcName[:-4] + "_{}".format(j) + ".mrc"
                    with mrcfile.new(dstDir + dstName) as mrcUnstack:
                        mrcUnstack.set_data(npdata[:, :, j])
            else:
                print("what I have to unstack? {}".format(npdata.shape))
                quit()
        print('# unstacker : {} -> {} number of {} [{}/{}] {:.2%}'.format(npdata.shape, unstackNum, unstackSize,
                                                                          i, len(input_list), i/len(input_list)), file=sys.stderr, end='\r')
print("All jobs done.")

# python /cdata/workspace/3DstackUnstacker.py /cdata/EMPIAR-10164_ORIGINAL/frames/ /cdata/EMPIAR-10164_ORIGINAL/framesEach/