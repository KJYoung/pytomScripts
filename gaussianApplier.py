import argparse
import os
import numpy as np
import mrcfile
from os.path import exists
import torch
import torchvision

# Parsing arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Gaussian noise applier for MRC file")
    # parser.add_argument('-m', '--merge', dest='isMerge', action='store_true', help='Merge if flaged, Fragment otherwise')

    parser.add_argument('-i', '-d', '--data', required=True, help="Path to dataset")
    # parser.add_argument('-id', '--iden', type=str, default="", help="Fragment : output prefix. Merge : output name.")
    # argument for fragment.
    # parser.add_argument('-s', '--size', type=int, default=128, help="Fragment : output size.")
    # argument for merge.
    # parser.add_argument('-wi', '--width', type=int, default=32, help="Merge : input patch's w number.")
    # parser.add_argument('-he', '--height', type=int, default=32, help="Merge : input patch's h number.")
    # parser.add_argument('-is', '--fragsize', type=int, default=128, help="Merge : input patch's size.")
        # assume w size = h size
    return check_args(parser.parse_args())

# Checking arguments
def check_args(args):
    # --epoch
    withError = False
    if args.isMerge:
        # Merge Args chekc!
        try:
            assert (args.iden == "" or args.iden.endswith(".mrc"))
        except:
            withError = True
            print("Output file name should be endswith .mrc")
    else:
        # Fragmentation Args check!
        try:
            assert args.data.endswith('.mrc')
        except:
            withError = True
            print('Input file should be the mrc file')
        
        try:
            if not args.isMerge:
                assert exists(args.data)
        except:
            withError = True
            print('Input file does not exists : {}'.format(args.data))
        # --batch_size
        try:
            assert args.size >= 16
        except:
            withError = True
            print('Size must be larger than or equal to 16')
    return withError, args

#################### MAIN ####################
def applier_Adder(cleanNumpy, noiseNumpy):
    return cleanNumpy + noiseNumpy
    
if __name__ == '__main__':
    print("------------ Gaussian Noise Applier started... ------------")
    cleanExample = "/cdata/WrapperTEM/Micrographs/clean4/D4_1_1_p657_df3512-06_10_21.mrc"
    
    with mrcfile.open(cleanExample) as mrcPatch:
        mrcNumpy        = mrcPatch.data
        # print(mrcNumpy.shape) # (1, 4096, 4096)
        size            = mrcNumpy.shape[1]
    mrcAvg          = np.average(mrcNumpy)
    mrcStd          = np.std(mrcNumpy)
    print("input mrc file's avg {} and std {}".format(mrcAvg, mrcStd))
    #gaussianNoise = torch.normal(0,1, size=(size, size))
    gaussianNoise   = np.random.normal(size=(size, size))
    applied_patch   = applier_Adder(mrcNumpy, gaussianNoise).astype(np.float32)

    with mrcfile.new("test.mrc", overwrite=True) as mergedMRC:
        mergedMRC.set_data(applied_patch)
    torchvision.utils.save_image(torch.from_numpy(applied_patch), "test.png")
    print("saved to test.png")
    
    
    
    # assert False
    # withError, args = parse_args()
    # if withError:
    #     pass
    # else:
    #     #print("Output size will be : ", args.size)
    #     #print("Input file is : ", args.data)
    #     if args.isMerge:
    #         print("-- Merge mode.")
    #         mergeMRC(args)
    #     else:
    #         print("-- Fragment mode.")
    #         fragmentMRC(args.data, args.size, args.iden)
    
    # print("------------ Gaussian Noise Applier ended... ------------")


# python fragmentMRC.py -s 128 -d /cdata/WrapperTEM/Micrographs/noisy3/FIN1_particle644_df3733_minDist111.9403_Date05_09_18_28.mrc -id /cdata/temp/FIN1_half