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
    print("------------ Grid Visualizer started... ------------")
    extractNoiseDir = "/cdata/db1/noiseExtract/04/noisePatch/"
    input_list=os.listdir(extractNoiseDir)

    container = torch.zeros(64, 1, 320, 320)
    for img in range(10):
        for i in range(64):
            with mrcfile.open(extractNoiseDir + input_list[i + img * 64]) as extNoise:
                numpyData = np.copy(extNoise.data)
                container[i, 0, :, :] = torch.Tensor(numpyData)
        grid = torchvision.utils.make_grid(container, normalize=True)
        torchvision.utils.save_image(grid, 'grid_{}.png'.format(img))
    print("------------ Grid Visualizer ended... ------------")