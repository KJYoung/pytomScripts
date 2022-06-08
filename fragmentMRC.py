import argparse
import os
import numpy as np
import mrcfile
from os.path import exists

# Parsing arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Fragmentation of MRC file")

    parser.add_argument('-s', '--size', type=int, default=128, help="output size.")
    parser.add_argument('-i', '-d', '--data', required=True, help='path to dataset')
    parser.add_argument('-id', '--iden', type=str, default="", help='output prefix')
    return check_args(parser.parse_args())

# Checking arguments
def check_args(args):
    # --epoch
    withError = False
    try:
        assert args.data.endswith('.mrc')
    except:
        withError = True
        print('Input file should be the mrc file')
    
    try:
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
def fragmentMRC(inputFile, size, ident):
    with mrcfile.open(inputFile) as mrc:
        print("Input file shape is : ", mrc.data.shape)
        wfactor, hfactor = (mrc.data.shape[1] // size), (mrc.data.shape[2] // size)
        print("Output patches will be generated : {} x {} patches".format(wfactor, hfactor))
        if ident == "":
            identifier = inputFile[:len(inputFile) - 4]
        else:
            identifier = ident
        for w in range(wfactor):
            for h in range(hfactor):
                patchName = "{}_w{}_h{}.mrc".format(identifier, w , h)
                with mrcfile.new(patchName) as mrcPat:
                    mrcPat.set_data(mrc.data[:, w * size : (w+1) * size , h * size : (h+1) * size ])
                    # print("  patch Name : {}".format(patchName))
                    # print("    patches : ({},{}) ~ ({},{}) ".format(w * size, h * size, w * size + size - 1, h * size + size -1))

if __name__ == '__main__':
    print("------------ fragmentMRC started... ------------")
    withError, args = parse_args()
    if withError:
        pass
    else:
        #print("Output size will be : ", args.size)
        #print("Input file is : ", args.data)
        fragmentMRC(args.data, args.size, args.iden)
    print("------------ fragmentMRC ended... ------------")


# python fragmentMRC.py -s 128 -d /cdata/WrapperTEM/Micrographs/noisy3/FIN1_particle644_df3733_minDist111.9403_Date05_09_18_28.mrc -id /cdata/temp/FIN1_half