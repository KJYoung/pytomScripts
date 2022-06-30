import argparse
import os
import numpy as np
import mrcfile
from os.path import exists

# Parsing arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Fragmentation of MRC file")

    parser.add_argument('-m', '--merge', dest='isMerge', action='store_true', help='Merge if flaged, Fragment otherwise')

    parser.add_argument('-i', '-d', '--data', required=True, help="Path to dataset.")
    parser.add_argument('-id', '--iden', type=str, default="", help="Output directory.")
    # argument for fragment.
    parser.add_argument('-s', '--size', type=int, default=128, help="Fragment : output size.")
    # argument for merge.
    parser.add_argument('-wi', '--width', type=int, default=32, help="Merge : input patch's w number.")
    parser.add_argument('-he', '--height', type=int, default=32, help="Merge : input patch's h number.")
    parser.add_argument('-is', '--fragsize', type=int, default=128, help="Merge : input patch's size.")
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
        if args.data.endswith('/'):
            pass
        else:
            args.data = args.data + "/"
        if args.iden.endswith('/'):
            pass
        else:
            args.iden = args.iden + "/"

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
def fragmentMRC(inputFile, size, inputIdent, outputIdent):
    with mrcfile.open(inputFile) as mrc:
        print("Input file shape is : ", mrc.data.shape)
        assert mrc.data.shape[1] % size == 0
        assert mrc.data.shape[2] % size == 0
        wfactor, hfactor = (mrc.data.shape[1] // size), (mrc.data.shape[2] // size)
        print("Output patches will be generated : {} x {} patches".format(wfactor, hfactor))
        
        identifier = inputIdent[:len(inputFile) - 4]

        for w in range(wfactor):
            for h in range(hfactor):
                patchName = "{}_w{}_h{}.mrc".format(identifier, w , h)
                with mrcfile.new(outputIdent + patchName) as mrcPat:
                    mrcPat.set_data(mrc.data[:, w * size : (w+1) * size , h * size : (h+1) * size ])
                    # print("  patch Name : {}".format(patchName))
                    # print("    patches : ({},{}) ~ ({},{}) ".format(w * size, h * size, w * size + size - 1, h * size + size -1))

def mergeMRC(args):
    # use args.identifier, w, f, output identifier.
    patchNamePattern = args.data # input file name pattern.
    if args.iden == "":
        outputName        = args.patchNamePattern
    else:
        outputName        = args.iden
    print("Input file pattern : ", patchNamePattern)
    print("Output file name   : ", outputName)

    mergedContainer = np.zeros((args.fragsize * args.width, args.fragsize * args.height),dtype=np.float32)
    for w in range(args.width):
        for h in range(args.height):
            patchName = patchNamePattern.replace("{w}", str(w), 1).replace("{h}", str(h), 1)
            with mrcfile.open(patchName) as mrcPatch:
                #print(mrcPatch.data.shape)
                mergedContainer[w * args.fragsize : (w+1) * args.fragsize , h * args.fragsize : (h+1) * args.fragsize] = mrcPatch.data

    with mrcfile.new(outputName) as mergedMRC:
        mergedMRC.set_data(mergedContainer)
    
if __name__ == '__main__':
    print("------------ fragmentMRC started... ------------")
    withError, args = parse_args()
    if withError:
        pass
    else:
        #print("Output size will be : ", args.size)
        #print("Input file is : ", args.data)
        if args.isMerge:
            print("-- Merge mode.")
            assert False
            mergeMRC(args)
        else:
            print("-- Fragment mode.")
            input_list = os.listdir(args.data)
            for i in input_list:
                fileName = args.data + i
                if fileName.endswith('.mrc'):
                    print(fileName)
                    fragmentMRC(fileName, args.size, i, args.iden)
    print("------------ fragmentMRC ended... ------------")


# python fragmentWrapper.py -s 256 -d /cdata/WrapperTEM/Micrographs/noisy4/ -id /cdata/db1/fragmentedTEM/
# python /cdata/workspace/fragmentWrapper.py -s 512 -d /cdata/WrapperTEM/Micrographs/noisy4/ -id /cdata/db1/size512/fragments/noisy4
# python /cdata/workspace/fragmentWrapper.py -s 512 -d /cdata/WrapperTEM/Micrographs/clean4/ -id /cdata/db1/size512/fragments/clean4