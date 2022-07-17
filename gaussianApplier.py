import argparse,os,time,mrcfile
import numpy as np
from os.path import exists
import sys

# Parsing arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Gaussian noise applier for MRC file")
    parser.add_argument('-i', '--input', required=True, help="Path to dataset")
    parser.add_argument('-o', '--output', required=True, help="Path to gaussian noised output")
    parser.add_argument('-a', '--augment', type=int, default=1, help="The number of the output from each input mrc")
    parser.add_argument('--std', type=float, default=1.0, help="Parameter for the std multiply")
    return check_args(parser.parse_args())

# Checking arguments
def check_args(args):
    # --epoch
    withError = False
    if not exists(args.input):
        withError = True
        print('Input Directory should exist : {}'.format(args.data))
    elif not args.input.endswith("/"):
        print("WARNING : input directory should endswith /")
        args.input = args.input + "/"
    
    if not exists(args.output):
        withError = True
        print('Output Directory should exist : {}'.format(args.output))
    elif not args.output.endswith("/"):
        print("WARNING : output directory should endswith /")
        args.output = args.output + "/"
    
    return withError, args

#################### MAIN ####################
def applier_Adder(cleanNumpy, noiseNumpy):
    return cleanNumpy + noiseNumpy
    


if __name__ == '__main__':
    withError, args = parse_args()
    if withError:
        pass
    else:
        script_start = time.time()
        input_list=os.listdir(args.input)

        avgTime = None
        for inputNum, inputName in enumerate(input_list, start=1):
            file_start = time.time()
            # print("Job#{} on {}".format(inputNum, inputName))
            with mrcfile.open(args.input + inputName) as mrcPatch:
                mrcNumpy = mrcPatch.data        # (1, Size, Size)
            size     = mrcNumpy.shape[1]    # Size
            mrcAvg   = np.average(mrcNumpy)
            mrcStd   = np.std(mrcNumpy)
            # print(" - input mrc file's avg {} and std {}".format(mrcAvg, mrcStd))

            if args.augment == 1:
                outputName = inputName
                gaussianNoise   = np.random.normal(loc= mrcAvg, scale= args.std * mrcStd, size=(size, size))
                applied_patch   = applier_Adder(mrcNumpy, gaussianNoise).astype(np.float32)
                with mrcfile.new(args.output + outputName) as mrcApplied:
                    mrcApplied.set_data(applied_patch)
                    # print(" - saved to {}".format(args.output + outputName))
            else:
                for augNum in range(args.augment):
                    outputName = inputName[:-4] + "_" + str(augNum) + ".mrc"
                    gaussianNoise   = np.random.normal(loc= mrcAvg, scale= args.std * mrcStd, size=(size, size))
                    applied_patch   = applier_Adder(mrcNumpy, gaussianNoise).astype(np.float32)
                    with mrcfile.new(args.output + outputName) as mrcApplied:
                        mrcApplied.set_data(applied_patch)
                        # print(" - saved to {}".format(args.output + outputName))
            file_end = time.time()
            if avgTime == None:
                avgTime = file_end - file_start
            else:
                avgTime = (file_end - script_start) / (inputNum)
            print('# gaussApplier : [{}/{}] {:.2%} | ETA : {}'.format(inputNum, len(input_list), inputNum/len(input_list), avgTime * (len(input_list) - inputNum)), file=sys.stderr, end='\r')
        print("Total elapsed time for gaussApplier : {}".format(time.time() - script_start))

# python /cdata/workspace/gaussianApplier.py -i /cdata/WrapperTEM/Micrographs/clean4/ -o /cdata/db1/tempGauss/ -a 2 --std 10
# python /cdata/workspace/gaussianApplier.py -i /cdata/tomsimDIR/0714_EMAN2_2wrj_cleanMRC -o /cdata/tomsimDIR/0714_EMAN2_2wrj_gaussMRC -a 2 --std 10