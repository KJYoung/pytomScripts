import sys,os

srcDir    = sys.argv[1]  # micrographs
dstDir    = sys.argv[2]  # raw micrographs
augNum    = int(sys.argv[3])  # number of workers.

input_list= os.listdir(srcDir)

for i in input_list:
    srcName = i
    for j in range(augNum):
        dstName = srcName[:-4] + "_{}".format(j) + ".mrc"
        command = "cp {}{} {}{}".format(srcDir, srcName, dstDir, dstName)
        os.system(command)
    print("Jobs for {}".format(srcName))
    # print('# bulkRenamer : [{}/{}] {:.2%}'.format(i, len(input_list), i/len(input_list)), file=sys.stderr, end='\r')
print("All jobs done.")

# python /cdata/workspace/bulkRenamer.py /cdata/db1/size512/fragments/clean4 /cdata/db1/size512/fragments/cleanPair 2