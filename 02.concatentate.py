''' CONCATENATE SAMPLES
This script takes a folder containing paired-end trimmed fastq files, finds multiple forward and multiple reverse pairs for each
sample and concatenates them into one forward and one reverse file per sample.
It is important that the order of files concatenated is preserved between forward and reserve files for each sample.
'''
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen, list2cmdline
from itertools import product
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder of trimmed fastq files to concatenate")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def exec_in_row(cmds):
    ''' Exec commands one after the other until finished.'''
    if not cmds:
        return  # empty list
    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)
    for task in cmds:
        print task
        p = Popen(task, shell=True)
        p.wait()
    if done(p):
            if success(p):
                print "done!"
            else:
                fail()

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("paired.fastq")]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infiles = list_folder(args.infolder)
    print "\nNumber of infiles (fastq) =", len(infiles)

    #Creates dictionary of forward and reverse pairs
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        if name.split("_")[-4] == "R1": #depends on naming of files; this should be something uniqe to the forward or reverse files (eg R1/R2)
            uniqid = name.split("R1")[0] #uniqid should include a uniq identifier for the lane
        elif name.split("_")[-4] == "R2": #depends on naming of files
            uniqid = name.split("R2")[0] #uniqid should include a uniq identifier for the lane
        filedictionary[uniqid].append(infile)
    print "Number of pairs of infiles (fastq) =", len(filedictionary)

    #Check each sample has a forward and a reverse
    for sample in filedictionary:
        check = []
        for s in filedictionary[sample]:
            uniqid = s.split("_")[-4] #depends on naming of files; the uniqid here should be something uniqe to the forward or reverse files (eg R1/R2)
            check.append(uniqid)
        if len(check) == 2 and str("R1") in check and str("R2") in check:
                pass
        else:
            print "ERROR - sample", sample, "missing forward or reverse file"
            break

    #Creates dictionary of forward and reverse pairs for each sample - ensures file order is preserved!
    samplefiledictionary = defaultdict(list)
    for file in filedictionary:
        sample = file.split("_")[-1]
        for f in filedictionary[file]:
            splitter = "_"+sample+"_"
            uniqid = sample+"_"+f.split(splitter)[1].split("_")[0]
            samplefiledictionary[uniqid].append(f)
    print "Number of samples =", len(samplefiledictionary)/2

    #Check that order is preserved
    checksamplefiledictionary = defaultdict(list)
    for sample in samplefiledictionary:
        sample_id = sample.split("_")[0]
        checksamplefiledictionary[sample_id].append(samplefiledictionary[sample])
    for sample_id in checksamplefiledictionary:
        splitter = "_"+sample_id+"_"
        checkforward = []
        checkreverse = []
        for c in checksamplefiledictionary[sample_id]:
            for file in c:
                lane = file.split(splitter)[0].split("_")[-1]
                direction = file.split(splitter)[1].split("_")[0]
                if direction == 1:
                    checkforward.append(lane)
                elif direction == 2:
                    checkreverse.append(lane)
        if checkforward == checkreverse:
            pass
        else:
            print "ERROR - lane order not preserved"
            print samplefiledictionary
            break

    #Make concatenate commands
    makeconcatenaterun = []
    for sample in samplefiledictionary:
        files = samplefiledictionary[sample]
        direction = sample.split("_")[-1]
        if direction == "1":
            outfile = os.path.dirname(files[0])+"/"+os.path.basename(files[0]).split("_")[0]+"_"+sample+"_forward_paired.fastq"
        elif direction == "2":
            outfile = os.path.dirname(files[0])+"/"+os.path.basename(files[0]).split("_")[0]+"_"+sample+"_reverse_paired.fastq"

        infilenames = files[0]
        for file in samplefiledictionary[sample][1:]:
            infilenames = infilenames+" "+file
        
        concatenatecommand = ["cat "+str(infilenames)+" > "+outfile]
        makeconcatenaterun.append(concatenatecommand)

    #Run concatenate
    print "Number of commands to run =", len(makeconcatenaterun)
    exec_in_row(makeconcatenaterun)

if __name__ == "__main__":
    main()