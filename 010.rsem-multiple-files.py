'''RSEM
This script will take a folder containing PE fastq reads and run RSEM calculate-expression
against a prep-referenced for RSEM
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
                    help="Infolder of fastq.gz files (trimmed) for RSEM")
parser.add_argument("pathtoRSEM", type=str,
                    help="This should be the pathway to RSEM: rsem-calculate-expression")
parser.add_argument("rsemreference", type=str,
                    help="This is the file that your rsem reference is based from, using the rsem-prep-reference command")

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
    '''Exec commands one after the other until finished'''
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith("fastq.gz")]
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infiles = list_folder(args.infolder)
    print "\nNumber of infiles (fastq.gz) =", len(infiles)

    #Creates dictionary of forward and reverse pair for each sample
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        uniqid = name.split("_")[-6] #depends on naming of files; the uniqid should be something that is unique to each sample
        filedictionary[uniqid].append(infile)
    print "Number of samples =", len(filedictionary) #check that this number is half the total number of infiles printed above as each sample has two files, a forward and a reverse

    #Check that each sample has a forward and a reverse
    for sample in filedictionary:
        check = []
        for s in filedictionary[sample]:
            uniqid = s.split("_")[-4] #depends on naming of files; the uniqid here should be something uniqe to the forward or reverse files (eg R1/R2; forward/reverse)
            check.append(uniqid)
        if len(check) == 2 and str("R1") in check and str("R2") in check: #depends on naming of files
                pass
        else:
            print "ERROR - sample", sample, "missing forward or reverse file"
            break

    #Prepare RSEM commands
    print "Path to RSEM =", args.pathtoRSEM
    makeRSEMrun = []
    for sample in filedictionary:
        for s in filedictionary[sample]:
            if s.split("_")[-4] == str("R1"): #depends on naming of files
                forward = s
            elif s.split("_")[-4] == str("R2"): #depends on naming of files
                reverse = s
        sample_forward = forward
        sample_reverse = reverse
        sample_individ = forward.split("_")[-6]
        rsemcommand = [args.pathtoRSEM+" --paired-end "+sample_forward+" "+sample_reverse+" "+
                        " --bowtie2 "+ args.rsemreference+" "+sample_individ]
        makeRSEMrun.append(rsemcommand)

    #Run RSEM
    print "Number of commands to run =", len(makeRSEMrun), "\n"
    exec_in_row(makeRSEMrun)
    print "\n"

if __name__ == "__main__":
    main()

