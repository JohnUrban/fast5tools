#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import f5class
from f5ops import *
import argparse
from glob import glob


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s, return fasta, fastq, qual, or intqual for all fast5s found.

John Urban (2015)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="input",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', 'input'.
Default: input.
''')

parser.add_argument('-o', '--outdir', type=str, default="",
                    help = '''If single fast5 specified, it will be reported to stdout.
If multiple fast5s are specified, they will be saved to files in the working dir by default.
This flag allows you to specify a different output directory.
Filenames will be the the name of the fast5 file with .events.txt appended.''')



args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
assert args.readtype[0] in "tci"
if args.readtype[0] == "t":
    args.readtype = "template"
elif args.readtype[0] == "c":
    args.readtype = "complement"
elif args.readtype[0] == "i":
    args.readtype = "input"
if args.outdir:
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    if args.outdir[-1] != "/":
        args.outdir += "/"



            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    f5list = f5class.Fast5List(args.fast5)
    if len(f5list) == 1:
        for f5 in f5list:
            if f5.has_read(args.readtype):
                print f5.get_events_string(args.readtype)
    elif len(f5list) > 1:
        for f5 in f5list:
            if f5.has_read(args.readtype):
                out = open(args.outdir + f5.get_base_info_name() + "." + args.readtype + ".events.txt", 'w')
                out.write(f5.get_events_string(args.readtype))
                out.close()

        


