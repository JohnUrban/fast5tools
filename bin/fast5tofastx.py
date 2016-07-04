#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
from fast5tools.f5class import *
from fast5tools.f5ops import *
import argparse
from glob import glob


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s, return fasta, fastq, qual, or intqual for all fast5s found.

For files that are corrupt or empty, for now it silently skips them.
As an alternative, fast5stats will tell you all files skipped (in stderr or to specified file).

John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="mol",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d', 'molecule', 'all'.
Default: molecule.
There is no need to write full word for options - can do: t, c, 2, m, a.
Molecule returns single fasta for each fast5 by following rules:
if 2d present, return 2d.
elif complement present with no 2d, return longer of template or complement.
elif only template present, return template.''')

parser.add_argument('-o', '--outtype', type=str, default="fasta",
                    help = '''Choices: fasta, fastq, qual, intqual, details.
Default: fasta.
If details, sequence not reported, but name, seqlen, and meanq are.''')

parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')

parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')

parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')

parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')


args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
assert args.outtype in ("fasta", "fastq", "qual", "intqual", "details")
assert args.readtype[0] in "tc2ma"
if args.readtype[0] == "t":
    args.readtype = "template"
elif args.readtype[0] == "c":
    args.readtype = "complement"
elif args.readtype[0] == "2":
    args.readtype = "2d"
elif args.readtype[0] == "m":
    args.readtype = "molecule"
elif args.readtype[0] == "a":
    args.readtype = "all"



#################################################
### uses output functions from f5ops.py
### fasta(), fastq(), qual(), intqual()
### 
#################################################



#################################################
## fast5tofastx specific "output" functions
#################################################

def details(f5, readtype):
    readstats = []
    readstats.append( f5._get_pore_info_name(readtype) )
    readstats.append( f5.get_seq_len(readtype) )
    readstats.append( f5.get_mean_qscore(readtype) )
    readstats.append( f5.get_num_events(readtype) )
    try:
        readstats.append( f5.get_num_called_events(readtype) )
    except:
        readstats.append("-")
    try:
        readstats.append( f5.get_num_skips(readtype) )
    except:
        readstats.append("-")
    try:
        readstats.append( f5.get_num_stays(readtype) )
    except:
        readstats.append("-")        
    return ("\t").join([str(e) for e in readstats])



    
#################################################
####### argument processing functions ###########
#################################################

def get_fast5tofastx_fxns(args):
    ### get outtype fxn ###
    if args.outtype == "fasta":
        output = fasta
    elif args.outtype == "fastq":
        output = fastq
    elif args.outtype == "qual":
        output = qual
    elif args.outtype == "intqual":
        output = intqual
    elif args.outtype == "details":
        output = details
    ### get readtype fxn ###
    if args.readtype == "template":
        getread = get_template_read
    elif args.readtype == "complement":
        getread = get_complement_read
    elif args.readtype == "2d":
        getread = get_2d_read
    elif args.readtype == "molecule":
        getread = get_molecule_read
    elif args.readtype == "all":
        getread = get_all_reads
    return output, getread



            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    output, getread = get_fast5tofastx_fxns(args)

    for f5 in Fast5List(args.fast5):
        if f5.is_not_corrupt() and f5.is_nonempty:
            read = getread(f5, args.minlen, args.maxlen, args.minq, args.maxq, output)
            if read:
                print read


