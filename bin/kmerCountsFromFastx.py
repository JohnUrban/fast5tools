#!/usr/bin/env python2.7
import h5py, os, sys, csv
import argparse
from glob import glob
import string
from cStringIO import StringIO
from collections import defaultdict
from Bio import SeqIO
from math import log10, log
import numpy as np

## R stuff
try:
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
except:
    sys.stderr.write('Differential representation')

# Plotting
import matplotlib.pyplot as plt

# Fast5Tools
from fast5tools.f5class import *
from fast5tools.f5ops import *
from fast5tools.helperops import *
from fast5tools.fileListClass import *
from fast5tools.plotops import *

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s,
- get kmer count tables
- see differentially-represented kmers (given reference set)




Note:
IF you have the right rpy2 setup,
you can calculate differentially represented kmers using EdgeR from bioconductor.
- http://bioconductor.org/packages/2.5/bioc/html/edgeR.html
The EdgeR analysis is derived from Brad Chapman's example at:
https://github.com/chapmanb/bcbb/blob/master/stats/count_diffexp.py

You can always do the median-normalization analysis.
Instead of TMM, the median is subtracted from each count,
    followed by normalizing it to the MAD.



John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-k', '--k', type=int, default=6,
                    help = '''Kmer-size. Default: 6.''')

parser.add_argument('--revcomp', default=False, action='store_true',
                    help = '''Also count kmers from reverse complement of each sequence.''')



parser_not_fast5 = parser.add_mutually_exclusive_group()
parser_not_fast5.add_argument('--fasta', '-fa', action='store_true', default=False,
                           help='''Looking at a FASTA file or list of FASTA files, not FAST5s''')
parser_not_fast5.add_argument('--fastq', '-fq', action='store_true', default=False,
                           help='''Looking at a FASTQ file or list of FASTQ files, not FAST5s''')

#nargs='+',
parser.add_argument('--reference', 
                    type=str, default=None, help='''All files after this flag and before the next, are interpreted as Reference fastA/Q/5 files.
NOTE: Unlike the default datasets that can contain as many files/dirs/fofns as you'd like,
this can only be pointed at 1 object UNLESS you put everything between quotation marks, which allows you to
specify as many reference files/dirs/FOFNs as you'd like.
E.g. --reference "fast5dir1/ fast5dir2/" ''')

parser_ref_not_fast5 = parser.add_mutually_exclusive_group()
parser_ref_not_fast5.add_argument('--reffasta', '-rfa', action='store_true', default=False,
                           help='''The reference dataset is a FASTA file or list of FASTA files, not FAST5s''')
parser_ref_not_fast5.add_argument('--reffastq', '-rfq', action='store_true', default=False,
                           help='''The reference dataset is a FASTQ file or list of FASTQ files, not FAST5s''')



parser.add_argument('-r', '--readtype', default='template',
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d', 'molecule', 'all', 'MoleQual'.
Default: template.
Molecule returns single fasta for each fast5 by following rules:
    if 2d present, return 2d.
    elif complement present with no 2d, return longer of template or complement.
    elif only template present, return template.
'MoleQual' is similar to molecule.
    It differs only in choosing between template and complement when a 2D is not present.
    Instead of choosing the longer one, it chooses the one with a higher quality mean quality score.''')




parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')

parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')

parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')

parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')





parser.add_argument('-n', '--nfiles', type=int, default=1000000000000,
                    help = '''This defaults to 1000000000000 in order to use all files (will you ever need to look at more than that?).
However, you can downsample with this option by adjusting this number to get info from the first N files.
Use --random to select N at random from the list.
Aim this script at a specific file for that file's contents.''')

parser.add_argument('-R', '--random', action='store_true', default=False,
                    help = '''Randomize what files are looked at.''')

parser.add_argument('-S', '--randomseed', type=int, default=False,
                    help = '''Randomize what files are looked at, but use given seed for reproducibility.''')

parser.add_argument('--filesused', type=str, default='qual_v_pos', help='''
''')

parser.add_argument('-o', '--outdir', type=str, default="./",
                    help = '''....''')

parser.add_argument('--filename', type=str, default=None, help='''
For output. Default None (stdout).
If a filename is given, filesused will be reported in a similarly named file ending with .filesused.fofn''')





parser.add_argument('--notarlite', action='store_true', default=False, help=''' The default methof (called tarlite) extracts 1 file from a given tarchive at a time, processes, and deletes it.
This options says to turn tarlite off resulting in extracting entire tarchive before proceeding (and finally deleting).
It is possible that --notarlite is faster, but at the expense of exceeding file number limits or disk storage quotas.
Nonetheless, the difference in speed is a lot smaller than the difference in space needed.
For example, not using tarlite will require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
A corollary is that tarlite just needs to be allowed to form 1 (or a few) files compared to what could be thousands to millions.
''')

parser.add_argument('--tarlite', action='store_true', default=False, help='''This legacy option is outdated.
However, it is kept here to avoid breaking pipelines that make use of it.
The tarlite approach is now default. Specifying this will not change that default behavior.
It will just prevent pipelines from breaking.
However, not specifying this will still also result in the tarlite approach.
Use --notarlite to turn it off.''')

args = parser.parse_args()


## run_kmer_counting() used to be defined here. I moved it to helperops in the kmer section.
## That required adding imports including fast5tools imports.
## If that interferes w/ other scripts down the road, then I will change back or make a specific kmerops file.

if __name__ == "__main__":
    
    # Process Args
    args.outdir = process_outdir(args.outdir)
    outfile = args.outdir + args.filename if (args.filename is not None) else None



    ## Execute
    kmerdict, filesused = run_kmer_counting(initial_list=args.fast5, \
                                 k=args.k, \
                                 readtype=args.readtype, \
                                 revcomp=args.revcomp, \
                                 nfiles=args.nfiles, \
                                 random=args.random, \
                                 randomseed=args.randomseed, \
                                 notarlite=args.notarlite, \
                                 fasta=args.fasta, \
                                 fastq=args.fastq, \
                                 minlen=args.minlen, \
                                 maxlen=args.maxlen, \
                                 minq=args.minq, \
                                 maxq=args.maxq)


    ## Write
    writekmer(kmerdict, outfile)
    
    ## Files used
    process_filesused(trigger=args.filename, filesused=filesused, outdir=args.outdir)

    
    ## Reference?
    if args.reference is not None:
        ## Out
        refoutfile = args.outdir + 'reference_' + args.filename if (args.filename is not None) else None
        ## Convert into list:
        args.reference = args.reference.strip().split()
        refdict, refsused = run_kmer_counting(initial_list=args.reference, \
                                 k=args.k, \
                                 readtype=args.readtype, \
                                 revcomp=args.revcomp, \
                                 nfiles=args.nfiles, \
                                 random=args.random, \
                                 randomseed=args.randomseed, \
                                 notarlite=args.notarlite, \
                                 fasta=args.fasta, \
                                 fastq=args.fastq, \
                                 minlen=args.minlen, \
                                 maxlen=args.maxlen, \
                                 minq=args.minq, \
                                 maxq=args.maxq)
        ## Write
        writekmer(refdict, refoutfile)
        
        ## Files used
        process_filesused(trigger=args.filename, filesused=refsused, outdir=args.outdir+'reference_')




## PLOTTING
#singleTableKmerPlot(kmerdict)
#singleTableKmerHist(kmerdict)
#if args.reference is not None:
    #singleTableKmerPlot(refdict)
    #twoTableKmerPlot(kmerdict, refdict)
    #twoTableKmer_MA_Plot(kmerdict, refdict)
    #twoTableKmer_MA_Plot(kmerdict, refdict, scale_to_ref=True, logplot=True)



        
    
