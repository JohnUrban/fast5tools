#!/usr/bin/env python2.7
import h5py, os, sys, csv
import argparse
# General
from glob import glob
import string
from cStringIO import StringIO
from collections import defaultdict
from Bio import SeqIO
from math import log10, log
import numpy as np

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
- read length histogram
    - bins have counts
- cumulative read length histogram
    - bins have cumulative counts
- read length sum histogram
    - bins have sums of readlengths in that bin
- cumulative read length histogram
    - bins have cumulative sums of readlengths in each bin

- RL v Q

TODO:
If reference dataset given, the above plots for both, and:

    - Subtraction plot:
        - Proprtional_test_bin - proportional_ref_bin

    - Log Fold Change plot:
        - log(Proprtional_test_bin / proportional_ref_bin)


John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')



parser_not_fast5 = parser.add_mutually_exclusive_group()
parser_not_fast5.add_argument('--fasta', '-fa', action='store_true', default=False,
                           help='''Looking at a FASTA file or list of FASTA files, not FAST5s''')
parser_not_fast5.add_argument('--fastq', '-fq', action='store_true', default=False,
                           help='''Looking at a FASTQ file or list of FASTQ files, not FAST5s''')


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

parser.add_argument('--bin-start', dest='bin_start', type=str, default='0', help='''Start binning at given read length. Default: 0 bp.
Can give number as bp, which requires no units - e.g. 500 = 500 bp.
Can also give as kb or Mb, but it requires the unit -- e.g. 500kb or 0.5Mb.
No space between the number and unit.''')
parser.add_argument('--bin-end', dest='bin_end', type=str, default='1Mb', help='''End binning at given read length.
Default: 1Mb.
Note that sometimes providing smaller end values (e.g. 100kb) makes for better plots.
Note that this plot purposely is independent of max read length as read length distributions trpically have long sparse tails, which make these plots front-heavy and end-sparse.
Other avenues (such as reporting the max length) might be better.
Nonetheless, feel free to give huge end values.''')
parser.add_argument('--bin-width', dest='bin_width', type=str, default='1kb', help='''Specify the width of bins in bp. Default: 1000 bp.''')



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

parser.add_argument('--filename', type=str, default='kmer_counts.txt', help='''
For output. Default: kmer_counts.txt. (Formerly defaulted to None).
If a filename is given, filesused will be reported in a similarly named file ending with .filesused.fofn
When --reference used, files will have similar name with reference_ prepended.''')

parser.add_argument('--plotfilesuffix', type=str, default=None, help='''
Suffix and extension for output plots. Default None (PDFs output in outdir using hard-coded prefixes).
Plots will be in specified outdir.
The minimum information to give is the extension (no dot needed) -- e.g. png, jpg, pdf.
Example1: myexperiment.png ((plots will be named plotprefix_myexperiment.png))
Example2: .jpg ((plots will be named plotprefix_.jpg))
Example3: jpg ((plots will be named plotprefix.jpg))
Example4: when None (default), plots will be named plotprefix.pdf''')




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




## for fa in fasta or fq in fastq or f5 in fast5 -- get all lengths, and meanQ
## make DP objs for RL and Q



if __name__ == "__main__":
    
    # Process Args
    args.outdir = process_outdir(args.outdir)
    outfile = args.outdir + args.filename if (args.filename is not None) else None
    start = interpret_base_length(args.bin_start)
    end = interpret_base_length(args.bin_end)
    bw = interpret_base_length(args.bin_width)


    ## Execute
    test_lengths, filesused = run_collect_lengths(initial_list=args.fast5, \
                                 readtype=args.readtype, \
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


    ## Get Partition
    test = DataPartition(test_lengths, start, end, bw)

    
    ## Write
    ## TODO: writeout partition
    
    ## Files used
    process_filesused(trigger=args.filename, filesused=filesused, outdir=args.outdir)

    
    ## Reference?
    do_comparison = False
    if args.reference is not None:
        do_comparison = True        
        ## Out
        refoutfile = args.outdir + 'reference_' + args.filename if (args.filename is not None) else None

        ## Convert into list:
        args.reference = args.reference.strip().split()

        ## Get object
        ref_lengths, refsused = run_collect_lengths(initial_list=args.reference, \
                                 readtype=args.readtype, \
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
        ## Get Partition
        ref = DataPartition(ref_lengths, start, end, bw)

        ## Write
        ## TODO: writeout partition
        ## Files used
        trigger = 'reference_'+args.filename if args.filename is not None else None
        process_filesused(trigger=trigger, filesused=refsused, outdir=args.outdir)



'''
Given path(s) to fast5 file(s) and/or directories of fast5s,
- read length histogram
    - bins have counts
- cumulative read length histogram
    - bins have cumulative counts
- read length sum histogram
    - bins have sums of readlengths in that bin
- cumulative read length histogram
    - bins have cumulative sums of readlengths in each bin

- RL v Q

TODO:
If reference dataset given, the above plots for both, and:

    - Subtraction plot:
        - Proprtional_test_bin - proportional_ref_bin

    - Log Fold Change plot:
        - log(Proprtional_test_bin / proportional_ref_bin)

'''
        
## COMPARITIVE ANALYSES
if do_comparison:

    ## DETERMINE PLOT SUFFIX 
    if args.plotfilesuffix is None:
        sfx = '.pdf'
    else:
        if '.' in args.plotfilesuffix:
            ## assumes format words.ext, makes sfx = _words.ext
            sfx = '_' + args.plotfilesuffix
        else:
            ## assumes image type specified (e.g. pdf, jpg, png)
            sfx = '.' + args.plotfilesuffix
    make_name = make_name_function(pfx=args.outdir, sfx=sfx)

    ## PLOTTING

    
    ## TODO: Add subtraction foldchange stuff
    
    ## Counts
    general_barplot(x=test.get_breaks()[:-1], height=test.get_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_read_length_bin_counts'))
    ## CumulativeCounts
    general_barplot(x=test.get_breaks()[:-1], height=test.get_cum_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_cumulative_read_length_bin_counts'))
    ## Data
    general_barplot(x=test.get_breaks()[:-1], height=test.get_heights(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_read_length_data_per_bin'))
    ## Cumlative Data
    general_barplot(x=test.get_breaks()[:-1], height=test.get_cum_heights(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_cumulative_read_length_data_per_bin'))

    ## Proportional Counts
    general_barplot(x=test.get_breaks()[:-1], height=test.get_proportion_total_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_proportional_read_length_bin_counts'))
    ## Proportional CumulativeCounts
    general_barplot(x=test.get_breaks()[:-1], height=test.get_proportion_cum_total_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_proportional_cumulative_read_length_bin_counts'))
    ## Proportional Data
    general_barplot(x=test.get_breaks()[:-1], height=test.get_proportion_total_height(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_proportional_read_length_data_per_bin'))
    ## Proportional Cumlative Data
    general_barplot(x=test.get_breaks()[:-1], height=test.get_proportion_cum_total_height(), width=bw, edgecolor='k', align='edge', saveas=make_name('test_proportional_cumulative_read_length_data_per_bin'))
    
    if args.reference is not None:
        ## Counts
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_read_length_bin_counts'))
        ## CumulativeCounts
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_cum_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_cumulative_read_length_bin_counts'))
        ## Data
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_heights(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_read_length_data_per_bin'))
        ## Counts
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_cum_heights(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_cumulative_read_length_data_per_bin'))

        ## Proportional Counts
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_proportion_total_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_proportional_read_length_bin_counts'))
        ## Proportional CumulativeCounts
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_proportion_cum_total_counts(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_proportional_cumulative_read_length_bin_counts'))
        ## Proportional Data
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_proportion_total_height(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_proportional_read_length_data_per_bin'))
        ## Proportional Cumlative Data
        general_barplot(x=ref.get_breaks()[:-1], height=ref.get_proportion_cum_total_height(), width=bw, edgecolor='k', align='edge', saveas=make_name('ref_proportional_cumulative_read_length_data_per_bin'))
   

        ##SUBTRACTION
        general_barplot(x=ref.get_breaks()[:-1], height=test.get_difference_in_proprtional_counts(ref), width=bw, edgecolor='k', align='edge', saveas=make_name('difference_in_read_length_bin_counts'))
        general_barplot(x=ref.get_breaks()[:-1], height=test.get_difference_in_proprtional_heights(ref), width=bw, edgecolor='k', align='edge', saveas=make_name('difference_in_read_length_data_per_bin'))

        ## LOG FOLD CHANGE
        general_barplot(x=ref.get_breaks()[:-1], height=test.get_log_fold_change_of_proportional_counts(ref), width=bw, edgecolor='k', align='edge', saveas=make_name('logfoldchange_of_read_length_bin_counts'))
        general_barplot(x=ref.get_breaks()[:-1], height=test.get_log_fold_change_of_proportional_heights(ref), width=bw, edgecolor='k', align='edge', saveas=make_name('logfoldchange_of_read_length_data_per_bin'))
  
        
    
