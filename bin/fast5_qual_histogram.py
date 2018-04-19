#!/usr/bin/env python2.7


import h5py, os, sys
from fast5tools.f5class import *
from fast5tools.f5ops import *
from fast5tools.helperops import process_outdir
from fast5tools.plotops import qualhist
import argparse
from glob import glob
from collections import defaultdict
import numpy as np
import matplotlib
## may need following line for remote jobs (e.g. submitting batch scripts)
## matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from random import shuffle, seed

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s,
give histogram of base qualities encountered.

John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-Q', '--mean-quality-scores', dest='mean_quality_scores', default=False, action='store_true',
                    help='''Instead of extracting all base-qualities from each read, only extract its mean quality score (Q).''')

parser.add_argument('-B', '--bin-range', type=str, dest='bin_range', default=None,
                    help='''To be consistent between plots, the same bins should be populated. Give comma-sep trio of integers describing minRange, maxRange, step_size.
If all base qualities, default: 0,32,1
If obly mean qualities, default: 0,16,1''')

parser.add_argument('-D', '--density', action='store_true', default=False,
                    help='''Plot density instead of frequency.''')

parser.add_argument('-C', '--cumulative', action='store_true', default=False,
                    help='''Plot cumulative distribution instead of regular.''')

parser.add_argument('-N', '--no-plot', dest='no_plot', action='store_true', default=False,
                    help='''Instead of plotting, just print to stdout as 2-columns of text: qual_score, frequency/density/etc.
Recommend to use with --filename and making the --filename end with .txt, not an image extension.''')

##parser.add_argument('--bin-width', dest='bin_width', default=1000, type=int, help='''The width of bins (default: 1000 bp).''')

##parser.add_argument('-z', '--zscores', default=False, action='store_true', help='''For each read, its quals are converted to z-scores such that each position's score is relative to the distribution of qual scores on that read.
##This prevents issues where there appears to be a drop off in quality as reads get longer that arises simply b/c there are some long reads with low quality across the entire read.
##In that sceario, using Z-scores shows that the quality does not drop as the read gets longer.
##Can also filter reads for mean quality score before plotting to prevent this effect.''')
##
##parser.add_argument('-Z', '--robust-zscores', dest='robust_zscores', default=False, action='store_true', help='''Same as --zscores, only Median and MAD used here.''')

parser.add_argument('-n', '--nfiles', type=int, default=100,
                    help = '''Often when this command is run, looking at 100 files will suffice.
Therefore, this defaults to 1, and looks at the first file in the FileList. Aim this script at a specific file for that file's contents.
Adjust this number to get info from the first N files.
Use --random to select N at random from the list.
Set this to super high number for all fast5files  -- e.g. 1000000000000''')

parser.add_argument('-R', '--random', action='store_true', default=False,
                    help = '''Randomize what files are looked at.''')

parser.add_argument('-S', '--randomseed', type=int, default=False,
                    help = '''Randomize what files are looked at, but use given seed for reproducibility.''')

parser.add_argument('-o', '--outdir', type=str, default="./",
                    help = '''....''')

parser.add_argument('--filename', type=str, default=None, help='''This will use the extension given to determine outputfile type.
Examples: qualhist.jpg, qualhist.png, qualhist.png.
Default: None
Default will have a window pop up with the plot instead. The file can be saved from there too.
If you choose this option, the filesused will be reported to stderr.
If a filename is given, filesused will be reported in a similarly named file ending with .filesused.fofn''')

parser.add_argument('--filesused', type=str, default='qual_v_pos', help='''
''')

parser.add_argument('-r', '--readtype', default="molecule",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d', 'molecule', 'all', 'MoleQual'.
Default: molecule.
There is no need to write full word for options - can do: t, c, 2, m, a, M.
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

#################################################
## This can be done from f5stats
#################################################
## qual v readlen scatter
##if read_lengths:
##    logger.info("Constructing quality vs. readlen scatter plot...")
##    plt.scatter(read_lengths, mean_scores)
##    plt.xlabel("Read Length (bp)")
##    plt.ylabel("Mean Quality score")
##    plot_out(args, "meanqual-Vs-readlen")
##    plt.close()
##else: logger.warn("No reads that meet criteria: cannot construct quality vs. readlen scatter plot...")

##    if args.qualVsLen:
##        x = list()
##        y = list()
##                elif args.qualVsLen:
##                    qscores = []
##                    for q in fq.qual:
##                        qscores.append(ord(q)-33)
##                    qmean = np.mean(qscores)
##                    qLen = len(qscores)
##                    x.append(qLen)
##                    y.append(qmean)
##                    
##                else:
##                    for q in fq.qual:
##                        ctr += 1
##                        qualpos[1+int(ctr//bin_width)].append(ord(q)-33)
##

#################################################
##
#################################################


if __name__ == "__main__":
    
    # Process Args
    if args.filename is not None:
        filesused_h_ = '.'.join(args.filename.split('.')[:-1]) + '.filesused.fofn'
        filesused_h = args.outdir + filesused_h_ if args.outdir.endswith('/') else args.outdir + '/' + filesused_h_
        outfile = args.outdir + args.filename if args.outdir.endswith('/') else args.outdir + '/' + args.filename
    else:
        outfile = None
        
    if args.bin_range is None:
        if args.mean_quality_scores:
            args.bin_range = '0,16,1'
        else:
            args.bin_range = '0,32,1'
    minrange, maxrange, step = (int(e) for e in args.bin_range.split(','))
    process_outdir(args.outdir)
    
    # Initialize
    quals = list()
    filesused = ''

    ## Iterate over fast5s
    for f5 in get_fast5_list(args.nfiles, args.fast5, args.random, args.randomseed, args.notarlite, filemode='r', sort=True):
        if meets_all_criteria(f5, args.readtype, args.minlen, args.maxlen, args.minq, args.maxq):
            filesused += f5.abspath + '\n'
            readtype = define_read_type(f5, args.readtype)
            if args.mean_quality_scores:
                quals.append( f5.get_mean_qscore(readtype) )
            else:
                ## TODO: make a direct way in f5class to get intquals (even if this terrible way)
                quals += [int(e) for e in (' '.join(f5.get_quals_as_int(readtype).split('\n')[1:])).split()]
 
    ##  Plot
    qualhist(quals, filename=outfile, minrange=minrange, maxrange=maxrange, step=step, density=args.density, cumulative=args.cumulative, text_only=args.no_plot)
    if args.filename is not None:
        with open(filesused_h, 'w') as fofnout:
            fofnout.write(filesused)
    else:
        sys.stderr.write(filesused)




