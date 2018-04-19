#!/usr/bin/env python2.7

import h5py, os, sys
from fast5tools.f5class import *
from fast5tools.f5ops import *
from fast5tools.helperops import process_outdir, process_filesused
from fast5tools.plotops import update_qualpos, qualposplot
import argparse
from glob import glob
from collections import defaultdict
import numpy as np
from random import shuffle, seed

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s,
give box plot of the distribution of quality scores over each position.
Also reports all files used to generate boxplot.
Even when you aim it at a bunch of directories and try to use all, not all will be used due corruption, etc.
The plot can then be re-created by passing the output file-of-file-names (.fofn) file.

Note that quality tends to stay the same across nanopore read lengths.
This type of plot can create an illusion that quality drops with length.
However, this is b/c there tend to be very long reads that have low quality across their lengths.
If you filter for a minimum mean quality score (Q), this effect typically goes away.

This also allows you to plot Z-scores over each position to mitigate that effect.
It calculates Z-scores for the quals each read separate such that the resulting qualscores over
each position reflect what that position looked like relative to all positions in that read.

John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

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
Examples: qual_v_pos.jpg, qual_v_pos.png, qual_v_pos.png.
Default: None
Default will have a window pop up with the plot instead. The file can be saved from there too.
If you choose this option, the filesused will be reported to stderr.
If a filename is given, filesused will be reported in a similarly named file ending with .filesused.fofn''')


parser.add_argument('-r', '--readtype', default="template",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', 'all'.
Default: template.''')

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


def staypos_fxn(f5, args):
    ## assumes basecalled
    ## f5 is hdf5 file connection
    name = get_basename(f5)
    if name.endswith("_strand"):
        name = name[:-6]
    events = store_template_events(f5, basecalled=True, getbasecallinfo=True)
    stay_locations_BED(events, name+"template")
    if has_complement(f5):
        events = store_complement_events(f5, basecalled=True, getbasecallinfo=True)
        stay_locations_BED(events, name+"complement")



def run(parser, args):
    for fast5 in Fast5FileSet(args.files):
        f5 = fast5.hdf5file
        staypos_fxn(f5, args)
        fast5.close()

if __name__ == "__main__":
    
    # Process Args
    args.outdir = process_outdir(args.outdir) ## adds / if needed
    outfile = args.outdir + args.filename if (args.filename is not None) else None
        
    # Initialize
    filesused = ''


    ## Iterate over fast5s
    for f5 in  Fast5List(args.fast5, keep_tar_footprint_small=(not args.notarlite), filemode='r', downsample=args.nfiles, random=args.random, randomseed=args.randomseed):
        if meets_all_criteria(f5, args.readtype, args.minlen, args.maxlen, args.minq, args.maxq):
            filesused += f5.abspath + '\n'
            readtype = define_read_type(f5, args.readtype)
            ## TODO: make a direct way in f5class to get intquals (even if this terrible way)
            abspath = f5.abspath
            name = f5.get_pore_info_name(readtype)
            #pos = f5.map_stay_events_to_read(readtype)
            pos = f5.map_stay_event_coverage_in_read(readtype)
            for e in sorted(pos.keys()):
                print '\t'.join([str(e) for e in [name, e, e+1, pos[e], abspath]])



    ## Files used
    #process_filesused(trigger=args.filename, filesused=filesused, outdir=args.outdir)


