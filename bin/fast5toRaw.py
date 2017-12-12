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

Given path(s) to fast5 file(s) and/or directories of fast5s, return raw data points.

Be careful. This will only work on recent fast5 file versions that provide raw data.


John Urban (2015, 2016, 2017)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')


parser.add_argument('-s', '--singleline', action='store_true', default=False, help='''The default output is one data point per line like a single-column table.
If you'd like a single-row table instead, use this option.
Use this option in conjunction with -d/--delimiter to determine the delimiter between datapoints.
The default delimiter is a space.''')

parser.add_argument('-d', '--delimiter', type=str, default=" ", help='''Only meaningful when used with -s/--singleline or -S/--segmented.
Default is a space. Write the word "tab" if you want tabs. Otherwise, whatever string is provided is used.''')

parser.add_argument('-S', '--segmented', action='store_true', default=False, help='''The default output is one data point per line like a single-column table.
If you'd like the data to be segmented as it is in the events.
It will be one event per line.''')

parser.add_argument('-E', '--eventdelimiter', type=str, default="\n", help='''Only meaningful when used with -S/--segmented.
Default is a newline. Write the word "tab" if you want tabs. Otherwise, whatever string is provided is used.''')

parser.add_argument('-N', '--noclips', action='store_true', default=False, help='''Only meaningful when used with -S/--segmented.
Raw data often has 5' and or 3' datapoints that are not part of events.
These are automatically included as first and last line (if they exist).
This says not to include them.''')

parser.add_argument('-x', '--eventstats', action='store_true', default=False, help='''
Will give Event number, mean, stdv, start, length on raw signal -- segmented according to events.
Useful to spot check events....''')

parser.add_argument('-M', '--mediannorm', action='store_true', default=False, help='''
Median normalize the signal before returning it or operating on it....
Each data point = 100*dp/median''')

parser.add_argument('-o', '--outdir', type=str, default="",
                    help = '''If single fast5 specified, it will be reported to stdout.
If multiple fast5s are specified, they will be saved to files in the working dir by default.
This flag allows you to specify a different output directory.
Filenames will be the the name of the fast5 file with .rawsignal.txt appended.''')

parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
Tarlite may become the default method after some testing if it performs at similar speeds.''')



args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################


    
if args.outdir:
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    if args.outdir[-1] != "/":
        args.outdir += "/"

if not args.singleline and not args.segmented:
    delimiter = "\n"
else:
    if args.delimiter == 'tab':
        delimiter = '\t'
    else:
        delimiter = args.delimiter

if args.eventdelimiter == 'tab':
    eventdelim = '\t'
else:
    eventdelim = args.eventdelimiter

clips = not args.noclips

            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    f5list = Fast5List(args.fast5, keep_tar_footprint_small=args.tarlite)

    if args.segmented:
        if len(f5list) == 1: ## if only one f5, print to stdout
            for f5 in f5list:
                print f5.get_segmented_raw_signal_string(includeclips=clips, datadelim=delimiter, eventdelim=eventdelim, readtype='template', median_normalized=args.mediannorm)
        elif len(f5list) > 1: ## if more than one, print each f5 to own text file
            for f5 in f5list:
                out = open(args.outdir + f5.filebasename + "." + args.readtype + ".rawsignal.txt", 'w')
                out.write(f5.get_segmented_raw_signal_string(includeclips=clips, datadelim=datadelim, eventdelim=delimiter, readtype='template', median_normalized=args.mediannorm))
                out.close()
    elif args.eventstats:
        if len(f5list) == 1: ## if only one f5, print to stdout
            for f5 in f5list:
                print f5.get_segmented_raw_signal_stats_string(includeclips=clips, readtype='template', median_normalized=args.mediannorm)
        elif len(f5list) > 1: ## if more than one, print each f5 to own text file
            for f5 in f5list:
                out = open(args.outdir + f5.filebasename + "." + args.readtype + ".rawsignal.txt", 'w')
                out.write('')
                out.close()
    else:
        if len(f5list) == 1: ## if only one f5, print to stdout
            for f5 in f5list:
                print f5.get_raw_signal_string(delimiter, median_normalized=args.mediannorm)
        elif len(f5list) > 1: ## if more than one, print each f5 to own text file
            for f5 in f5list:
                out = open(args.outdir + f5.filebasename + "." + args.readtype + ".rawsignal.txt", 'w')
                out.write(f5.get_raw_signal_string(delimiter, median_normalized=args.mediannorm))
                out.close()
                





