#!/usr/bin/env python
## JOHN URBAN (2015,2016)

import h5py, os, sys

from fast5tools.f5class import *
from fast5tools.hmm_class import *
from fast5tools.parsedEventsClass import *

import argparse
from glob import glob
import matplotlib
## may need following line for remote jobs (e.g. submitting batch scripts)
##matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
##from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file, return plot of events in various ways.




John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')


parser.add_argument('-t', '--title', type=str, default=None,
                    help='''Provide title.''')

parser.add_argument('-p', '--plottype', type=str, default='scatter',
                    help='''Specify: scatter, line ....  Default: scatter.''')

parser.add_argument('-c', '--color', type=str, default='b',
                    help='''Colors as understood by pyplot... default: b''')


parser.add_argument('-s', '--start', type=int, default=0,
                    help='''Start index of data points to show. Default: 0. i.e. from the beginning.''')

parser.add_argument('-e', '--end', type=int, default=None,
                    help='''End index of data points to show. Default: number of data points. i.e. to the end.''')

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
## deal with some of the arguments
#################################################


def plotter(x,y,start,end,d, col='b', plottype='scatter'):
    ''' d = dict'''
    if plottype == "scatter":
        plt.scatter(d[x][start:end], d[y][start:end], c=col, edgecolor='k')
    elif plottype == "line":
        plt.plot(d[x][start:end], d[y][start:end], c=col)



#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################



if __name__ == "__main__":
    if args.title is None:
        args.title = "" 
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=(not args.notarlite)):
        x = []
        y = []
        raw = f5.get_raw_signal()
        nraw = range(f5.get_raw_duration())
        ## will need to process those with "-"
        plt.title(args.title)
        plt.xlabel('Raw Data Point Number')
        plt.ylabel('Level')
        if args.end is None:
            end =  f5.get_raw_duration()
        else:
            end = args.end
        plotter('num', 'raw', args.start, end, {'num':nraw, 'raw':raw}, col=args.color, plottype=args.plottype)

##        if args.outdir:
##            plt.savefig(os.path.join(args.outdir, '%s.png' % f5.filebasename))

        plt.show()

