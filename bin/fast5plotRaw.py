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
from matplotlib import collections  as mc
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

parser.add_argument('-P', '--prefix', type=str, default=None,
                    help='''Prefix for output file. Default: prefix of fast5 + .rawsignal.(extension)''')


parser.add_argument('--ext', type=str, default='pdf',
                    help='''Provide extension (png|pdf|jpg|etc). Default: pdf''')

parser.add_argument('--cex', type=float, default=0.2,
                    help='''Cex point size for scatter. Default=0.2''')

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


parser.add_argument('-S', '--segments', action='store_true', default=False,
                    help='''This plots the input event segmentation mean +/- 2stdv information.
                            It interprets -b and -e as event indexes. Raw data corresponding to these events is shown.''')

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


def plotter(x,y,start,end,d, col='b', plottype='scatter', events=None):
    ''' d = dict'''
    plt.title(args.title)
    plt.xlabel('Raw Data Point Number')
    plt.ylabel('Level')
    if plottype == "scatter":
        plt.scatter(d[x][start:end], d[y][start:end], c=col, edgecolor='k', s=0.1)
    elif plottype == "line":
        plt.plot(d[x][start:end], d[y][start:end], c=col,zorder=1)

def segmentedRawPlotter(events, raw, start, end, col, title):
    rstart = sum(events['length'][0:start])
    rend = rstart + sum(events['length'][start:end])
    rawplot = raw[rstart:rend]

    sumlengths = np.cumsum(events['length'])
    y=events['mean'][start:end]
    sdup = events['mean'][start:end] + 2*events['stdv'][start:end]
    sddown = events['mean'][start:end] - 2*events['stdv'][start:end]
    if start == 0:
        x0 = [0] + list(sumlengths[start:(end-1)]) #- sumlengths[start]
    elif start > 0:
        x0 = sumlengths[(start-1):(end-1)]
    x1=sumlengths[(start):(end)] #- sumlengths[start]
    ##print len(y), len(sdup), len(sddown), len(x0), len(x1)
    ##assert len(y) == len(sdup) == len(sddown) == len(x0) == len(x1)
    mus = mc.LineCollection([[(x0[i], y[i]), (x1[i], y[i])] for i in range(len(y))], colors='r', linewidths=1)
    up = mc.LineCollection([[(x0[i], sdup[i]), (x1[i], sdup[i])] for i in range(len(y))], colors='b', linewidths=0.1)
    down = mc.LineCollection([[(x0[i], sddown[i]), (x1[i], sddown[i])] for i in range(len(y))], colors='b', linewidths=0.1)
    fig, ax = matplotlib.pyplot.subplots()
    plt.title(title)
    plt.xlabel('Raw Data Point Number')
    plt.ylabel('Level')
    #w=white, k=black, b=blue, r=red, 'none' = nocolor
    ax.scatter(range(rstart,rend), rawplot, c='none', edgecolor='k', s=0.5, linewidths=0.1, marker='o')
    ax.add_collection(mus)
    ax.add_collection(up)
    ax.add_collection(down)
    ax.autoscale()
##    ax.margins(0.1)
    

#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

##rstart = sum(events['length'][0:start])
##rend = rstart + sum(events['length'][start:end])
##plot x=range(sum(events['length'][start:end])), y=raw[rstart:rend]
##sumlengths = np.cumsum(events['length'][0:end])
##y=events['mean'][start:end]
##sdup = events['mean'] + 2*events['stdv'][start:end]
##sddown = events['mean'] - 2*events['stdv'][start:end]
##segment x0=sumlengths[(start-1):end], x1=sumlengths(start:(end+1)), y0=y, y1=y
##segment x0=sumlengths[(start-1):end], x1=sumlengths(start:(end+1)), y0=sdup, y1=sdup
##segment x0=sumlengths[(start-1):end], x1=sumlengths(start:(end+1)), y0=sddown, y1=sddown


if __name__ == "__main__":
    if args.title is None:
        args.title = "" 
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=(not args.notarlite)):
        if args.prefix is None:
            if args.segments:
                args.prefix = f5.filebasename + '.rawsignalsegments.'
            else:
                args.prefix = f5.filebasename + '.rawsignal.'
        outname = args.prefix + args.ext
        outname = '.'.join([e for e in outname.split('.') if e])
        x = []
        y = []
        raw = f5.get_raw_signal()
        nraw = range(f5.get_raw_duration())
        if args.segments:
            events = f5.get_events('input')
            
        else:
            events = None
        ## will need to process those with "-"
        start = args.start
        if args.end is None:
            if args.segments:
                end = len(events['mean'])
            else:
                end =  f5.get_raw_duration()
        else:
            end = args.end

        ##PLOTTING
        if args.segments:
            segmentedRawPlotter(events, raw, start, end, col=args.color, title=args.title)
        else:
            plotter('num', 'raw', start, end, {'num':nraw, 'raw':raw}, col=args.color, plottype=args.plottype, events=events)

##        if args.outdir:
##            plt.savefig(os.path.join(args.outdir, '%s.png' % f5.filebasename))
        plt.savefig(outname)
##        plt.show()

