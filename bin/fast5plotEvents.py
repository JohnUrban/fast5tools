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

parser.add_argument('-x', '--x', type=str, default='num',
                    help='''Provide whether x-axis should correspond to event 'num', 'mean', 'stdv', 'length', 'start'.''')

parser.add_argument('-y', '--y', type=str, default='mean',
                    help='''Provide whether y-axis should correspond to event 'num', 'mean', 'stdv', 'length', 'start'.''')

parser.add_argument('-t', '--title', type=str, default=None,
                    help='''Provide title.''')

parser.add_argument('-r', '--readtype', type=str, default='input',
                    help='''Specify: template, complement, or input. Default: input.''')

parser.add_argument('-p', '--plottype', type=str, default='scatter',
                    help='''Specify: scatter, line, ....  Default: scatter.''')

parser.add_argument('-a', '--adapters', action='store_true', default=False,
                    help='''Highlight adapter sequences. Default: False.''')

parser.add_argument('-o', '--outdir', default=None,
                    help='''Directory of where to save plots. If none given, plots will not be saved.''')

parser.add_argument('-AF', '--adapterfxn', type=int, default=1,
                    help='''When highlighting adapter sequences, find adapters with fxn 1 (older hairpins) or 2 (newer hairpins). Default: 1. Fxn2 is not well-tested - experimental.''')


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


safe_keys = ['num', 'mean', 'stdv', 'length', 'start']
assert args.x in safe_keys
assert args.y in safe_keys

names = {}
names['num'] = "Event Number"
names['mean'] = "Event Means"
names['stdv'] = "Event Standard Deviations"
names['length'] = "Event Lengths"
names['start'] = "Event Start Times"

def make_title(x,y, names):
    return names[x] + " Vs. " + names[y]
    

def plotter(x,y,start,end,events_dict, col='b'):
    if args.plottype == "scatter":
        plt.scatter(events_dict[x][start:end], events_dict[y][start:end], c=col, edgecolor='k')
    elif args.plottype == "line":
        plt.plot(events_dict[x][start:end], events_dict[y][start:end], c=col)

def parser(f5, events=None, hpfxn=1):
##    if f5.basecalling_detected(): ## TODO: just get coordinates of adapters from fast5 file
##            parsed_events = ParsedEvents(f5=f5) ## for now just filling in with adapter finding fxn
##    else: ## find coordinates 
##        parsed_events = ParsedEvents(events=events)
##    parsed_events = ParsedEvents(events=events)
    if hpfxn==1:
        halfsize=40
    elif hpfxn==2:
        halfsize=100
    parsed_events = ParsedEvents(f5=f5, hp_half_size=halfsize, hpfxn=hpfxn)
    return parsed_events

#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################



if __name__ == "__main__":
    if args.title is None:
        args.title = make_title(args.x, args.y, names=names)
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=(not args.tarlite)):
        x = []
        y = []
        events_dict = f5.get_events_dict(args.readtype)
        nevents = len(events_dict['mean'])
        events_dict['num'] = range(1,nevents+1)
        ## will need to process those with "-"
        plt.title(args.title)
        plt.xlabel(names[args.x])
        plt.ylabel(names[args.y])

        plotter(args.x, args.y, 0, nevents, events_dict)
        if args.adapters:
            parsed_events = parser(f5, hpfxn=args.adapterfxn)
            plotter(args.x, args.y, 0, parsed_events.template_start-1, events_dict, col='r')
            if parsed_events.hairpin_detected:
                plotter(args.x, args.y, parsed_events.hpstart, parsed_events.complement_start, events_dict, col='r')
        
        if args.outdir:
            plt.savefig(os.path.join(args.outdir, '%s.png' % f5.filebasename))

        plt.show()

