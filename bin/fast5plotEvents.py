#!/usr/bin/env python
## JOHN URBAN (2015,2016)

import h5py, os, sys

from fast5tools.f5class import *
from fast5tools.hmm_class import *
from fast5tools.adapters import *
from fast5tools.parsedEventsClass import *

import argparse
from glob import glob
import matplotlib
## may need following line for remote jobs (e.g. submitting batch scripts)
##matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
##from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

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
                    help='''Specify: scatter, ....  Default: scatter.''')

parser.add_argument('-a', '--adapters', action='store_true', default=False,
                    help='''Highlight adapter sequences. Default: False.''')

parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
Tarlite may become the default method after some testing if it performs at similar speeds.''')


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

def parser(f5, events=None):
##    if f5.basecalling_detected(): ## TODO: just get coordinates of adapters from fast5 file
##            parsed_events = ParsedEvents(f5=f5) ## for now just filling in with adapter finding fxn
##    else: ## find coordinates 
##        parsed_events = ParsedEvents(events=events)
##    parsed_events = ParsedEvents(events=events)
    parsed_events = ParsedEvents(f5=f5)
    return parsed_events

#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################



if __name__ == "__main__":
    if args.title is None:
        args.title = make_title(args.x, args.y, names=names)
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=args.tarlite):
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
            parsed_events = parser(f5)
            plotter(args.x, args.y, 0, parsed_events.template_start-1, events_dict, col='r')
            if parsed_events.hairpin_detected:
                plotter(args.x, args.y, parsed_events.hpstart, parsed_events.complement_start, events_dict, col='r')
        
        plt.show()
    

