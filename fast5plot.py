#!/usr/bin/env python
## JOHN URBAN (2015,2016)

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import f5class
from f5ops import *
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

Given path(s) to fast5 file(s) and/or directories of fast5s, return desired plot given x and y.


John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-x', '--x', type=int, default=5,
                    help='''Provide integer corresponding to what information is on x-axis.''')

parser.add_argument('-y', '--y', type=int, default=8,
                    help='''Provide integer corresponding to what information is on y-axis.''')

parser.add_argument('-t', '--title', type=str, default="2D Seq len vs Mean Q-score",
                    help='''Provide title.''')

args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################

num_f5cmds = len(f5fxn.keys())



def get_fast5_data(f5cmd, f5):
    try:
        return float(f5fxn[f5cmd](f5))
    except:
        return None
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################


##TODO:
## Also make plotting from fast5totable files
if __name__ == "__main__":
    x = []
    y = []
    for f5 in f5class.Fast5List(args.fast5):
        x.append( get_fast5_data(args.x, f5) )
        y.append( get_fast5_data(args.y, f5) )
    print x
    print y
    ## will need to process those with "-"
    plt.title(args.title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.scatter(x,y)
    plt.show()
    

