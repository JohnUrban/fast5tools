#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import f5class
from f5ops import *
from f5tableops import *
import argparse
from glob import glob
from collections import defaultdict

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 table file (from fast5stats.py), return summary of indicated columns.

John Urban (2015)

1 = base_info_name
2 = molecule length
3 = has complement
4 = has 2d
5 = 2d seq len
6 = template seq len
7 = complement seq len
8 = 2d mean q score
9 = template mean q score
10 = complement mean q score
11 = num input events
12 = num template events
13 = num complement events
14 = num called template events
15 = num called complement events
16 = num skips in template
17 = num skips in complement
18 = fast5 filename (path as given)
19 = fast5 filename (absolute path)


    """, formatter_class = argparse.RawTextHelpFormatter)

parser.add_argument('table', metavar='table', nargs='+',
                   type= str, 
                   help='''Paths to one or as many table files as you want. Tables need to have same design -- i.e. same columns.''')
parser.add_argument('-c', '--custom', type=str, 
                    help='''Provide comma-separated list of integers corresponding to what information is in table and order it appears in.''')
parser.add_argument('-a', '--all', action="store_true", default=False,
                    help='''The "All" option (in fast5stats) was used to generate table.''')
parser.add_argument('-s', '--standard', action="store_true", default=False,
                    help='''The "standard" option (specifying options 1-17) was used in fast5stats to generate the table.''')

parser.add_argument('-d', '--delimiter', type=str, default='\t',
                    help='''Provide delimiter used to generate table. Default: tab.
For new line, type 'newline'.
Tab is default, but can also specify with "tab"''')

parser.add_argument('-w', '--which', type=str, default="",
                    help='''Specify which columns of incoming table to summarize by comma-separated integers. Default is that all columns are summarized.''')



args = parser.parse_args()


##VARIABLES
#f5cmds is the list of cmds used (in order) for data in each col of table
#col2f5cmd holds info in table column --> command
#col2list has col num as key and list of variables encountered as value

#IDEA
# read in lines and build col2list
# use col2f5cmd to know exactly what to do with given list

#################################################
## deal with some of the arguments
#################################################
num_f5cmds = len(f5fxn.keys())
f5cmds=[]
col2list = defaultdict(list)

if args.custom:
    f5cmds = [int(e) for e in args.custom.strip().split(",")]
    col2f5cmd = {k:f5cmd[k-1] for k in range(1,len(f5cmds)+1)} # 1-based col
    col2f5cmd = {k:f5cmd[k] for k in range(len(f5cmds))} # 0-based col
if args.standard:
    f5cmds = range(1,18)
    col2f5cmd = {k:k for k in f5cmds} # 1-based col
    col2f5cmd = {k-1:k for k in f5cmds} # 0-based col
if args.all:
    f5cmds = range(1,num_f5cmds+1)
    col2f5cmd = {k:k for k in f5cmds} # 1-based col
    col2f5cmd = {k-1:k for k in f5cmds} # 0-based col

if args.delimiter == "newline":
    args.delimiter = "\n"
elif args.delimiter == "tab":
    args.delimiter = "\t"


if not args.which:
    which = f5cmds
else:
    which = [int(e) for e in args.which.strip().split(",")]


#################################################
#### FUNCTIONS @@@@@@@@@@@@
#################################################





#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    for table in args.table:
        for line in open(table):
            line = line.strip().split()
            for i in range(num_f5cmds):
                try:
                    col2list[i].append( float(line[i]) )
                except:
                    pass

    f5summary[3](col2list[2])
    f5summary[4](col2list[3])
    
    f5summary[2](col2list[1],x=[25,50,75])
    f5summary[5](col2list[4],x=[25,50,75])
    f5summary[6](col2list[5],x=[25,50,75])
    f5summary[7](col2list[6],x=[25,50,75])

    f5summary[8](col2list[7])
    f5summary[9](col2list[8])
    f5summary[10](col2list[9])

