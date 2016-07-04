#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
from fast5tools.f5class import *
from fast5tools.f5ops import *
import argparse
from glob import glob

## JOHN URBAN (2015,2016)
#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s, return desired tab-delimited information from each fast5 file.

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

John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-c', '--custom', type=str, 
                    help='''Provide comma-separated list of integers corresponding to what information to return''')
parser.add_argument('-a', '--all', action="store_true", default=False,
                    help='''Return all information.''')
parser.add_argument('-s', '--standard', action="store_true", default=False,
                    help='''Return all information from options 1-17.''')

parser.add_argument('-d', '--delimiter', type=str, default='\t',
                    help='''Provide choice of delimiter. Default: tab.
For new line, type 'newline'.
Tab is default, but can also specify with "tab"''')

parser.add_argument('-e', '--errfile', type=str, default="",
                    help='''All successfully base-called files can be opened (not corrupt) and have template information.
Files that are corrupt or that do not have template information are normally listed in stderr.
Use this flag to re-direct filenames of error containing files and were therefore not reported in the stats output
to an error file.
Error output for no template is: abs path to file, number events, Log information found.
Error output for corrupt file that cannot be opened: abs path to file, -1, errmsg
''')
parser.add_argument('--verbose', type=str, default=False,
                    help='''Spit out information to progress file or stderr. Specify filename or 'stderr'. ''')


args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
f5cmds = set([])
num_f5cmds = len(f5fxn.keys())
if args.custom:
    f5cmds = f5cmds.union([int(e) for e in args.custom.strip().split(",")])
if args.standard:
    f5cmds = f5cmds.union(range(1,18))
if args.all:
    f5cmds = f5cmds.union(range(1,num_f5cmds+1))

if args.delimiter == "newline":
    args.delimiter = "\n"
elif args.delimiter == "tab":
    args.delimiter = "\t"


if not args.errfile:
    errfile = sys.stderr
else:
    errfile = open(args.errfile, 'w')

## rename fast5 to table
## allow different delimiters -- i.e. any delimiter (e.g. comma, newline, etc)

def get_fast5_stats(f5cmds, f5, delim="\t"):
    return (delim).join([ str(f5fxn[f5cmd](f5)) for f5cmd in f5cmds ])

def get_error_fast5_stats(f5, delim="\t"):
    return (delim).join([ str(f5fxn[f5cmd](f5)) for f5cmd in [19,11,20] ])

#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################



if __name__ == "__main__":
    if args.verbose:
        if args.verbose == 'stderr':
            err = sys.stderr
	else:
            err = open(args.verbose,'w')
    for f5 in Fast5List(args.fast5):
	if args.verbose:
	    err.write(f5.filename + "\n")
        if f5.is_not_corrupt() and f5.is_nonempty():
	    if f5.has_reads():
		print get_fast5_stats(f5cmds, f5, args.delimiter)
	    else:
		errfile.write( get_error_fast5_stats(f5, args.delimiter) + "\n" )
	else:
	    errfile.write( f5fxn[19](f5) + "\t" + str(-1) + "\tCould-not-open,maybe-corrupt-or-empty.\n" )

if args.errfile:
    errfile.close()
if args.verbose:
    if args.verbose != "stderr":
        err.close()


