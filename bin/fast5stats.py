#!/usr/bin/env python2.7

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
3 = molecule mean q score
4 = has complement
5 = has 2d
6 = 2d seq len
7 = template seq len
8 = complement seq len
9 = 2d mean q score
10 = template mean q score
11 = complement mean q score
12 = num input events
13 = num template events
14 = num complement events
15 = num called template events
16 = num called complement events
17 = num skips in template
18 = num skips in complement
19 = fast5 filename (path as given)
20 = fast5 filename (absolute path)


Note:
This was initially written in the early days of the MinION.
Then, hairpin adapters were used to give template strand, complement strand, and 2D reads.
Now only the 1D template strands are common.

John Urban (2015, 2016,2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-c', '--custom', type=str, 
                    help='''Provide comma-separated list of integers corresponding to what information to return''')
parser.add_argument('-a', '--all', action="store_true", default=False,
                    help='''Return all information. For files where complement and 2D information are relevant. Otherwise, see --template.''')
parser.add_argument('-s', '--standard', action="store_true", default=False,
                    help='''Return all information from options 1-17. For files where complement and 2D information are relevant. Otherwise, see --template''')

parser.add_argument('-t', '--template', action="store_true", default=False,
                    help='''Returns template-relevant information only: 1,7,10,13,15.
To also get path information, just use "--template --custom 19" or "--template --custom 20".
This would be the same as doing: "--custom 1,7,10,13,15,19" or "--custom 1,7,10,13,15,20"''')



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
if args.template:
    f5cmds = f5cmds.union([1,7,10,13,15])

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
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=(not args.tarlite)):
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


