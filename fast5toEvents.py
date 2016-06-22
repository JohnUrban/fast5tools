#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import f5class
from f5ops import *
import argparse
from glob import glob


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s, return table of events for each.

Be careful. The columns of events tables change with fast5 version.


Example of very early basecalled fast5 files circa October 2014:
input events:
mean, stddev, start, length (duration)

template/complement events:
mean, start, stddev, length, model_state, model_level, move, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T, 


Example of more recent (circa Nov 2015) basecalled fast5 (notice extra "weights" in template/complement):
input events (same):
mean, stddev, start, length (duration)

template/complement events:
mean, start, stddev, length, model_state, model_level, move, weights, p_model_state, mp_state, p_mp_state, p_A, p_C, p_G, p_T, 


There is also a difference between non-basecalled and base-called input events. For example, for an early comparison, see input events for basecalled file above. For a non-basecalled file, the input events look like:
start, length, mean, variance
Not only a diffent order, but variance is shown instead of std deviation.

I originally had poreminion correct for this and output the same as a basecalled file. 
Currently, fast5tools does not do that in order to be more flexible with all the changes made in fast5 versions and changes to come.
Instead (TODO), it should just optionally print a header -- and optionally convert variance to std (though that can be done with awk or something).


John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="input",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', 'input'.
Default: input.
''')

parser.add_argument('-H', '--header', default=False,
                   action='store_true', 
                   help='''Print only header and exit.''')

parser.add_argument('-W', '--withheader', default=False,
                   action='store_true', 
                   help='''Print events with header line.''')


parser.add_argument('-T', '--headertable', default=False,
                   action='store_true', 
                   help='''Instead of printing tab-separated header, print 2 columns:
1. filename. 2. comma-separated header names. NOTE: this over-rides the default behavior when more than 1 fast5 of writing a text file for each fast5.
It might make more sense to do this when all you want is to see headers of a set of different fast5s for comparison.''')

parser.add_argument('-o', '--outdir', type=str, default="",
                    help = '''If single fast5 specified, it will be reported to stdout.
If multiple fast5s are specified, they will be saved to files in the working dir by default.
This flag allows you to specify a different output directory.
Filenames will be the the name of the fast5 file with .events.txt appended.''')



args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
assert args.readtype[0] in "tci"
if args.readtype[0] == "t":
    args.readtype = "template"
elif args.readtype[0] == "c":
    args.readtype = "complement"
elif args.readtype[0] == "i":
    args.readtype = "input"
if args.outdir:
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    if args.outdir[-1] != "/":
        args.outdir += "/"



            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    f5list = f5class.Fast5List(args.fast5)
    if len(f5list) == 1: ## if only one f5, print to stdout
        for f5 in f5list:
            if f5.has_read(args.readtype):
                if args.header:
                    if args.headertable:
                        print ("\t").join([f5.filename, (",").join( e for e in f5.get_events_header(args.readtype))])
                    else:
                        print f5.get_events_header_string(args.readtype)
                else:
                    if args.withheader:
                        print f5.get_events_header_string(args.readtype)
                    print f5.get_events_string(args.readtype)
    elif len(f5list) > 1: ## if more than one, print each f5 to own text file
        for f5 in f5list:
            if f5.has_read(args.readtype):
                if args.header:
                    if args.headertable:
                        print ("\t").join([f5.filename, f5.get_events_header(args.readtype)])
                    else:
                        out = open(args.outdir + f5.filebasename + "." + args.readtype + ".eventsheader.txt", 'w')
                        out.write(f5.get_events_header_string(args.readtype))
                        out.close()
                else:
                    out = open(args.outdir + f5.filebasename + "." + args.readtype + ".events.txt", 'w')
                    if args.withheader:
                        out.write(f5.get_events_header_string(args.readtype))
                    out.write(f5.get_events_string(args.readtype))
                    out.close()

        ## TODO TEST


