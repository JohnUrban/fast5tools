#!/usr/bin/env python2.7

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

Given path(s) to fast5 file(s) and/or directories of tombo-resquiggled fast5s,

return desired features from Tombo.

For files that are corrupt or empty, for now it silently skips them.
As an alternative, fast5stats will tell you all files skipped (in stderr or to specified file).

John Urban (2015, 2016, 2017)


    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser_get = parser.add_mutually_exclusive_group(required=True)
parser_get.add_argument('--getmap', action='store_true', default=False, help='''Get map positions.''')
parser_get.add_argument('--getevents', action='store_true', default=False, help='''Get tombo events.''')
parser_get.add_argument('--getparams', action='store_true', default=False, help='''Get stored tombo parameters.''')

##parser.add_argument('-r', '--readtype', default="mol",
##                   type= str, 
##                   help='''Choose type of fasta to get.
##Choices: 'template', 'complement', '2d', 'molecule', 'all', 'MoleQual'.
##Default: molecule.
##There is no need to write full word for options - can do: t, c, 2, m, a, M.
##Molecule returns single fasta for each fast5 by following rules:
##if 2d present, return 2d.
##elif complement present with no 2d, return longer of template or complement.
##elif only template present, return template.
##'MoleQual' is similar to molecule.
##It differs only in choosing between template and complement when a 2D is not present.
##Instead of choosing the longer one, it chooses the one with a higher quality mean quality score.''')


##parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')
##
##parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')
##
##parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')
##
##parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
##Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')

parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
Tarlite may become the default method after some testing if it performs at similar speeds.''')




args = parser.parse_args()


                
            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":

    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=args.tarlite):
        if f5.is_not_corrupt() and f5.is_nonempty and f5.tombo_exists():
            if args.getmap and f5.tombo_aln_exists():
                print f5.get_tombo_map_position_string() + '\t' + f5.filename
            elif args.getevents and f5.tombo_events_exist():
                print f5.get_tombo_events_string() ## this needs to write to file when more than one f5
            elif args.getparams and f5.tombo_params_exist():
                if f5.tombo_successful():
                    print f5.filename + '\t' + f5.get_tombo_parameter_string() + '\t' + f5.get_tombo_parameter_header_string(',')
                else:
                    print f5.filename + '\t' +  ('_').join(f5.get_tombo_parameter_string().split()) + '\t' + ('\t').join([str(e) for e in ['*']*7])  + '\t' + f5.get_tombo_parameter_header_string(',')
        else:
            if not f5.tombo_exists():
                print "No_Tombo_traces_found...\t" + f5.filename


