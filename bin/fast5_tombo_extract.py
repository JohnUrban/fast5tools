#!/usr/bin/env python2.7

import h5py, os, sys, tarfile
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
parser_get.add_argument('--getgenomicevents', action='store_true', default=False, help='''Get genome-anchored tombo events.''')
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


##TODO
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
Filenames will be the the name of the fast5 file with .tombo_events.txt appended.''')

parser.add_argument('--targzout', type=str, default=False,
                    help=''' For now, only relevant to --getevents with >1 fast5 file as input.
                    This will extract the events into a file, then add the file to a growing tarchive,
                    followed by deleting the file. This is a good option to prevent exceeding file number quotas, etc.
                    Provide a tarchive name - this script insists on using .tar.gz at the end, and will add it if absent.
                    This will put the tarchive into the specified outdir whether or not the outdir is included as part of the given name here.''')

args = parser.parse_args()


                

def tarpit(tarchive, fn, arcname=None):
    ''' Assumes tarhive is an already existing tarfile object for 'w|gz'
        fn is filename to add.
        Will add file to the arhive, then delete the original...
        Be careful.
        arcname is an alternate name to give the file inside the tarchive.'''
    if arcname is None:
        arcname = fn
    tarchive.add(name=fn, arcname=arcname)
    os.remove(fn)

#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    if args.outdir:
        if not os.path.isdir(args.outdir):
            os.mkdir(args.outdir)
        if args.outdir[-1] != "/":
            args.outdir += "/"

    if args.targzout:
        if not args.targzout.endswith('.tar.gz'):
            args.targzout += '.tar.gz'
        tarpath = args.targzout
        if not tarpath.startswith(args.outdir):
            tarpath = args.outdir + tarpath
        tarchive = tarfile.open(tarpath, 'w|gz')

    f5list = Fast5List(args.fast5, keep_tar_footprint_small=(not args.notarlite))
    f5listlen = len(f5list)
    only_one_file = f5listlen == 1

##    if f5listlen == 1 or args.getmap or args.getparams: ## If only 1 file to extract events from OR getting map positions or parameters for reads (which have readnames as a column in output), report to stdout 
    if args.getmap or args.getparams: 
        for f5 in f5list:
            if f5.is_not_corrupt() and f5.is_nonempty and f5.tombo_exists():
                if args.getmap and f5.tombo_aln_exists():
                    print f5.get_tombo_map_position_string() + '\t' + f5.filename
                    ## No longer want to allow different behavior for 1 file vs many since this may have unintended effects
##                elif (args.getevents or args.getgenomicevents) and f5.tombo_events_exist():
##                    if args.getevents:
####                        if args.
##                        print f5.get_tombo_events_string() ## this needs to write to file when more than one f5
##                    elif args.getgenomicevents:
##                        print f5.get_tombo_genomic_events_string()
                elif args.getparams and f5.tombo_params_exist():
                    if f5.tombo_successful():
                        print f5.filename + '\t' + f5.get_tombo_parameter_string() + '\t' + f5.get_tombo_parameter_header_string(',')
                    else:
                        print f5.filename + '\t' +  ('_').join(f5.get_tombo_parameter_string().split()) + '\t' + ('\t').join([str(e) for e in ['*']*7])  + '\t' + f5.get_tombo_parameter_header_string(',')
            else:
                if not f5.tombo_exists():
                    print "No_Tombo_traces_found...\t" + f5.filename
##    elif f5listlen > 1:
    elif (args.getevents or args.getgenomicevents): #getevents or getgenomicevents
        ## Only matters when extracting events
        for f5 in f5list:
            if f5.tombo_events_exist():
                if args.header:
                    if args.headertable:
                        if args.getevents:
                            print ("\t").join([f5.filename, f5.get_tombo_events_header()])
                        else: #assume genomic
                            print ("\t").join([f5.filename, f5.get_tombo_genomic_events_header()])
                    else:
                        fn = args.outdir + f5.filebasename + ".tombo_events_headers.txt"
                        out = open(fn, 'w')
                        if args.getevents:
                            out.write(f5.get_tombo_events_header_string())
                        else:
                            out.write(f5.get_tombo_genomic_events_header_string())
                        out.close()
                        if args.targzout:
                            arcfn = f5.filebasename + ".tombo_events_headers.txt"
                            tarpit(tarchive, fn, arcfn)
                else:
                    fn = args.outdir + f5.filebasename + ".tombo_events.txt"
                    out = open(fn, 'w')
                    if args.withheader:
                        if args.getevents:
                            out.write(f5.get_tombo_events_header_string() + "\n")
                        else:
                            out.write(f5.get_tombo_genomic_events_header_string() + "\n")
                    if args.getevents:
                        out.write(f5.get_tombo_events_string())
                    else:
                        out.write(f5.get_tombo_genomic_events_string())
                    out.close()
                    if args.targzout:
                        arcfn = f5.filebasename + ".tombo_events.txt"
                        tarpit(tarchive, fn, arcfn)

    if args.targzout:
        tarchive.close()


## TODO - do reference mapped events by getting map position and arithmeticing the start and end across events...
