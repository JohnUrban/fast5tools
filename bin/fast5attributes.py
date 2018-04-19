#!/usr/bin/env python2.7

import h5py, os, sys
from fast5tools.f5class import *
from fast5tools.f5ops import *
import argparse
from glob import glob
from random import shuffle

#import poreminion

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s,
grab a lot of information from each.

Look at the basename.fast5info.txt output.

John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')


parser.add_argument('-n', '--nfiles', type=int, default=1,
                    help = '''Often when this command is run, looking at a single file will suffice.
Therefore, this defaults to 1, and looks at the first file in the FileList. Aim this script at a specific file for that file's contents.
Adjust this number to get info from the first N files.
Use --random to select N at random from the list.
Set this to <= 0 to get info for all fast5 files -- not typically needed/done.''')

parser.add_argument('-r', '--random', action='store_true', default=False,
                    help = '''Randomize what files are looked at.''')

parser.add_argument('-o', '--outdir', type=str, default="./",
                    help = '''....''')

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


def run_get_attributes(f5, args):
    getattrscript = '/'.join(sys.argv[0].split('/')[:-1]) + '/../fast5tools/scripts/getAttributes.sh'
    outfile = args.outdir + f5.filebasename + '.fast5info.txt'
    f5file = f5.abspath
    f5.close()
    cmd = getattrscript + " " + f5file + " > " + outfile + ' 2>/dev/null'
    #sys.stderr.write(cmd + '\n')
    os.system(cmd)
    
if __name__ == "__main__":
    # Process Args
    if not args.outdir.endswith('/'):
        args.outdir += '/'
    if not os.path.exists(args.outdir):
        os.system('mkdir ' + args.outdir)
    if args.nfiles <= 0:
        args.nfiles = len(args.fast5)
    if args.random:
        shuffle(args.fast5) ## This only shuffles initial targets, not final files -- but can help simplify target expansion

    # Downsample as necessary
    args.fast5 = args.fast5[:args.nfiles] ## This only shrinks number of targets to simplify target expansion

    # Expand targets to get initial list
    f5list = Fast5List(args.fast5, keep_tar_footprint_small=(not args.notarlite), filemode='r')

    ## Iterate over fast5s
    for f5 in f5list.get_sample(n=args.nfiles, random=args.random, sort=True): ## shuffling and downsampling actually happen here on indiv fast5s
        #run_get_attributes(f5, args)
        outfile = args.outdir + f5.filebasename + '.fast5info.txt'
        with open(outfile, 'w') as f:
            f.write(f5.get_all_attributes())

