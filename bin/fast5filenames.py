#!/usr/bin/env python
## JOHN URBAN (2015,2016)

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import fast5tools.f5class
import argparse
from glob import glob

parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s, return FILENAME for all fast5s found.

For example, use to create file of file names (.fofn).

Right now, if tars are included, it reports the filename as "tarball.tar{.gz}/filename-in-tar".


John Urban (2015, 2016)
    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')



args = parser.parse_args()


if __name__ == "__main__":
    for f in f5class.Fast5List(args.fast5, tar_filenames_only=True):
        print f

