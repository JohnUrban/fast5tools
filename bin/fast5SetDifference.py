#!/usr/bin/env python

import h5py, os, sys
from fast5class.f5class import *
import argparse


##JOHN URBAN (2015,2016)
#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to 2 directories of fast5 file(s), return files unique to each.
Output is 2 columns format with first column given directory it is unique to and second column given file name.

Assumes that if a file appears in both directories, it has the same name.


John Urban (2015, 2016)
    """, formatter_class = argparse.RawTextHelpFormatter)



parser.add_argument('-a', '--dir1', type=str, 
                    help='''Path to directory#1 with fast5s.''')
parser.add_argument('-b', '--dir2', type=str,
                    help='''Path to directory#2 with fast5s.''')


args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################


#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    A = Fast5List(args.dir1)
    B = Fast5List(args.dir2)
    abase = set(A.get_basenames())
    bbase = set(B.get_basenames())
    adir = list(set(A.get_dirnames()))
    bdir = list(set(B.get_dirnames()))
    assert len(adir) == 1 and len(bdir) == 1
    
    u2A = list(abase.difference(bbase))
    u2B = list(bbase.difference(abase))

    for e in u2A:
        print adir[0] + "\t" + e
    for e in u2B:
        print adir[0] + "\t" + e
    
