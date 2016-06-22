#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import f5class

if __name__ == "__main__":
    f5 = f5class.Fast5(sys.argv[1])
    f5.test()
