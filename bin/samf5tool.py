#!/usr/bin/env python2.7
import argparse
from collections import defaultdict
from fast5tools.samclass import *
##from fast5tools.f5class import *
##from fast5tools.f5ops import *


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Designed specifically for one case (described below), but may have use elsewhere.

    Use case:
    Started by mapping long reads to Ecoli genome that has been doubled to allow clean
        read mappings across the "cut site" where the circular genome was
        linearized.
    Now you want to readjust all starting positions to be inside the first copy
        of the two copy genome.
    For example, the genome size of K12 MG1655 is:
        NC_000913.3	4641652
    For all alignments that start at positions > 4641652, subtract 4641652 from them.
    
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--sam', '-s',
                   type= str, default=False,
                   help='''Input file in SAM format.''')



## FOR NOW, MUST BE SAM -- NOT BAM -- but can be STDIN SAM
##parser_input.add_argument('--bam', '-b',
##                   type= str, default=False,
##                   help='''Input file in BAM format.''')

args = parser.parse_args()






##sam = SamSplitAlnAggregator(args.sam)
##for read in sam:
##    print record

def alignment_summary(sam):
    for read in sam:
        print ("\t").join([str(e) for e in [read.get_qname_field(), read.get_rname_field(), read.get_pos_field(), read.get_read_len(), read.get_AS_field(), read.get_SEQ_len(), read.get_SEQ_len_without_clipped_regions(), read.get_reference_aln_len(), read.get_edit_dist_field(), read.get_edit_dist_with_clipping(), read.get_clipping_dist(), read.get_fast5_field()]])

sam = Sam(args.sam)

alignment_summary(sam)






