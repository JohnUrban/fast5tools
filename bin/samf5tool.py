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








def alignment_summary(sam):
    '''sam is Sam object'''
    for read in sam:
        print ("\t").join([str(e) for e in [read.get_qname_field(), read.get_rname_field(), read.get_pos_field(), read.get_read_len(), read.get_AS_field(), read.get_SEQ_len(), read.get_SEQ_len_without_clipped_regions(), read.get_reference_aln_len(), read.get_edit_dist_field(), read.get_edit_dist_with_clipping(), read.get_clipping_dist(), read.get_fast5_field()]])

def get_genomic_windows(sam):
    ''' sam is SamSplitAlnAggregator object'''
    for read in sam:
##        read = SplitReadSamRecord(read)
        if read.has_alignments():
            ans = read.get_genomic_window(flank=0.25, merge_dist=0, majority=0.5, require_order=False, require_strand=False)
            if ans is not None:
                print ("\t").join( [read.get_read_name()]+[str(e) for e in ans]+[read.get_fast5_info()])
                ## if reference fasta provided, then print out sequence in last field of table
                ## should any other SAM info be included?
                ## include edit dist? num matches? %id? of longest only? all alns? all_alns vs readlen?
            else:
                print ans
        else:
            pass #print to no alignments file....


def get_pct_identity_v_Q(sam):
    ''' sam is SamSplitAlnAggregator object'''
    for read in sam:
        print read.get_edit_dist_fields()
        print ("\t").join( [read.get_read_name(), read.get_fast5_info()])





##sam = Sam(args.sam)
##alignment_summary(sam)

sam = SamSplitAlnAggregator(args.sam)
##get_genomic_windows(sam)
get_pct_identity_v_Q(sam)




