#!/usr/bin/env python2.7
import argparse
from collections import defaultdict
from fast5tools.samclass import *
from fast5tools.samops import *


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Given a SAM file (with F5:Z: info attached) that is sorted by read name:
        - get the alignment or set of splitread alignments for each read
        - determine most likely genomic region read came from (assuming no structural variation)
            - if one alignment, assume it comes from there
            - if multiple alignments,
                check for overlap of their individual genomic windows (alignment adjusted for clipping on each side + flank/buffer)
                if no merges,
                    use majority or longest alignment (majority is longest alignment that also meets a majority threshold)
                if there is a single merge -- i.e. they all come from same genomic region (and perhaps required to be ordered and stranded - see options) -
                    use merged result from merged genomic windows
                if there is 1 or more merges (but still more than 1 genomic region)
                    see if longest merge has a 'spanning alignment' longer than longest/majority alignment
                    if so use that, if not use the longest/majority alignment
                - report on alignments and merges in all cases
        - get coordinates for a window that surrounds that chosen genomic region
            - this is the chosen genomic window for that read
            - coordinates for genomic window should be proportional to read length + some extra buffering/flanking sequence
        - print out gw coordinates, notes on choice, F5 info, and perhaps genomic sequence chosen


    flank=0.25, merge_dist=0, majority=0.5, require_order=False, require_strand=False, reference=False

        flank = Add buffer/flank lengths to each side of a genomic window in two ways:
            (1) int > 1 adds/subtracts that int.
            (2) float [0,1] adds/subtracts that proportion of read length
                    NOTE: 1 defaults to 100% of read length, not 1 bp
            merge_dist:
                allows a gap up to d between intervals to still be an overlap - default 0
            majority
                threshold to exceed to be considered a majority.
            require_order
                when True, alignments must be ordered as they appear in the read to be considered a valid merge.
                Defaults to False as noisy alignments could easily break this. Status is reported in output anyway.
            require_strand
                when True, alignments must ALL be on the same strand to be considered a valid merge.
                Defaults to False as noisy alignments could easily break this. Status is reported in output anyway.
                

   """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--sam', '-s',
                   type= str, default=False,
                   help='''Input file in SAM format.''')

## FOR NOW, MUST BE SAM -- NOT BAM -- but can be STDIN SAM
##parser_input.add_argument('--bam', '-b',
##                   type= str, default=False,
##                   help='''Input file in BAM format.''')

parser.add_argument('--flank', '-f', type=float, default=0.25,
                    help=''' ''')

parser.add_argument('--merge_dist', '-m', type=int, default=0,
                    help=''' ''')

parser.add_argument('--majority', '-M', type=float, default=0.5,
                    help=''' ''')

parser.add_argument('--require_order', '-ro', action='store_true', default=False,
                    help=''' ''')

parser.add_argument('--require_strand', '-rs', action='store_true', default=False,
                    help=''' ''')

parser.add_argument('--reference', '-r', type=str, default=False,
                    help=''' Path to reference genome file to be used to extract sequences corresponding to genomic windows identified.
                            Optional. Sequences will be tagged on to an additional end column if provided.''')


parser.add_argument('--getF5info', '-f5', action='store_true', default=False,
                    help='''Return F5:Z: field from fast5tools in output.
This is from extracting fasta/fastq using fast5tofastx.py with --comments and --samflag''')

parser.add_argument('--getBCinfo', '-BC', action='store_true', default=False,
                    help=''' Return BC:Z: field from fast5tools in output.
This is from creating fasta/fastq from output of fast5_sw_bardecoder.py specified with --sequence/--quals,
and merging all desired barcode info into string following BC:Z:''')

args = parser.parse_args()




get_genomic_windows(samfilepath=args.sam, flank=args.flank, merge_dist=args.merge_dist, majority=args.majority, require_order=args.require_order, require_strand=args.require_strand, reference=args.reference, getF5field=args.getF5info, getBCfield=args.getBCinfo)



