#!/usr/bin/env python2.7
import argparse
from collections import defaultdict
from fast5tools.samclass import *



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

parser_action = parser.add_mutually_exclusive_group(required=True)
parser_action.add_argument('--replace', '-r', type=int, default=False,
                           help='''This integer replaces integer is POS field,
                                    given that your requirement is met.''')
parser_action.add_argument('--add', '-a', type=int, default=False,
                           help='''This integer is added to the POS field,
                                    given that your requirement is met.
                                    Negative integers will be subtracted.
                                    E.g. Use -4641652 for Doubled Ecoli genome of length 4641652 with limti of 4641652 (below).''')

parser.add_argument('--operator', '-op', type=str, default='gt',
                           help='''Relationship of POS field to limit argument below. Options: gt, ge, lt, le, eq, ne. Default: gt''')

parser.add_argument('--limit', '-l', type=int, default=4641652,
                            help='''Boundar/limit for the relation ship with POS field.
                                    e.g. If operator is gt and this is 4641652 and add=-4641652,
                                        then positions > 4641652 will have 4641652 subtracted.''')

args = parser.parse_args()




# EXECUTE
sam = Sam(args.sam)
    
for read in sam:
    if read.pos_field_is(args.operator, args.limit):
        if args.add:
            read.update_pos_field(add=args.add)
        elif args.replace:
            read.update_pos_field(replace=args.replace)
    print read




