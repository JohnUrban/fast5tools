#!/usr/bin/env python2.7
import argparse, sys
from collections import Counter, defaultdict
from string import maketrans
import re

from fast5tools.samclass import *


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Input = samtools mpileup output.
    Given sorted BAM/SAM do:
        samtools mpileup -f ref.fa reads.bam | samRefPosAlnComposition.py -
        samtools mpileup -f ref.fa reads.bam | samRefPosAlnComposition.py stdin
    
    Output = tab-delimited file describing all positions or each non-perfect position of reference (where a perfect position is 100% agreement w/ ref on both strands).
        Reports the alignment coordinates only (clipping ignored) and number of mismatches therein.
        Columns:
            chr, start, end, number mismatches, read name

    
    John Urban (2015, 2016, 2017, 2018)
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-m', '--mpileup', 
                   type= str, default='stdin', required=False,
                   help='''Paths to SAMtools mpileup output file. Can be "-" or "stdin" for standard input.''')
parser.add_argument('-S', '--stranded', default=False, action='store_true',
                    help='''Add stranded information. This will also give counts, probs, and composition specific to each strand.
                        If the library is not strand-specific, this might be somewhat meaningless.
                        Use and interpret with caution.''')

parser.add_argument('-d', '--mindepth', default=0, type=float,
                    help='''Look at only positions with at least this much depth. Default: 0.''')

parser.add_argument('-t', '--threshold', default=1.1, type=float,
                    help='''A probability threshold to partition positions in genome into pass/fail.
                            Use this with --stringency and optionally --above.
                            Default: Everything is returned by using 1.1 as the value, which triggers it to return "False" for >= threshold.
                            When used without --above, it should pass/keep all positions below threshold.
                            When used with --above, it should pass/keep all above or equal to threshold. (None by default)''')
parser.add_argument('-a', '--above', default=False, action='store_true',
                    help='''Used with --threshold (and --stringency).
                            Not specifying --above should pass/keep all positions below threshold.
                            Specfying --above should pass/keep all above or equal to threshold.''')
parser.add_argument('-r', '--reject_undefined', default=False, action='store_true',
                    help='''Default: False. Undefined values are reported.
                            IMPORTANT: you should reject the undefined values in some analyses -- can do this "ex post facto" w/ awk (for/against "NA") if need be.
                            Used with --threshold (and --stringency).
                            Not specifying --reject_undefined should pass/keep all positions where the probability was "NA".
                            Specfying --reject_undefined should fail/reject them.
                            An example where the probability is undefined in the case of using stringency=1 when there are no matches or mismatches, only deletions.
                            In such a case, marginalizing out deletions and Ns would result in 0/(0+0).
                            NOTE: This only rejects/keeps NA values when encountered as part of the stringency level.
                            Thus if you are using stringency=1 and get an 'NA' as described above, it will act on it.
                            However, if stringency=2 or stringency=3 was used, that 'NA' would not be seen as part of the threshold test.
                            Thus, some columns can still have 'NA' values even if you use this option.
                            Again- since this can also be handled ex post facto with AWK, it might be better to let all NA values through and worry later.
                            ''')

parser.add_argument('-s', '--stringency', default=1, type=int,
                    help='''If thresholding is applied, this changes the stringency of filtering.
                        Options are (1, 2, 3).
                        3 means p(match) = 1. There can only be matches over the position - no mismatch, no del, no N.
                        2 means marg_out_N_p(match) = 1. No mismatches, no dels (Ns are allowed/ignored).
                        1 means marg_out_N_D_p(match) = 1. No mismatches (Ns and dels are allowed/ignored).
                        Default stringency level = 1.

                        If you have 4 positions with 10 reads.
                        Pos 1 has 7 matches, 1 mismatch, 1 del, and 1 N.
                        Pos 2 has 8 matches, 1 mismatch, 1 del.
                        Pos 3 has 9 matches, 1 mismatch.
                        Pos 4 has 10 matches.
                        Pos    p(match)    marg_out_N_p(match)    marg_out_N_D_p(match)
                        Pos 1: 0.7    0.778    0.875
                        Pos 2: 0.8    0.8    0.889
                        Pos 3: 0.9    0.9    0.9
                        Pos 4: 1.0    1.0    1.0

                        Using any stringency with threshold = 1 -- would partition Pos 4 away from Pos1,2,3
                        Using threshold=0.85:
                           Stringency=1 would keep/reject all
                           Stringency=2 would partition Pos1,2 from Pos3,4
                           Stringency=3 would do same as 2 here...''')
parser.add_argument('-no', '--no-header', dest='no_header', default=False, action='store_true', help='''No header in output.''')
args = parser.parse_args()


###################################################
'''TODO'''
###################################################
'''
    Work base and map qualities in as an option.
        This is a low priority since there are many great programs (incl SAMtools) to look at SNPs/etc that already do this.
'''

###################################################
'''notes to self'''
###################################################
'''
The mpileup documentation does not make clear how deletions are represented.
It says it will be something like -1G, and the following lines will be represented as "*"...

So does the deletion start on the line with -1G or is it only in the following lines?
- Mismatches occur on the position they occur..
- As far as my toy analyses here show, only the following lines with '*' represent the deletions...
        - If I delete the 14th base, -1X shows up on the 13th position and '*' at the 14th
	- Similarly if I delete 5 bases 24-28 inclusive, '-5XXXXX' shows up on line 23 and '*' shows up on lines 24-28 inclusive
- Insertions show up on the line preceding the position of the insertion in the read -- i.e. sequence is inserted after this base.
	- If insert something after pos 14, line 14 says +1X (no '*' stuff)
	- If insert 5 bases after line 24, line 24 says +5XXXXX (no '*' stuff on following lines)
- So the deletions mimic the insertion stuff in a way saying an insertion/deletion follow this line/position.
- But for my purposes there is no sense in keeping track of the '-1X' parts, only the '*'
- Also: 
	- reverse complement reads show differences on the fwd strand in lower-case -- i.e. if you see ..a.. that means there was a T in the revcomp
            - But importantly - it means that it is reflecting composition on the fwd strand (it is not 't', it is 'a')
	- ^[ precedes a read's match/mm/indel status
	- $ follows a read's match/mm/indel status
	- the read match/mm/indel statuses are in order found in BAM (but this likely gets much trickier to track as coverage gets more complicated)
	- the base qualities in subsequent column only seem to track with matches and mismatches -- of course the read has no qual for a deletion and I guess insertion qualities are ignored in this format.
		- there are no qualities for the start and end symbols in the qual column
		- however, the symbol following start symbol is the mapq score for the read

Made reads with known mutations -- all reads are from the same sequence starting at the beginning of the lambda genome sequence.
>wt-1
>wt-2
>wt-3
>wt-4
>wt-5
>4A-1
>4A-2
>14D-1
>24DDDDD-1
>14I-one-insertion-after14th-base
>24IIIII-five-insertions-after-24th-base
>wt-rc-1
>4A-rc-1
>14D-1-rc
>14I-rc-one-insertion-after14th-base

## With seqs
>wt-1
GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>wt-2
GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>wt-3
GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>wt-4
GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>wt-5
GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>4A-1
GGGAGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>4A-2
GGGAGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>14D-1
GGGCGGCGACCTCCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>24DDDDD-1
GGGCGGCGACCTCGCGGGTTTTCTTATGAAAATTTTCCGGTTT
>14I-one-insertion-after14th-base
GGGCGGCGACCTCGACGGGTTTTCGCTATTTATGAAAATTTTCCGGTTT
>24IIIII-five-insertions-after-24th-base
GGGCGGCGACCTCGCGGGTTTTCGAAAAACTATTTATGAAAATTTTCCGGTTT
>wt-rc-1
AAACCGGAAAATTTTCATAAATAGCGAAAACCCGCGAGGTCGCCGCCC
>4A-rc-1
AAACCGGAAAATTTTCATAAATAGCGAAAACCCGCGAGGTCGCCTCCC
>14D-1-rc
AAACCGGAAAATTTTCATAAATAGCGAAAACCCGGAGGTCGCCGCCC
>14I-rc-one-insertion-after14th-base
AAACCGGAAAATTTTCATAAATAGCGAAAACCCGTCGAGGTCGCCGCCC

#Map reads
bwa mem idx/lambdaGenomeSequence manual-reads.fa | samtools sort > manual-reads.bam

# Looking at BAM
samtools view *bam | cut -f 1,2,3,4,5,6,12,13,14,15,16

#read    flag    ref     pos     mapq    cigar   editdist        MD      AS      XS
wt-1	0	lambdaGenome	1	60	48M	NM:i:0	MD:Z:48	AS:i:48	XS:i:0
wt-2	0	lambdaGenome	1	60	48M	NM:i:0	MD:Z:48	AS:i:48	XS:i:0
wt-3	0	lambdaGenome	1	60	48M	NM:i:0	MD:Z:48	AS:i:48	XS:i:0
wt-4	0	lambdaGenome	1	60	48M	NM:i:0	MD:Z:48	AS:i:48	XS:i:0
wt-5	0	lambdaGenome	1	60	48M	NM:i:0	MD:Z:48	AS:i:48	XS:i:0
4A-1	0	lambdaGenome	1	60	48M	NM:i:1	MD:Z:3C44	AS:i:44	XS:i:0
4A-2	0	lambdaGenome	1	60	48M	NM:i:1	MD:Z:3C44	AS:i:44	XS:i:0
14D-1	0	lambdaGenome	1	60	13M1D34M	NM:i:1	MD:Z:13^G34	AS:i:40	XS:i:0
24DDDDD-1	0	lambdaGenome	1	59	23M5D20M	NM:i:5	MD:Z:23^GCTAT20	AS:i:32	XS:i:0
14I-one-insertion-after14th-base	0	lambdaGenome	1	60	14M1I34M	NM:i:1	MD:Z:48	AS:i:41	XS:i:0
24IIIII-five-insertions-after-24th-base	0	lambdaGenome	1	48	24M5I24M	NM:i:5	MD:Z:48	AS:i:37	XS:i:0
wt-rc-1	16	lambdaGenome	1	60	48M	NM:i:0	MD:Z:48	AS:i:48	XS:i:0
4A-rc-1	16	lambdaGenome	1	60	48M	NM:i:1	MD:Z:3C44	AS:i:44	XS:i:0
14D-1-rc	16	lambdaGenome	1	60	13M1D34M	NM:i:1	MD:Z:13^G34	AS:i:40	XS:i:0
14I-rc-one-insertion-after14th-base	16	lambdaGenome	1	60	14M1I34M	NM:i:1	MD:Z:48	AS:i:41	XS:i:0

#Look at mpileup of alignments
samtools mpileup -f idx/lambdaGenomeSequence.fa manual-reads.bam | less


lambdaGenome	1	G	15	^].^].^].^].^].^].^].^].^\.^].^Q.^],^],^],^],	~~~~~~~~~~~~~~~
lambdaGenome	2	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	3	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	4	C	15	.....AA....,a,,	~~~~~~~~~~~~~~~
lambdaGenome	5	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	6	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	7	C	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	8	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	9	A	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	10	C	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	11	C	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	12	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	13	C	15	........-1G...,,,-1g,	~~~~~~~~~~~~~~~
lambdaGenome	14	G	15	.......*..+1A.,,*,+1a	~~~~~~~~~~~~~~~
lambdaGenome	15	C	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	16	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	17	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	18	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	19	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	20	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	21	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	22	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	23	C	15	.........-5GCTAT..,,,,	~~~~~~~~~~~~~~~
lambdaGenome	24	G	15	........*..+5AAAAA,,,,	~~~~~~~~~~~~~~~
lambdaGenome	25	C	15	........*..,,,,	~~~~~~~~~~~~~~~
lambdaGenome	26	T	15	........*..,,,,	~~~~~~~~~~~~~~~
lambdaGenome	27	A	15	........*..,,,,	~~~~~~~~~~~~~~~
lambdaGenome	28	T	15	........*..,,,,	~~~~~~~~~~~~~~~
lambdaGenome	29	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	30	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	31	A	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	32	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	33	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	34	A	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	35	A	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	36	A	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	37	A	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	38	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	39	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	40	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	41	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	42	C	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	43	C	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	44	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	45	G	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	46	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	47	T	15	...........,,,,	~~~~~~~~~~~~~~~
lambdaGenome	48	T	15	.$.$.$.$.$.$.$.$.$.$.$,$,$,$,$	~~~~~~~~~~~~~~~

'''


###################################################
'''execute'''
###################################################
mpileup = Mpileup(args.mpileup, args.stranded)

if not args.no_header:
    print mpileup.header()

for rec in mpileup:
    if rec.has_min_depth(args.mindepth):
        if rec.passes_threshold(threshold=args.threshold, stringency=args.stringency, strand=0, above=args.above, reject_undefined=args.reject_undefined): ## Want option of only printing imperfect positions
            print rec
mpileup.close()




