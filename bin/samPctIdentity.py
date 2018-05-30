#!/usr/bin/env python2.7
import argparse
from collections import defaultdict
from fast5tools.samclass import *
from fast5tools.samops import *



parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Input = SAM file sorted by Read Name using "samtools sort -n".
    -- FOR NOW, MUST BE SAM -- NOT BAM -- but can be STDIN SAM
    -- Assumes reads can have split alignments
    
    Output = Tab-delimited
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group(required=True)
parser_input.add_argument('--sam', '-s',
                   type= str, default=False,
                   help='''Input file in SAM format.''')

parser_input.add_argument('--prefix', '-p',
                   type=str, default=False,
                   help='''To name some out files.... e.g. readsPrefix''')

## FOR NOW, MUST BE SAM -- NOT BAM -- but can be STDIN SAM
##parser_input.add_argument('--bam', '-b',
##                   type= str, default=False,
##                   help='''Input file in BAM format.''')

args = parser.parse_args()



sam = SamSplitAlnAggregator(args.sam)



def perfectAlignment(read, pctidobj, numRec):
    #if read.get_read_length() == pctidobj['match']:
    #    return 1
    ## If entire alignment length is just matches and read is not split into more than 1 aln
    if (pctidobj['MDIU'] == pctidobj['match']) and numRec == 1:
        return 1
    return 0

header = ['read', 'numSplitRecords', 'perfect', 'match', 'mismatch', 'deletion', 'insertion', 'unaligned', 'pct_id_1', 'pct_id_2', 'pct_id_3', 'pct_id_4', 'pct_id_5']
print '\t'.join(header)
total = 0
totalAligned = 0
totalPerfect = 0
for read in sam:
    total += 1
    if read.has_alignments():
        totalAligned += 1
        pctid =  read.get_pct_identity()
        #out = [read.get_qname_field(), read.get_rname_field(), read.get_pos_field()]
        numRec = read.get_num_aln()
        perfect = perfectAlignment(read, pctid, numRec)
        totalPerfect += perfect
        out = [read.get_read_name(), numRec, perfect]
        out += [pctid['match'], pctid['mismatch'], pctid['del'], pctid['ins'], pctid['unaligned']]
        out += pctid['pctid']
        print '\t'.join([str(e) for e in out])
##        print pctid, abs(pctid[1]-pctid[0]), abs(pctid[1]-pctid[3]), abs(pctid[3]-pctid[4]), abs(pctid[1]-pctid[0]) > abs(pctid[1]-pctid[3]),  abs(pctid[1]-pctid[0]) > abs(pctid[3]-pctid[4])
##        print read.get_per_base_align_status_for_read(), read.get_number_bases_in_read_aligned(), read.get_number_bases_in_read_not_aligned(), read.get_pct_of_read_aligned(), read.get_pct_of_read_not_aligned(), pctid[0], pctid[1], pctid[2], read.get_pct_identity_proxy()


totalImperfect = totalAligned-totalPerfect
perfectMsg = '''Total Reads = %d
Total Aligned Reads = %d
Total Perfect Alignments = %d
Total Imperfect Alignments = %d
Proportion Perfect Alignments = %f
Proportion Imperfect Alignments = %f
''' % (total, totalAligned, totalPerfect, totalImperfect, float(totalPerfect)/totalAligned, float(totalImperfect)/totalAligned)


if not args.prefix:
    sys.stderr.write(perfectMsg)
else:
    with open(args.prefix + '.summary.txt') as f:
        f.write(perfectMsg)

