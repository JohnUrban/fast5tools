from Bio import SeqIO
from fast5tools.samclass import *

def alignment_summary(sam):
    '''sam is Sam object'''
    for read in sam:
        print ("\t").join([str(e) for e in [read.get_qname_field(), read.get_rname_field(), read.get_pos_field(), read.get_read_len(), read.get_AS_field(), read.get_SEQ_len(), read.get_SEQ_len_without_clipped_regions(), read.get_reference_aln_len(), read.get_edit_dist_field(), read.get_edit_dist_with_clipping(), read.get_clipping_dist(), read.get_fast5_field()]])

def get_genomic_windows(samfilepath, flank=0.25, merge_dist=0, majority=0.5, require_order=False, require_strand=False, reference=False):
    ''' samfilepath is path to sam file -- used to create "sam" SamSplitAlnAggregator object
        output is BED-like, but is 1-based start and closed end..

        flank =
        merge_dist =
        majority =
        require_order =
        require_strand =
        reference = path to fasta file used to extra referenc genomic window sequences from (optional)'''

    sam = SamSplitAlnAggregator(samfilepath)
    if reference:
        #Load into dictionary
        refseq = {}
        for fa in SeqIO.parse(reference, 'fasta'):
            refseq[str(fa.name)] = str(fa.seq)
    for read in sam:
##        read = SplitReadSamRecord(read)
        if read.has_alignments():
            ans = read.get_genomic_window(flank=flank, merge_dist=merge_dist, majority=majority, require_order=require_order, require_strand=require_strand)
            if ans is not None:
                seq = []
                if reference:
                    chrom = ans[0]
                    start = ans[1]-1
                    if start < 0:
                        start = 0
                    end = ans[2]
                    seq = [refseq[chrom][start:end]]
                print ("\t").join( [read.get_read_name()] + [str(e) for e in ans] + [read.get_fast5_info()] + seq)
                ## if reference fasta provided, then print out sequence in last field of table
                ## should any other SAM info be included?
                ## include edit dist? num matches? %id? of longest only? all alns? all_alns vs readlen?
            else:
                ## shouldn't come here -- can clean out eventually
                print "NONETYPE", ans
        else:
            seq = []
            if reference:
                seq = ['*']
            print ("\t").join( [read.get_read_name()] + ['*']*3 + ['no_alignments'] + [read.get_fast5_info()] + seq )



def get_aln_pct_identity_v_Q(sam):
    for read in sam:
        cigar = read.cigar
##        print read.get_pct_identity()
        print ("\t").join( [str(e) for e in [read.get_edit_dist_field(), read.get_read_len(), read.get_SEQ_len(), read.get_SEQ_len_without_clipped_regions(), read.get_fast5_field(), read.get_cigar_counts()]])
        
def get_pct_identity_v_Q(sam):
    ''' sam is SamSplitAlnAggregator object'''
    for read in sam:
        print read.get_edit_dist_fields(), sum(read.get_edit_dist_fields()), read.get_read_length()
        print ("\t").join( [read.get_read_name(), read.get_fast5_info()])



