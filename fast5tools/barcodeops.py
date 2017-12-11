import swalign, re, nwalign
from Bio import SeqIO
import numpy as np
from scipy.special import comb as nchoosek


def read_in_barcodes(barcodetable=False, barcodefasta=False, barcolumns='3,5'):
    ''' barcode table is path to tab-sep file with barcode names and seqs. Assumes columns for those are 3 and 5 respectively.
        barcolumns - default 3,5. Otherwise provide comma-sep pair.
        barcodefasta = path to fasta with barcode names and seqs
        Can only pick table to fasta.'''
    assert (barcodetable and not barcodefasta) or (barcodefasta and not barcodetable)
    barcodes = {}
    if barcodetable:
        namecol, seqcol = (int(e)-1 for e in barcolums.split(','))
        with open(barcodetable, 'r') as bfile:
            for row in bfile:
                if row[0] != '#':
                    row = row.strip().split()
                    barcodes[ row[namecol] ] = row[seqcol]
    elif barcodefasta:
        for fa in SeqIO.parse(barcodefasta, 'fasta'):
            barcodes[str(fa.name)] = str(fa.seq)
    return barcodes


def create_sw_scoring_matrix(match=2, mismatch=-2):
    return swalign.NucleotideScoringMatrix(2, -2)

def get_alignment_object(match, mismatch, gap_open, gap_ext, gap_decay, verbose=False, globalalign=False, full_query=False):
    ## See https://github.com/mbreese/swalign/blob/master/bin/swalign for help/docs
    ''' scoring is a scoring matrix object'''
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(
        scoring,
        gap_open, gap_ext,
        gap_decay, verbose=verbose, globalalign=globalalign, full_query=full_query)
    return sw


def get_barcode_search_length(barcodes, barlength=False, barflank=False, minbarlen=80):
    if barlength is False:
        barlength = 0
        for name, seq in barcodes.iteritems():
            seqlen = len(seq)
            if seqlen > barlength:
                barlength = seqlen
    if barflank is False: #this way ensures 0 is not interpreted as False
        barflank = barlength
    search_len = max(minbarlen, barlength+barflank)
    return search_len
    

def get_smallest_barcode_length(barcodes):
    barlength = float('inf')
    for name, seq in barcodes.iteritems():
        seqlen = len(seq)
        if seqlen < barlength:
            barlength = seqlen
    return barlength


def get_alignments(barcodes, sw, sequence, ref_name='', use_entire_barcode=True):
    ''' barcodes is dict of barcodes
        sequence is given -- usually first X bp of a read.
        use_entire_barcode -- if False, it will trim barcodes to length of shortest.'''
    alignments = {}
    if not use_entire_barcode:
        minbarlen = get_smallest_barcode_length(barcodes)
    for barcode_name, barcode_seq in barcodes.iteritems():
##        print barcode_name
        if not use_entire_barcode:
            barcode_seq = barcode_seq[:minbarlen]
        alignment = sw.align(query=barcode_seq, ref=sequence, query_name=barcode_name, ref_name=ref_name)
        alignments[barcode_name] = alignment
    return alignments


def get_max_scoring_barcode(alignments, match_prob=0.75):
    mismatch_prob = 1 - match_prob
    maxscore = float('-inf')
    scores = []
    for barcode, alignment in alignments.iteritems():
        scores.append( float(alignment.score) )
        ## Below - was just exploring other metrics...
##        D = sum([int(e) for e in re.findall('(\d)D',alignment.cigar_str)])
##        I = sum([int(e) for e in re.findall('(\d)I',alignment.cigar_str)])
##        DI = D+I
##        p = np.e**(-1*alignment.score)
##        n = len(alignment.orig_ref)*len(alignment.orig_query)
##        evalue = n*p
##        N = len(alignment.orig_query) #min(len(alignment.orig_ref), len(alignment.orig_query))
##        mm = N-alignment.matches
##        lap =  (match_prob**alignment.matches) * (mismatch_prob**mm)
##        alt_lap =  (match_prob**alignment.matches) * ((mismatch_prob*0.66)**(mm*0.5)) * ((mismatch_prob*0.34)**(mm*0.5))
##        binom = nchoosek(N, alignment.matches) * lap
##        alt_bin = nchoosek(N, alignment.matches) * alt_lap
##        print barcode, alignment.score, alignment.q_pos, alignment.q_end, alignment.r_pos, alignment.r_end, alignment.matches, alignment.mismatches, n, p,  evalue, lap, binom, N, alignment.matches, mm, DI,  alt_lap, alt_bin
##        get_formatted_pairwise_alignment(alignment) ##delete
        if alignment.score > maxscore:
            maxscore = alignment.score
            maxbar = barcode
    scores = np.array(scores)
    scorestring = (',').join([str(e) for e in scores])
    meanscore = scores.mean()
    epowers = [np.e**score for score in scores]
    epowersum = float(sum(epowers))
    probs_string = (',').join([str(epower/epowersum) for epower in epowers])
    maxprob = np.e**maxscore / epowersum
    # could also give cigar: alignments[maxbar].cigar_str
##    D = sum([int(e) for e in re.findall('(\d)D',alignments[maxbar].cigar_str)])
##    I = sum([int(e) for e in re.findall('(\d)I',alignments[maxbar].cigar_str)])
##    M = sum([int(e) for e in re.findall('(\d)M',alignments[maxbar].cigar_str)])
##    print "FINAL:"
##    get_formatted_pairwise_alignment(alignments[maxbar])
    ##  alignments[maxbar].identity, alignments[maxbar].matches, alignments[maxbar].mismatches, I, D, M,
    return maxbar, maxprob, maxscore, meanscore, scorestring, probs_string, alignments[maxbar].q_pos, alignments[maxbar].q_end, alignments[maxbar].r_pos, alignments[maxbar].r_end, alignments[maxbar].cigar_str #, alignments[maxbar].orig_query,  alignments[maxbar].orig_ref 


def choose_barcode(barcodes, sw, read_seq, search_start, search_end, ref_name='', use_entire_barcode=True):

    search_seq = read_seq[search_start:search_end]
##    print len(search_seq)
    # get alignments
    alignments = get_alignments(barcodes, sw, search_seq, ref_name, use_entire_barcode)
    #get best barcode info
    return get_max_scoring_barcode(alignments)


def get_formatted_pairwise_alignment(alignment, blocksize=100):
    n_m = 0
    n_mm = 0
    n_d = 0
    n_i = 0
    ref = ''
    query = ''
    sticks = ''
##    refseq = alignment.orig_ref[alignment.r_pos:alignment.r_end]
    r_i = alignment.r_pos
    q_i = alignment.q_pos
    for e in alignment.cigar:
        if e[1] == 'M':
            newref = alignment.orig_ref[r_i:r_i+e[0]]
            newquer = alignment.orig_query[q_i:q_i+e[0]]
            ref += newref
            query += newquer
            r_i += e[0]
            q_i += e[0]
            for i in range(len(newref)):
                sticks += '|' if newref[i].upper() == newquer[i].upper() else ' '
                n_m += 1 if newref[i].upper() == newquer[i].upper() else 0
                n_mm += 0 if newref[i].upper() == newquer[i].upper() else 1
        if e[1] == 'D':
            ref += alignment.orig_ref[r_i:r_i+e[0]]
            query += '-'*e[0]
            sticks += ' '*e[0]
            r_i += e[0]
            n_d += 1
        elif e[1] == 'I':
            ref += '-'*e[0]
            query += alignment.orig_query[q_i:q_i+e[0]]
            sticks += ' '*e[0]
            q_i += e[0]
            n_i += 1
    assert len(ref) == len(query)
    print "Match:"+str(n_m), "Mismatch:"+str(n_mm), "Deletion:"+str(n_d), "Insertion:"+str(n_i), "PercentIdentity:"+str(100.0*n_m/sum([n_m, n_mm, n_d, n_i]))
##    print alignment.matches, alignment.mismatches, alignment.mismatches-n_mm-n_d-n_i, len(query), sum([n_m, n_mm, n_d, n_i]),  sum([n_mm, n_d, n_i])
    print "Ref (top):", alignment.r_name
    print "Query (bottom):", alignment.q_name
    for i in range(0, len(ref), blocksize):
        print ref[i:i+blocksize]
        print sticks[i:i+blocksize]
        print query[i:i+blocksize]
        print
    
