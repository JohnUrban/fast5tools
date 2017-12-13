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
    return maxbar, maxprob, maxscore, meanscore, scorestring, probs_string, alignments[maxbar].q_pos, alignments[maxbar].q_end, alignments[maxbar].r_pos, alignments[maxbar].r_end, alignments[maxbar].cigar_str, epowersum #, alignments[maxbar].orig_query,  alignments[maxbar].orig_ref 


def choose_barcode(barcodes, sw, read_seq, search_start, search_end, ref_name='', use_entire_barcode=True, return_alns=False):

    search_seq = read_seq[search_start:search_end]
##    print len(search_seq)
    # get alignments
    alignments = get_alignments(barcodes, sw, search_seq, ref_name, use_entire_barcode)
    #get best barcode info
    if return_alns:
        return get_max_scoring_barcode(alignments), alignments
    else:
        return get_max_scoring_barcode(alignments)

def get_formatted_pairwise_alignments(barcodechoice, getall=False):
    ''' barcode choice = output of choose_barcode() with  return_alns=True'''
    outstring = ''
    if getall:
        for barcode, alignment in barcodechoice[1].iteritems():
            outstring += get_formatted_pairwise_alignment(alignment,report_prob=(np.e**alignment.score)/barcodechoice[0][-1])
        return outstring
    else:
        best = barcodechoice[1][barcodechoice[0][0]]
        return get_formatted_pairwise_alignment(best, report_prob=barcodechoice[0][1])
            
    
def get_formatted_pairwise_alignment(alignment, blocksize=100, e_s=0.051, e_i=0.049, e_d=0.078, report_prob=False, report_evalue=False, with_unaligned_portion=False):
    ## TODO: printout unaligned portion as part of alignment viz
    ## For report_prob, you need to give the prob
    ## e_i, e_d, and e_s are insertion/deletion/substitution errors found in early MinION sequencing by Jain et al: Improved data analysis for the MinION nanopore sequencer
    ## Can compute prob of alignment by p_m=1-e_i-e_d-e_s; 
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
            n_d += e[0] #1
        elif e[1] == 'I':
            ref += '-'*e[0]
            query += alignment.orig_query[q_i:q_i+e[0]]
            sticks += ' '*e[0]
            q_i += e[0]
            n_i += e[0] #1
####    print alignment.matches, alignment.mismatches, alignment.mismatches-n_mm-n_d-n_i, len(query), sum([n_m, n_mm, n_d, n_i]),  sum([n_mm, n_d, n_i])
##    print "Ref (top):", alignment.r_name
##    print "Query (bottom):", alignment.q_name
##    print "Match:"+str(n_m), "Mismatch:"+str(n_mm), "Deletion:"+str(n_d), "Insertion:"+str(n_i), "PercentIdentity:"+str(100.0*n_m/sum([n_m, n_mm, n_d, n_i]))
##    for i in range(0, len(ref), blocksize):
##        print ref[i:i+blocksize]
##        print sticks[i:i+blocksize]
##        print query[i:i+blocksize]
##        print
    qbases = n_m + n_mm + n_i
    rbases = n_m + n_mm + n_d
    n_u = len(alignment.orig_query) - qbases ## these were bases not in the alignment. Since barcode is query... i.e. pieces of barcode not found in read
    assert qbases == alignment.q_end-alignment.q_pos and rbases == alignment.r_end-alignment.r_pos
    assert len(ref) == len(query)
    assert len(ref) == sum([n_m, n_mm, n_d, n_i])
    outstring = alignment.q_name + '\n'
    outstring += "Ref (top): " + alignment.r_name + ' ' + str(alignment.r_pos) + '-' + str(alignment.r_end) + ' r_bases_aligned:' + str(rbases) + ' pct_r_bases_aligned:' + str(100.0*rbases/len(alignment.orig_ref)) + ' refLen:' + str(len(alignment.orig_ref)) + ' bp\n'
    outstring += "Query (bottom): " + alignment.q_name + ' ' + str(alignment.q_pos) + '-' + str(alignment.q_end) + ' q_bases_aligned:' + str(qbases) + ' pct_q_bases_aligned:'+str(100.0*qbases/len(alignment.orig_query)) + ' queryLen:' + str(len(alignment.orig_query)) + ' bp\n'
    stats = ["AS:" + str(alignment.score), "Match:"+str(n_m), "Mismatch:"+str(n_mm), "Deletion:"+str(n_d), "Insertion:"+str(n_i), "AlignmentLength:"+str(len(ref)), "PercentIdentity:"+str(100.0*n_m/sum([n_m, n_mm, n_d, n_i])),  "Barcode_Unaligned:"+str(n_u), "PercentIdentity_with_unaligned:"+str(100.0*n_m/sum([n_m, n_mm, n_d, n_i, n_u]))]
    outstring += (' ').join(stats) + '\n'
    p = np.e**(-1*alignment.score)
    n = len(alignment.orig_ref)*len(alignment.orig_query)
    evalue = n*p
    stats = ['n:' + str(n), 'p:' + str(p), 'e_value:' + str(evalue)]
    outstring += (' ').join(stats) + '\n'
    p_m = 1 - e_i - e_d - e_s
    p_minion_aln =  (p_m**n_m) * (e_s**n_mm) * (e_d**n_d) * (e_i**n_i)
    Esum = e_s + e_i + e_d
    E_s = e_s/Esum
    E_i = e_i/Esum
    E_d = e_d/Esum
    p_minion_un = (e_s**(E_s*n_u)) * (e_i**(E_i*n_u)) *(e_d**(E_d*n_u)) ## since it is unclear if the unaligned were subs/dels/ins I am making use of all
    p_minion = p_minion_aln * p_minion_un
    N=sum([n_m, n_mm, n_i, n_d])
    N2 = N - n_m
    N3 = N2 - n_mm
    N4 = N3 - n_i
    N5 = N4 - n_d
    NcK = nchoosek(N, n_m) * nchoosek(N2, n_mm) * nchoosek(N3, n_i) * nchoosek(N4, n_d) 
    norm_p_minion_aln = NcK * p_minion_aln
    N_u = N + n_u
    norm_p_minion = norm_p_minion_aln * p_minion_un ## There seems to be no reason to multiply the unaligned by a nchoosek b/c they are fixed at the ends and correspond to some "single" unknown composite error
    stats = ['p_minion_aln:' + str(p_minion_aln), 'p_minion_un:' + str(p_minion_un), 'p_minion:' + str(p_minion), 'norm_p_minion_aln:'+str(norm_p_minion_aln), 'norm_p_minion:'+str(norm_p_minion)] ## p_minion not necessarily comparable when barcodes are different lengths - can divide by bc_len or q*r maybe
    outstring += (' ').join(stats) + '\n'
    
    ## actually since the read is the "ref" here.... n_d and n_i are dels/ins from/in read. Barcode is the real "reference" meaning for minion probs we need e_i**n_d and e_d**n_i 
    p_minion_aln =  (p_m**n_m) * (e_s**n_mm) * (e_d**n_i) * (e_i**n_d)
    p_minion = p_minion_aln * p_minion_un
    norm_p_minion_aln = NcK * p_minion_aln
    norm_p_minion = norm_p_minion_aln * p_minion_un
    stats = ['p_minion_aln:' + str(p_minion_aln), 'p_minion_un:' + str(p_minion_un), 'p_minion:' + str(p_minion), 'norm_p_minion_aln:'+str(norm_p_minion_aln), 'norm_p_minion:'+str(norm_p_minion)] ## p_minion not necessarily comparable when barcodes are different lengths - can divide by bc_len or q*r maybe
    outstring += (' ').join(stats) + '\n'
    
    binom_prob = nchoosek(sum([n_m, n_mm, n_d, n_i]), n_m) * (p_m**n_m) * ((1-p_m)**(n_mm+n_i+n_d))
    outstring += 'binom_prob_k_matches_in_alignment:' + str(binom_prob) + '\n'

    binom_prob = nchoosek(sum([n_m, n_mm, n_d, n_i, n_u]), n_m) * (p_m**n_m) * ((1-p_m)**(n_mm+n_i+n_d+n_u))
    outstring += 'binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode:' + str(binom_prob) + '\n'
 
    #in barcode only
    binom_prob = nchoosek(len(alignment.orig_query), n_m) * (p_m**n_m) * ((1-p_m)**(len(alignment.orig_query)-n_m))
    outstring += 'binom_prob_k_matches_in_barcode:' + str(binom_prob) + '\n'

    if report_prob is not False:
        outstring += "Marginalized Probability Given Barcode Set: " + str(report_prob) + "\n"
    for i in range(0, len(ref), blocksize):
        outstring += ref[i:i+blocksize] + '\n'
        outstring += sticks[i:i+blocksize] + '\n'
        outstring += query[i:i+blocksize] + '\n\n'
    return outstring
    
