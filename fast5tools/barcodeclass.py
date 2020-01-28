import swalign, re, nwalign
from Bio import SeqIO
import numpy as np
from scipy.special import comb as nchoosek



class BarcodeChoice(object):
    def __init__(self, barcodes, sw, read_seq, search_start, search_end, ref_name='', compute_all=False, e_s=0.051, e_i=0.049, e_d=0.078, use_entire_barcode=True):
        ## Goal: be able to access all alignments, all scores, all probs, all alt probs -- etc as dicts
        ##  This way it would be trivial to output all scores wanted... as in aln viz
        ''' barcodes = output of read_in_barcodes() from barcodeops.py
            sw = output of get_alignment_object() from barcodeops.py
            search_len = output of get_barcode_search_length() from barcodeops.py'''
        self.barcodes = barcodes
        self.sw = sw
        self.read_seq = read_seq
        self.search_start = search_start
        self.search_end = search_end
        self.search_seq = read_seq[search_start:search_end]
        self.ref_name = ref_name
        self.use_entire_barcode = use_entire_barcode
        self.alignments = None
        self.scores = None
        self.score_powers = None
        self._get_alignments()
        self._get_max_scoring_barcode()
        self._get_score_probabilities()
        self.counts = {}
        self.alnstrings = {}
        self.expectedvalues = {}
        self.minionprobs = {}
        self.binomprobs = {}
        self.margin_minion_probs = {}
        self.margin_binom_probs = {}
        self.minion_prob_sums = {}
        self.binomial_prob_sums = {}
        if compute_all:
            self._compute_all_formatted_alignment_strings()
            self._get_expected_values()
            self._get_marginalized_minion_probabilities() ## default error rates if not changed
            self._get_marginalized_binomial_probabilities() #default error rates if not changed

    ## Public
    def get_barcode_start(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].q_pos

    def get_barcode_end(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].q_end

    def get_read_window_start(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].r_pos

    def get_read_window_end(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].r_end

    def get_barcode_aln_len(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].q_end - self.alignments[barcode].q_pos

    def get_read_aln_len(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].r_end - self.alignments[barcode].r_pos

    def get_score(self, barcode):
        self._alignment_check()
        self._score_check()
        return self.alignments[barcode].score
    
    def get_score_probability(self, barcode):
        self._alignment_check()
        self._score_check()
        return self.score_probabilities[barcode]

    def get_cigar(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].cigar

    def get_cigar_string(self, barcode):
        self._alignment_check()
        return self.alignments[barcode].cigar_str

    def get_max_scoring_barcode(self):
        self._alignment_check()
        self._score_check()
        return self.maxscoringbarcode

    def get_max_score(self):
        self._alignment_check()
        self._score_check()
        return self.maxscore

    def get_max_score_probability(self):
        self._alignment_check()
        self._score_check()
        self._score_probability_check()
        return self.score_probabilities[self.maxscoringbarcode]
    
    def get_mean_score(self):
        self._alignment_check()
        self._score_check()
        return self.scores.mean()


    def get_expected_value(self, barcode, get='evalue'):
        ''' get = evalue, n, or p'''
        try:
            return self.expectedvalues[barcode][get]
        except:
            self._get_expected_value(barcode)
            return self.expectedvalues[barcode][get]

    def get_minion_probability(self, barcode, get='norm_p_minion', e_s=0.051, e_i=0.049, e_d=0.078):
        '''get in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']'''
        try:
            if [e_s, e_i, e_d] != self.minionprobs[barcode]['params']:
                self._get_minion_probability(barcode, e_s, e_i, e_d)
            return self.minionprobs[barcode][get]
        except: 
            self._get_minion_probability(barcode, e_s, e_i, e_d)
            return self.minionprobs[barcode][get]

    def get_margin_minion_probability(self, barcode, get='norm_p_minion', e_s=0.051, e_i=0.049, e_d=0.078):
        '''get in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']'''
        if not self.margin_minion_probs:
            self._get_marginalized_minion_probabilities(e_s, e_i, e_d)
        elif [e_s, e_i, e_d] != self.minionprobs[barcode]['params']:
            self._get_marginalized_minion_probabilities(e_s, e_i, e_d)
        return self.margin_minion_probs[barcode][get]

    def get_binomial_probability(self, barcode, get=1, e_s=0.051, e_i=0.049, e_d=0.078):
        '''get in [0,1,2,3] for:
            0 = binom_prob_k_matches_in_alignment
            1 = binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode
            2 = binom_prob_k_matches_in_barcode
            3 = p_m used'''
        try:
            if 1-e_s-e_i-e_d != self.binomprobs[barcode][-1]:
                self._get_binomial_probability(barcode, e_s, e_i, e_d)
            return self.binomprobs[barcode][get]
        except: 
            self._get_binomial_probability(barcode, e_s, e_i, e_d)
            return self.binomprobs[barcode][get]
        outstring = 'binom_prob_k_matches_in_alignment:' + str(self.binomprobs[barcode][0]) + '\n'
        outstring += 'binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode:' + str(self.binomprobs[barcode][1]) + '\n'
        outstring += 'binom_prob_k_matches_in_barcode:' + str(self.binomprobs[barcode][2]) + '\n'

    def get_margin_binomial_probability(self, barcode, get=1, e_s=0.051, e_i=0.049, e_d=0.078):
        '''get in [0,1,2,3] for:
            0 = binom_prob_k_matches_in_alignment
            1 = binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode
            2 = binom_prob_k_matches_in_barcode
            3 = p_m used'''
        if not self.margin_binom_probs:
            self._get_marginalized_binomial_probabilities(e_s, e_i, e_d)
        elif 1-e_s-e_i-e_d != self.binomprobs[barcode][-1]:
            self._get_marginalized_binomial_probabilities(e_s, e_i, e_d)
        return self.margin_binom_probs[barcode][get]

    def get_edit_distance(self, barcode):
        try:
            return self.counts[barcode]['mm'] + self.counts[barcode]['d'] + self.counts[barcode]['i']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.counts[barcode]['mm'] + self.counts[barcode]['d'] + self.counts[barcode]['i']

    def get_edit_distance_with_query_unaln(self, barcode):
        try:
            return self.get_edit_distance(barcode) + self.counts[barcode]['u']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.get_edit_distance(barcode) + self.counts[barcode]['u']

    def get_edit_distance_with_ref_unaln(self, barcode):
        try:
            return self.get_edit_distance(barcode) + self.counts[barcode]['ur']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.get_edit_distance(barcode) + self.counts[barcode]['ur']

    def get_edit_distance_with_total_unaln(self, barcode):
        try:
            return self.get_edit_distance(barcode) + self.counts[barcode]['utotal']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.get_edit_distance(barcode) + self.counts[barcode]['utotal']

        
            
        
    def get_score_string(self):
        self._alignment_check()
        self._score_check()
        return (',').join([str(e) for e in self.scores])

    def get_score_probabilities_string(self):
        self._alignment_check()
        self._score_check()
        self._score_probability_check()
        return (',').join([str(self.score_probabilities[barcode]) for barcode in self.score_probabilities.keys()])


    def get_barcode_aln_string(self, barcode):
        try:
            return self.alnstrings[barcode]['query']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.alnstrings[barcode]['query']

    def get_read_aln_string(self, barcode):
        try:
            return self.alnstrings[barcode]['ref']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.alnstrings[barcode]['ref']

    def get_match_aln_string(self, barcode):
        try:
            return self.alnstrings[barcode]['sticks']
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            return self.alnstrings[barcode]['sticks']


    def get_formatted_pairwise_alignment(self, barcode, blocksize=100, e_s=0.051, e_i=0.049, e_d=0.078):
        return self._get_formatted_pairwise_alignment(barcode, blocksize, e_s, e_i, e_d)

    def get_all_formatted_pairwise_alignments(self, blocksize=100, e_s=0.051, e_i=0.049, e_d=0.078):
        return self._get_all_formatted_pairwise_alignments(blocksize, e_s, e_i, e_d)


    def get_smallest_barcode_length(self):
        barlength = float('inf')
        for name, seq in self.barcodes.iteritems():
            seqlen = len(seq)
            if seqlen < barlength:
                barlength = seqlen
        return barlength

    def get_maxbar_list(self):
        self._alignment_check()
        self._score_check()
        self._score_probability_check()
        return self.get_max_scoring_barcode(), self.get_max_score_probability(), self.get_max_score(), self.get_mean_score(), self.get_score_string(), self.get_score_probabilities_string(), self.alignments[self.maxscoringbarcode].q_pos, self.alignments[self.maxscoringbarcode].q_end, self.alignments[self.maxscoringbarcode].r_pos, self.alignments[self.maxscoringbarcode].r_end, self.alignments[self.maxscoringbarcode].cigar_str, self.score_powers_sum

    def get_subset_maxbar_list(self, get=[0,1]):
        ##get is just list of indexes
        ans = self.get_maxbar_list()
        return [ans[i] for i in get]

    #Private

    def _alignment_check(self):
        if self.alignments is None:
            self._get_alignments()
            
    def _score_check(self):
        if self.scores is None:
            self._get_max_scoring_barcode()

    def _score_probability_check(self):
        if self.score_powers is None:
            self._get_score_probabilities()
            


    def _get_alignments(self): #barcodes, sw, sequence, ref_name='', use_entire_barcode=True):
        ''' barcodes is dict of barcodes
            sequence is given -- usually first X bp of a read.
            use_entire_barcode -- if False, it will trim barcodes to length of shortest.'''
        self.alignments = {}
        if not self.use_entire_barcode:
            minbarlen = self.get_smallest_barcode_length()
        for barcode_name, barcode_seq in self.barcodes.iteritems():
            if not self.use_entire_barcode:
                barcode_seq = barcode_seq[:minbarlen]
            alignment = self.sw.align(query=barcode_seq, ref=self.search_seq, query_name=barcode_name, ref_name=self.ref_name)
            self.alignments[barcode_name] = alignment



    def _get_max_scoring_barcode(self):
        self._alignment_check()
        maxscore = float('-inf')
        scores = []
        for barcode, alignment in self.alignments.iteritems():
            scores.append( float(alignment.score) )

            if alignment.score > maxscore:
                maxscore = alignment.score
                maxbar = barcode
        self.scores = np.array(scores)
        self.maxscore = maxscore
        self.maxscoringbarcode = maxbar

    
    def _get_score_probabilities(self):
        self._alignment_check()
        self._score_check()
        self.score_powers = {barcode:np.e**alignment.score for barcode, alignment in self.alignments.iteritems()}
        self.score_powers_sum = float(sum([self.score_powers[barcode] for barcode in self.score_powers.keys()]))
        self.score_probabilities = {barcode:self.score_powers[barcode]/self.score_powers_sum for barcode in self.score_powers.keys()}


    def _compute_formatted_alignment_strings_for(self, barcode):
        self._alignment_check()
        self.counts[barcode] = {'m':0, 'mm':0, 'd':0, 'i':0, 'u':0, 'qbases':0, 'rbases':0}
        self.alnstrings[barcode] = {'ref':'', 'query':'', 'sticks':''}
        r_i = self.alignments[barcode].r_pos
        q_i = self.alignments[barcode].q_pos
        for e in self.alignments[barcode].cigar:
            if e[1] == 'M':
                newref = self.alignments[barcode].orig_ref[r_i:r_i+e[0]] ## call on these in forloop below
                newquer = self.alignments[barcode].orig_query[q_i:q_i+e[0]]
                self.alnstrings[barcode]['ref'] += newref
                self.alnstrings[barcode]['query'] += newquer
                r_i += e[0]
                q_i += e[0]
                for i in range(len(newref)):
                    if newref[i].upper() == newquer[i].upper():
                        self.alnstrings[barcode]['sticks'] += '|'
                        self.counts[barcode]['m'] += 1
                    else:
                        self.alnstrings[barcode]['sticks'] += ' '
                        self.counts[barcode]['mm'] += 1
            if e[1] == 'D':
                self.alnstrings[barcode]['ref'] += self.alignments[barcode].orig_ref[r_i:r_i+e[0]]
                self.alnstrings[barcode]['query'] += '-'*e[0]
                self.alnstrings[barcode]['sticks'] += ' '*e[0]
                r_i += e[0]
                self.counts[barcode]['d'] += e[0] #1
            elif e[1] == 'I':
                self.alnstrings[barcode]['ref'] += '-'*e[0]
                self.alnstrings[barcode]['query'] += self.alignments[barcode].orig_query[q_i:q_i+e[0]]
                self.alnstrings[barcode]['sticks'] += ' '*e[0]
                q_i += e[0]
                self.counts[barcode]['i'] += e[0] #1

        self.counts[barcode]['qbases'] = self.counts[barcode]['m'] + self.counts[barcode]['mm'] + self.counts[barcode]['i']
        self.counts[barcode]['rbases'] = self.counts[barcode]['m'] + self.counts[barcode]['mm'] + self.counts[barcode]['d']
        self.counts[barcode]['refLen'] = len(self.alignments[barcode].orig_ref)
        self.counts[barcode]['queryLen'] = len(self.alignments[barcode].orig_query)
        
        self.counts[barcode]['u'] = self.counts[barcode]['queryLen'] - self.counts[barcode]['qbases'] ## these were bases not in the alignment. Since barcode is query... i.e. pieces of barcode not found in read
        self.counts[barcode]['ur'] = self.counts[barcode]['refLen'] - self.counts[barcode]['rbases']
        self.counts[barcode]['utotal'] = self.counts[barcode]['u'] + self.counts[barcode]['ur'] 
        self.counts[barcode]['alnlen'] = sum([self.counts[barcode][count] for count in ('m', 'mm', 'i', 'd')])
        self.counts[barcode]['PercentIdentity'] = 100.0*self.counts[barcode]['m']/self.counts[barcode]['alnlen']
        self.counts[barcode]['PercentIdentity_with_unaligned'] = 100.0*self.counts[barcode]['m']/sum([self.counts[barcode]['alnlen'], self.counts[barcode]['u']])
        self.counts[barcode]['PercentRbasesAligned'] = 100.0*self.counts[barcode]['rbases']/self.counts[barcode]['refLen']
        self.counts[barcode]['PercentQbasesAligned'] = 100.0*self.counts[barcode]['qbases']/self.counts[barcode]['queryLen']
        assert self.counts[barcode]['qbases'] == self.alignments[barcode].q_end-self.alignments[barcode].q_pos
        assert self.counts[barcode]['rbases'] == self.alignments[barcode].r_end-self.alignments[barcode].r_pos
        assert len(self.alnstrings[barcode]['ref']) == len(self.alnstrings[barcode]['query'])
        assert len(self.alnstrings[barcode]['ref']) == self.counts[barcode]['alnlen']

    def _compute_all_formatted_alignment_strings(self):
        for barcode in self.barcodes.keys():
            self._compute_formatted_alignment_strings_for(barcode)

    def _aln_string_check(self, barcode):
        try:
            self.counts[barcode]
            self.alnstrings[barcode]
        except:
            self._compute_formatted_alignment_strings_for(barcode)
            
    def _get_expected_value(self, barcode):
        self._alignment_check()
        self.expectedvalues[barcode] = {'n':0, 'p':0, 'evalue':0}
        self.expectedvalues[barcode]['n'] = len(self.alignments[barcode].orig_ref)*len(self.alignments[barcode].orig_query)
        self.expectedvalues[barcode]['p'] = np.e**(-1*self.alignments[barcode].score)
        self.expectedvalues[barcode]['evalue'] = self.expectedvalues[barcode]['n'] * self.expectedvalues[barcode]['p']
        

    def _get_expected_values(self):
        minscore = float('inf')
        for barcode in self.barcodes.keys():
            self._get_expected_value(barcode)
            if self.expectedvalues[barcode]['evalue'] < minscore:
                minscore = self.expectedvalues[barcode]['evalue']
                minbar = barcode
        self.lowestevalue = minscore
        self.lowestevaluebar = minbar
            

    def _get_formatted_pwaln_header(self, barcode):
        self._aln_string_check(barcode)
        try:
            return self.alnstrings[barcode]['header']
        except:
            outstring = barcode + ' to '  + self.alignments[barcode].r_name + '\n'
            outstring += "Ref (top): " + self.alignments[barcode].r_name + ' ' + str(self.alignments[barcode].r_pos) + '-' + str(self.alignments[barcode].r_end) + ' r_bases_aligned:' + str(self.counts[barcode]['rbases']) + ' pct_r_bases_aligned:' + str(self.counts[barcode]['PercentRbasesAligned']) + ' refLen:' + str(self.counts[barcode]['refLen']) + ' bp\n'
            outstring += "Query (bottom): " + self.alignments[barcode].q_name + ' ' + str(self.alignments[barcode].q_pos) + '-' + str(self.alignments[barcode].q_end) + ' q_bases_aligned:' + str(self.counts[barcode]['qbases']) + ' pct_q_bases_aligned:'+str(self.counts[barcode]['PercentQbasesAligned']) + ' queryLen:' + str(self.counts[barcode]['queryLen']) + ' bp\n'
            stats = ["AS:" + str(self.alignments[barcode].score), "Match:"+str(self.counts[barcode]['m']), "Mismatch:"+str(self.counts[barcode]['mm']), "Deletion:"+str(self.counts[barcode]['d']), "Insertion:"+str(self.counts[barcode]['i']), "AlignmentLength:"+str(self.counts[barcode]['alnlen']), "PercentIdentity:"+str(self.counts[barcode]['PercentIdentity']),  "Barcode_Unaligned:"+str(self.counts[barcode]['u']), "PercentIdentity_with_unaligned:"+str(self.counts[barcode]['PercentIdentity_with_unaligned'])]
            outstring += (' ').join(stats) + '\n'
            self.alnstrings[barcode]['header'] = outstring
            return self.alnstrings[barcode]['header']

    def _get_formatted_evalue_string(self, barcode):
        stats = ['n:' + str(self.get_expected_value(barcode, 'n')), 'p:' + str(self.get_expected_value(barcode, 'p')), 'e_value:' + str(self.get_expected_value(barcode, 'evalue'))]
        outstring = (' ').join(stats) + '\n'
        return outstring

    def _get_minion_probability(self, barcode, e_s=0.051, e_i=0.049, e_d=0.078):
        ## since the read is the "ref" here.... n_d and n_i are dels/ins from/in read.
        ## Barcode is the real "reference" when considering 'deletions' and 'insertions' in MinION read
        ## ...meaning for minion probs the number of insertions is what we have listed as numdels
        ## .... and number of dels is what we have listed as num insertions
        ## i.e. we need e_i**n_d and e_d**n_i instead of e_i**n_i and e_d**n_d
        ## NOTE2: ## p_minion not necessarily comparable when barcodes are different lengths - can divide by barcode_len or q*r maybe -- or can multiply by multinomial maybe...
        self._aln_string_check(barcode)
        ## Gather variables needed
        p_m = 1 - e_i - e_d - e_s
        Esum = e_s + e_i + e_d
        E_s = e_s/Esum
        E_i = e_i/Esum
        E_d = e_d/Esum
        N = self.counts[barcode]['alnlen']
        N2 = N - self.counts[barcode]['m']
        N3 = N2 - self.counts[barcode]['mm']
        N4 = N3 - self.counts[barcode]['i']
        N5 = N4 - self.counts[barcode]['d']
        NcK = nchoosek(N, self.counts[barcode]['m']) * nchoosek(N2, self.counts[barcode]['mm']) * nchoosek(N3, self.counts[barcode]['i']) * nchoosek(N4, self.counts[barcode]['d']) 
        N_u = N + self.counts[barcode]['u']
        #populate dict
        self.minionprobs[barcode] = {}
        self.minionprobs[barcode]['params'] = [e_s, e_i, e_d]
        self.minionprobs[barcode]['p_minion_aln'] =  (p_m**self.counts[barcode]['m']) * (e_s**self.counts[barcode]['mm']) * (e_d**self.counts[barcode]['i']) * (e_i**self.counts[barcode]['d']) ## See note above as to why i and d are switched here
        self.minionprobs[barcode]['p_minion_un'] = (e_s**(E_s*self.counts[barcode]['u'])) * (e_i**(E_i*self.counts[barcode]['u'])) *(e_d**(E_d*self.counts[barcode]['u'])) ## since it is unclear if the unaligned were subs/dels/ins I am making use of all
        self.minionprobs[barcode]['p_minion'] = self.minionprobs[barcode]['p_minion_aln'] * self.minionprobs[barcode]['p_minion_un']
        self.minionprobs[barcode]['norm_p_minion_aln'] = NcK * self.minionprobs[barcode]['p_minion_aln']
        self.minionprobs[barcode]['norm_p_minion'] = self.minionprobs[barcode]['norm_p_minion_aln'] * self.minionprobs[barcode]['p_minion_un'] ## There seems to be no reason to multiply the unaligned by a nchoosek b/c they are fixed at the ends and correspond to some "single" unknown composite error



    def _get_minion_probabilities(self, e_s=0.051, e_i=0.049, e_d=0.078):
        self.minion_prob_sums = {k:0 for k in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']}
        for barcode in self.barcodes.keys():
            self._get_minion_probability(barcode, e_s, e_i, e_d)
            for k in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']:
                self.minion_prob_sums[k] += self.minionprobs[barcode][k]
    
    def _get_marginalized_minion_probabilities(self, e_s=0.051, e_i=0.049, e_d=0.078):
        self._get_minion_probabilities(e_s, e_i, e_d)
        maxscore = {k:float('-inf') for k in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']}
        maxbar = {}
        for barcode in self.barcodes.keys():
            self.margin_minion_probs[barcode] = {}
            for k in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']:
                self.margin_minion_probs[barcode][k] = self.minionprobs[barcode][k]/self.minion_prob_sums[k]
                if self.margin_minion_probs[barcode][k] > maxscore[k]:
                    maxscore[k] = self.margin_minion_probs[barcode][k]
                    maxbar[k] = barcode
        self.maxminion = maxscore
        self.maxminionbar = maxbar
        
    def _marg_minion_check(self, e_s=0.051, e_i=0.049, e_d=0.078, force_redo=False):
        if not self.margin_minion_probs or force_redo:
            self._get_marginalized_minion_probabilities(e_s, e_i, e_d)

    def _get_formatted_minionprob_string(self, barcode, e_s=0.051, e_i=0.049, e_d=0.078):
        ''' '''
        stats = [k+':'+str(self.get_minion_probability(barcode, get=k, e_s=e_s, e_i=e_i, e_d=e_d)) for k in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']] 
        outstring = (' ').join(stats) + '\n'
        return outstring

    def _get_formatted_margminionprob_string(self, barcode, e_s=0.051, e_i=0.049, e_d=0.078):
        ''' '''
        self._marg_minion_check(e_s, e_i, e_d)
        stats = ['marg_'+k+':'+str(self.get_margin_minion_probability(barcode, get=k, e_s=e_s, e_i=e_i, e_d=e_d)) for k in ['p_minion_aln', 'p_minion_un', 'p_minion', 'norm_p_minion_aln', 'norm_p_minion']] 
        outstring = (' ').join(stats) + '\n'
        return outstring
    
    def _get_binomial_probability(self, barcode, e_s=0.051, e_i=0.049, e_d=0.078):
        p_m = 1 - e_i - e_d - e_s
        n_m = self.counts[barcode]['m']
        n_mm = self.counts[barcode]['mm']+self.counts[barcode]['i']+self.counts[barcode]['d']
        n_u = self.counts[barcode]['u']
        #binom_prob_k_matches_in_alignment
        binom_prob1 = nchoosek(self.counts[barcode]['alnlen'], n_m) * (p_m**n_m) * ((1-p_m)**(n_mm))
        #binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode
        binom_prob2 = nchoosek(self.counts[barcode]['alnlen']+n_u, n_m) * (p_m**n_m) * ((1-p_m)**(n_mm+n_u))
        #binom_prob_k_matches_in_barcode (w/ unaligned parts)
        binom_prob3 = nchoosek(self.counts[barcode]['queryLen'], n_m) * (p_m**n_m) * ((1-p_m)**(self.counts[barcode]['queryLen']-n_m))
        self.binomprobs[barcode] = [binom_prob1, binom_prob2, binom_prob3, p_m]
        

    def _get_binomial_probabilities(self, e_s=0.051, e_i=0.049, e_d=0.078):
        self.binomial_prob_sums = {k:0 for k in [0,1,2]}
        for barcode in self.barcodes.keys():
            self._get_binomial_probability(barcode, e_s, e_i, e_d)
            for k in [0,1,2]:
                self.binomial_prob_sums[k] += self.binomprobs[barcode][k]
    
    def _get_marginalized_binomial_probabilities(self, e_s=0.051, e_i=0.049, e_d=0.078):
        self._get_binomial_probabilities(e_s, e_i, e_d)
        maxscore = {0:float('-inf'), 1:float('-inf'), 2:float('-inf')}
        maxbar = {}
        for barcode in self.barcodes.keys():
            self.margin_binom_probs[barcode] = {}
            for k in [0,1,2]:
                self.margin_binom_probs[barcode][k] = self.binomprobs[barcode][k]/self.binomial_prob_sums[k]
                if self.margin_binom_probs[barcode][k] > maxscore[k]:
                    maxscore[k] = self.margin_binom_probs[barcode][k]
                    maxbar[k] = barcode
        self.maxbinom = maxscore
        self.maxbinombar = maxbar

    def _marg_binom_check(self, e_s=0.051, e_i=0.049, e_d=0.078, force_redo=False):
        if not self.margin_binom_probs or force_redo:
            self._get_marginalized_binomial_probabilities(e_s, e_i, e_d)

    def _get_formatted_binomprob_string(self, barcode, e_s=0.051, e_i=0.049, e_d=0.078):
        outstring = 'binom_prob_k_matches_in_alignment:' + str(self.get_binomial_probability(barcode, get=0, e_s=0.051, e_i=0.049, e_d=0.078)) + '\n'
        outstring += 'binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode:' + str(self.get_binomial_probability(barcode, get=1, e_s=0.051, e_i=0.049, e_d=0.078)) + '\n'
        outstring += 'binom_prob_k_matches_in_barcode:' + str(self.get_binomial_probability(barcode, get=2, e_s=0.051, e_i=0.049, e_d=0.078)) + '\n'
        return outstring

    def _get_formatted_margbinomprob_string(self, barcode, e_s=0.051, e_i=0.049, e_d=0.078):
        self._marg_binom_check(e_s, e_i, e_d)
        outstring = 'marg_binom_prob_k_matches_in_alignment:' + str(self.get_margin_binomial_probability(barcode, get=0, e_s=0.051, e_i=0.049, e_d=0.078)) + '\n'
        outstring += 'marg_binom_prob_k_matches_in_alignment_incl_unaligned_portions_of_barcode:' + str(self.get_margin_binomial_probability(barcode, get=1, e_s=0.051, e_i=0.049, e_d=0.078)) + '\n'
        outstring += 'marg_binom_prob_k_matches_in_barcode:' + str(self.get_margin_binomial_probability(barcode, get=2, e_s=0.051, e_i=0.049, e_d=0.078)) + '\n'
        return outstring

    def _get_formatted_pairwise_alignment(self, barcode, blocksize=100, e_s=0.051, e_i=0.049, e_d=0.078):
        ## TODO: printout unaligned portion as part of alignment viz
        ## e_i, e_d, and e_s are insertion/deletion/substitution errors found in early MinION sequencing by Jain et al: Improved data analysis for the MinION nanopore sequencer
        ## Can compute prob of alignment by p_m=1-e_i-e_d-e_s;
        outstring = self._get_formatted_pwaln_header(barcode)
        outstring += self._get_formatted_evalue_string(barcode)
        outstring += self._get_formatted_minionprob_string(barcode)
        outstring += self._get_formatted_margminionprob_string(barcode)
        outstring += self._get_formatted_binomprob_string(barcode)
        outstring += self._get_formatted_margbinomprob_string(barcode)
        outstring += "Marginalized Probability Given Barcode Set: " + str(self.get_score_probability(barcode)) + "\n"
        for i in range(0, self.counts[barcode]['alnlen'], blocksize):
            outstring += self.get_read_aln_string(barcode)[i:i+blocksize] + '\n'
            outstring += self.get_match_aln_string(barcode)[i:i+blocksize] + '\n'
            outstring += self.get_barcode_aln_string(barcode)[i:i+blocksize] + '\n\n'

        return outstring

    def _get_all_formatted_pairwise_alignments(self, blocksize=100, e_s=0.051, e_i=0.049, e_d=0.078):
        outstring = ''
        try:
            self.lowestevaluebar
        except:
            self._get_expected_values()
        for barcode in self.barcodes.keys():
            outstring += self._get_formatted_pairwise_alignment(barcode, blocksize, e_s, e_i, e_d)
        outstring += "Highest Marginalized Score Probability: " + self.maxscoringbarcode + ' ' + str(self.score_probabilities[self.maxscoringbarcode]) + "\n\n"
        outstring += "Lowest E-value: " + self.lowestevaluebar + ' ' + str(self.lowestevalue) + "\n\n"

        outstring += "Highest Marginialized MinION Probability of alignment: " + self.maxminionbar['p_minion_aln'] + ' ' + str(self.maxminion['p_minion_aln']) + "\n"
        outstring += "Highest Marginialized MinION Probability of alignment (including unaligned): " + self.maxminionbar['p_minion'] + ' ' + str(self.maxminion['p_minion']) + "\n"
        outstring += "Highest Marginialized MinION Probability of alignment (scaled by multinomial): " + self.maxminionbar['norm_p_minion_aln'] + ' ' + str(self.maxminion['norm_p_minion_aln']) + "\n"
        outstring += "Highest Marginialized MinION Probability of alignment (including unaligned; scaled by multinomial): " + self.maxminionbar['norm_p_minion'] + ' ' + str(self.maxminion['norm_p_minion']) + "\n\n"

        outstring += "Highest Marginialized Binomial Probability of k matches in alignment: " + self.maxbinombar[0] + ' ' + str(self.maxbinom[0])  + "\n"
        outstring += "Highest Marginialized Binomial Probability of k matches in alignment (including unaligned): " + self.maxbinombar[1] + ' ' + str(self.maxbinom[1])  + "\n"
##        outstring += "Highest Marginialized Binomial Probability of k matches in barcode: " + self.maxbinombar[2] + ' ' + str(self.maxbinom[2])  + "\n\n"
        return outstring                                                           
                                                              
    def _get_highest_barcode(self, barcodedict):
        '''Assumes everything needed is min place'''
        maxscore = float('-inf')
        for barcode, score in barcodedict:
            if score > maxscore:
                maxscore = score
                maxbar = barcode
        return maxbar, maxscore

    def _get_highest_barcode_string(self, barcodedict, label=''):
        ans = self.get_highest_barcode(barcodedict)
        return 'Highest ' + label + ': ' + str(ans)

    def _get_lowest_barcode(self, barcodedict):
        '''Assumes everything needed is min place'''
        minscore = float('inf')
        for barcode, score in barcodedict:
            if score < minscore:
                minscore = score
                minbar = barcode
        return minbar, minscore
    
    def _get_lowest_barcode_string(self, barcodedict, label=''):
        ans = self.get_lowest_barcode(barcodedict)
        return 'Lowest ' + label + ': ' + str(ans)
    
        # could also give cigar: alignments[maxbar].cigar_str
    ##    D = sum([int(e) for e in re.findall('(\d)D',alignments[maxbar].cigar_str)])
    ##    I = sum([int(e) for e in re.findall('(\d)I',alignments[maxbar].cigar_str)])
    ##    M = sum([int(e) for e in re.findall('(\d)M',alignments[maxbar].cigar_str)])
    ##    print "FINAL:"
    ##    get_formatted_pairwise_alignment(alignments[maxbar])
        ##  alignments[maxbar].identity, alignments[maxbar].matches, alignments[maxbar].mismatches, I, D, M,
##        return maxbar, maxprob, maxscore, meanscore, scorestring, probs_string, alignments[maxbar].q_pos, alignments[maxbar].q_end, alignments[maxbar].r_pos, alignments[maxbar].r_end, alignments[maxbar].cigar_str, epowersum #, alignments[maxbar].orig_query,  alignments[maxbar].orig_ref 



