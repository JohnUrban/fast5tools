import itertools
import numpy as np
import time
from collections import defaultdict
from profilehooks import profile
from scipy.stats import norm


from fast5tools.tools import *

STATES_ONEMERS = [''.join(e) for e in itertools.product("ACGT")]
STATES_DIMERS = [''.join(e) for e in itertools.product("ACGT","ACGT")]
STATES_TRIMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT")]
STATES_FOURMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT")]
STATES_FIVEMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT")]
STATES_SIXMERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT","ACGT")]




def get_emiss_probs_from_model(model, twoemits=False):
    ''' model is object returned from get_stored_model() in model_tools '''
    states = sorted(model[1].keys())
    num_states = len(states)
    t_emissions = np.zeros([2,num_states])
    c_emissions = np.zeros([2,num_states])
    if twoemits:
        t_emissions2 = np.zeros([2,num_states])
        c_emissions2 = np.zeros([2,num_states])
    for i in range(num_states):
        t_emissions[0,i] = model[1][states[i]][0]
        t_emissions[1,i] = model[1][states[i]][1]
        c_emissions[0,i] = model[2][states[i]][0]
        c_emissions[1,i] = model[2][states[i]][1]
        if twoemits:
            t_emissions2[0,i] = model[1][states[i]][2]
            t_emissions2[1,i] = model[1][states[i]][3]
            c_emissions2[0,i] = model[2][states[i]][2]
            c_emissions2[1,i] = model[2][states[i]][3]
    if twoemits:
        return t_emissions, c_emissions, t_emissions2, c_emissions2
    return t_emissions, c_emissions


class HMM(object):

    def __init__(self, states=None, emissions=None):
        # data
        self.states = states
        self.statepath = None
        self.initial_probs = None
        self.transition_probs = None
        self.emission_probs = None
        self.emissions = None
        self.num_states = None
        self.num_emits = None
        if emissions is None: # initialize with random values
            if self.states is None:
                self.states = STATES_FIVEMERS
            self.randomize_emission_probs()
            self.randomize_transition_probs()
            self.randomize_initial_probs()
            self.randomize_statepath()
            self.randomize_emissions()
        else:
            self.emissions = emissions
            

    def add_initial_probs(self, initial_probs):
        ## should be 1 x nstates np.array
        assert len(initial_probs) == len(self.states)
        self.initial_probs = initial_probs
        

    def add_transition_probs(self, transition_probs):
        ## must be nstate x nstate np.array
        self.transition_probs = transition_probs

    def add_emission_probs(self, emission_probs):
        ## should be be 2 x nstate np.array (if normal -- need means and stdevs)
        self.emission_probs = emission_probs

    def add_states(self, states):
        self.states = states

    def get_initial_probs(self):
        return self.initial_probs

    def get_transition_probs(self):
        return self.transition_probs

    def get_emission_probs(self):
        return self.emission_probs

    def randomize_initial_probs(self, uniform=True):
        nstates = len(self.states)
        if uniform:
            self.initial_probs = [1.0 / nstates] * nstates
        else:
            self.initial_probs = np.random.poisson(lam=10.0, size=nstates)
            self.initial_probs /= sum(self.initial_probs)

    def randomize_transition_probs(self, allow_gaps=True):
        """
        If allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
        """
        k = len(self.states[0])
        if k > 2:
            nstates = len(self.states)
            self.transition_probs = np.zeros([nstates, nstates])
            # make prefix-suffix dict -- overlaps of k-1 and k-2
            prefix = defaultdict(list)
            for i in xrange(nstates):
                pref = self.states[i][:k-1]
                prefix[pref].append(i)
                pref = self.states[i][:k-2]
                prefix[pref].append(i)
            # create transition probs -- can soft code the sampling parameters later if want
            for i in xrange(nstates):
                ## overlap by k-1 (move = 1)
                current_suffix = self.states[i][1:]
                poisson = np.random.poisson(lam=365.0, size=k-1)
                for t, j in enumerate(prefix[current_suffix]):
                    self.transition_probs[i,j] = poisson[t]
                if allow_gaps:
                    ## overlap by k-2 (move = 2) -- add additional counts
                    current_suffix = self.states[i][2:]
                    poisson = np.random.poisson(lam=4.0, size=(k-1)**2)
                    for t, j in enumerate(prefix[current_suffix]):
                        self.transition_probs[i,j] += poisson[t]
                    ## stay in place: add additional probability to staying in place (move = 0)
                    current_suffix = self.states[i]
                    poisson = np.random.poisson(lam=20.0, size=1)
                    self.transition_probs[i,i] += poisson
                ## normalize all counts by sum to create probs that sum to 1
                self.transition_probs[i,:] /= sum(self.transition_probs[i,:])

    def randomize_emission_probs(self, level=True):
        ## generates either level emissions or sd emissions
        # mu.mean and sigma.mean are the mean and std dev of the r7.3 state level means to be used to generate emission means 
        # mu.sd, sigma.sd -- same for std devs of signals
        if level:
            mu_mean = 65.57454
            sigma_mean = 6.497453
            mu_sd = 1.163836
            sigma_sd = 0.4116285
        else: ## sd emission
            mu_mean = 1.37316
            sigma_mean = 0.3144043
            mu_sd = 0.1761904
            sigma_sd = 0.06263217
        nstates = len(self.states)
        self.emission_probs = np.zeros([2, nstates])
        for i in xrange(nstates):
            self.emission_probs[0,i] = np.random.normal(mu_mean, sigma_mean)
            self.emission_probs[1,i] = abs(np.random.normal(mu_sd, sigma_sd))

    def randomize_statepath(self, length=10):
        nstates = len(self.states)
        self.statepath = [np.random.choice(nstates, p=self.initial_probs)]
        for _ in xrange(length-1):
            self.statepath.append(
                    np.random.choice(
                            nstates,
                            p=self.transition_probs[self.statepath[-1]]))


    def randomize_emissions(self):
        means = self.emission_probs[0, self.statepath]
        stdevs = self.emission_probs[1, self.statepath]
        self.emissions = np.random.normal(means, stdevs)


    def randomize_emissions_twoemits(self):
        pass


    def get_num_states(self):
        if self.num_states == None:
            self.num_states = len(self.states)
        return self.num_states

    def get_num_emits(self):
        if self.num_emits == None:
            self.num_emits = len(self.emissions)
        return self.num_emits
    
    def forward(self):
        nstates = self.get_num_states()
        nemits = self.get_num_emits()
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        Forward = np.zeros([nstates, nemits])
        scalefactors = np.zeros([2, nemits])
        # initial
        Forward[:, 0] = np.multiply(self.initial_probs, ep.pdf(self.emissions[0]))
        # scale to prevent underflow -- keep track of scaling
        scalefactors[0,0] = sum(Forward[:,0])
        scalefactors[1,0] = np.log(scalefactors[0,0])
        Forward[:,0] /= scalefactors[0,0]
        # iterate
        for k in xrange(1, nemits):
            emit = ep.pdf(self.emissions[k])
            Forward[:,k] = np.multiply(emit, np.dot(Forward[:,k-1], self.transition_probs))
            scalefactors[0,k] = sum(Forward[:,k])
            scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k-1]
            Forward[:,k] /= scalefactors[0,k]
        return Forward, scalefactors

    def backward(self):
        nstates = self.get_num_states()
        nemits = self.get_num_emits()
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        Backward = np.zeros([nstates, nemits])
        scalefactors = np.zeros([2, nemits])
        end = nemits - 1
        # initial
        Backward[:, end] = 1
        # scale to prevent underflow -- keep track of scaling
        scalefactors[0,end] = sum(Backward[:,end])
        scalefactors[1,end] = np.log(scalefactors[0,end])
        Backward[:,end] /= scalefactors[0,end]
        # iterate
        for k in xrange(end-1, -1, -1):
            emit = ep.pdf(self.emissions[k+1])
            a = np.multiply(Backward[:,k+1], emit).transpose()
            Backward[:,k] = np.dot(self.transition_probs, a).transpose()
            scalefactors[0,k] = sum(Backward[:,k])
            scalefactors[1,k] = np.log(scalefactors[0,k]) + scalefactors[1,k+1]
            Backward[:,k] /= scalefactors[0,k]
        return Backward, scalefactors

    def posterior_decoding(self, Forward, F_scales, Backward, B_scales):
        ##F and B are scaled long seq matrices -- the scales are scalefactors that come with them out of long fxns
        nstates = self.get_num_states()
        nemits = np.shape(Forward)[1]
        posterior_path = np.zeros(nemits, dtype=int)
        for i in xrange(nemits):
            fb = Forward[:,i] * Backward[:,i]
            posterior_path[i] = int(fb.argmax())
        return posterior_path

    def prob_data(Forward, scalefactors, nemits=None):
        if nemits == None:
            end = np.shape(Forward)[1]-1
        else:
            end = nemits-1
        return sum(Forward[:,end])*np.exp(scalefactors[1,end])

##    @profile

    def viterbi(self):
        np.seterr(divide='ignore')
        nstates = self.get_num_states()
        nemits = self.get_num_emits()
        initial_probs = np.log(self.initial_probs)
        tran_probs = np.log(self.transition_probs)
        ep = norm(self.emission_probs[0,:], self.emission_probs[1,:])
        pointer = np.zeros([nemits, nstates])
        Viterbi = np.zeros([nstates, nemits])  
        ## need to add log_probs instead of multiply probs to prevent underflow
        Viterbi[:,0] = self.initial_probs + ep.logpdf(self.emissions[0])
        pointer[0,:] = 1
        for j in range(1,nemits):
            selection = Viterbi[:,j-1] + tran_probs.transpose() 
            maxstates = np.apply_along_axis(max_and_index, 1, selection)
            Viterbi[:,j] = ep.logpdf(self.emissions[j]) + maxstates[:,1]
            pointer[j,:] = maxstates[:,0]
        end = nemits - 1
        #path init
        viterbi_path = np.zeros(nemits).astype(int)
        viterbi_path[end] = Viterbi[:,end].argmax()
        #prob
        viterbi_prob = Viterbi[viterbi_path[end], end]
        #path iter
        for j in range(end,0,-1):
            viterbi_path[j-1] = pointer[j,viterbi_path[j]]
        return viterbi_path, viterbi_prob


##    def viterbi2(self):
##        np.seterr(divide='ignore')
##        num_states = self.get_num_states()
##        num_emits = self.get_num_emits()
##        initial_probs = np.log(self.initial_probs)
##        tran_probs = np.log(self.transition_probs)
##        ep1 = norm(emission_probs[0,:], emission_probs[1,:])
##        ep2 = norm(emission_probs2[0,:], emission_probs2[1,:])
##        pointer = np.zeros([num_emits, num_states])
##        Viterbi = np.zeros([num_states, num_emits])  
##        ## need to add log_probs instead of multiply probs to prevent underflow
##        Viterbi[:,0] = initial_probs + ep1.logpdf(emitted_data[0]) + ep2.logpdf(emitted_data2[0])
##        pointer[0,:] = 1
##        for j in range(1,num_emits):
##            selection = Viterbi[:,j-1] + tran_probs.transpose() 
##            maxstates = np.apply_along_axis(max_and_index, 1, selection)
##            Viterbi[:,j] = ep1.logpdf(emitted_data[j]) + ep2.logpdf(emitted_data2[j]) + maxstates[:,1]
##            pointer[j,:] = maxstates[:,0]
##        end = num_emits - 1
##        #path init
##        viterbi_path = np.zeros(num_emits).astype(int)
##        viterbi_path[end] = Viterbi[:,end].argmax()
##        #prob
##        viterbi_prob = Viterbi[viterbi_path[end], end]
##        #path iter
##        for j in range(end,0,-1):
##            viterbi_path[j-1] = pointer[j,viterbi_path[j]]
##        return viterbi_path, viterbi_prob
##

    def compare_statepath(self, dst, src=None):
        if src is None:
            src = self.statepath
        ident = sum(a==b for a,b in itertools.izip(dst, src))
        edit_dist = len(src) - ident
        return edit_dist, ident, 100.0*ident/len(src)


    def baumwelch():
        pass

    def nwalign(self, s2, s1=None):
        return nw.global_align(s1, s2)
    def compare_seq_nwa(s1,s2):
        s1, s2 = nwalign(s1,s2)
        length = len(s1)
        dist = edit_dist(s1,s2,length)
        return dist, pct_id(length,dist)

    def combine_2_seq(s1,s2, length=None):
        '''Assumes length s1 == length s2 '''
        s1,s2 = nwalign(s1,s2)
        if length == None:
            length = len(s1)
        editdist = 0
        combinedseq = ''
        for i in range(length):
            if s1[i] == s2[i]:
                combinedseq += s1[i]
            elif s1[i] == "-":
                editdist += 1
                combinedseq += s2[i]
            elif s2[i] == "-":
                editdist += 1
                combinedseq += s1[i]
            else: ## mismatch -- arbitrarily go with complement
                editdist += 1
                combinedseq += s1[i]
        return combinedseq, editdist

    def compare_seq_nwa(s1,s2):
        s1, s2 = nwalign(s1,s2)
        length = len(s1)
        dist = edit_dist(s1,s2,length)
        return dist, pct_id(length,dist)

    def get_2D_seq(t,c):
        c = complement(c)
        return combine_2_seq(t,c)

    def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
        ## states are some type of kmer
        ## statepath is vector of numbers (indexes)
        path_length = len(statepath)
        moves = [0]*path_length ## first move is 0
        k = len(states[0])
        end = k-1
        if k == 1 or k == 2:
            return "This currently only works with 3-mers as smallest kmer."
        else:
            #init
            seq = states[statepath[0]]
            moves[0] = 0
            #iter
            for i in range(1,path_length):
                lastSuffix = states[statepath[i-1]][1:]
                currentPrefix = states[statepath[i]][:k-1]
                if statepath[i-1] == statepath[i]:
                    moves[i] = 0
                elif lastSuffix == currentPrefix:
                    seq += states[statepath[i]][end]
                    moves[i] = 1
                else:
                    lastSuffix = states[statepath[i-1]][2:]
                    currentPrefix = states[statepath[i]][:k-2]
                    if lastSuffix == currentPrefix:
                        seq += states[statepath[i]][end-1:]
                        moves[i] = 2
                    elif posterior_decoded:
                        seq += states[statepath[i]][end]
                        moves[i] = -1 
                        ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
                        ## it turns out that adding the base from the illegal move does not hurt the seq overall much
        return seq, moves

    ### Allow longer gaps(skips)
##    def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
##        ## states are some type of kmer
##        ## statepath is vector of numbers (indexes)
##        path_length = len(statepath)
##        moves = [0]*path_length ## first move is 0
##        k = len(states[0])
##        end = k-1
##        if k == 1 or k == 2:
##            return "This currently only works with 3-mers as smallest kmer."
##        else:
##            #init
##            seq = states[statepath[0]]
##            moves[0] = 0
##            #iter
##            for i in range(1,path_length):
##                lastSuffix = states[statepath[i-1]][1:]
##                currentPrefix = states[statepath[i]][:k-1]
##                if lastSuffix == currentPrefix:
##                    seq += states[statepath[i]][end]
##                    moves[i] = 1
##                elif statepath[i-1] == statepath[i]:
##                    ## by checking same state last, only heteropolymers affected
##                    ## homopolymers would be caught in first condition
##                    moves[i] = 0
##                    ## nothing is added to sequence
##                    ## could make another fxn that just spits out events and states line by line like 'template events' in f5
##         
##                else:
##                    lastSuffix = states[statepath[i-1]][2:]
##                    currentPrefix = states[statepath[i]][:k-2]
##                    if lastSuffix == currentPrefix:
##                        seq += states[statepath[i]][end-1:]
##                        moves[i] = 2
##                    else:
##                        lastSuffix = states[statepath[i-1]][3:]
##                        currentPrefix = states[statepath[i]][:k-3]
##                        if lastSuffix == currentPrefix:
##                            seq += states[statepath[i]][end-2:]
##                            moves[i] = 3
##                        else:
##                            lastSuffix = states[statepath[i-1]][4:]
##                            currentPrefix = states[statepath[i]][:k-4]
##                            if lastSuffix == currentPrefix:
##                                seq += states[statepath[i]][end-3:]
##                                moves[i] = 4
##                            else:
##                                ## skip 5
##                                seq += states[statepath[i]][end-4:]
##                                moves[i] = 5
##                       ## ELSE::: do what? ... in other one just added centroid seq regardless...
##    ##                elif posterior_decoded:
##    ##                    seq += states[statepath[i]][end]
##    ##                    moves[i] = -1 
##    ##                    ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
##                        ## it turns out that adding the base from the illegal move does not hurt the seq overall much
##        return seq, moves






##############################END HMM CLASS

##def generate_random_kmer_transition_probs(states, allow_gaps=True, unif=False):
##    ## if allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
##    ## can set nonzero.trans to any vector -- for DNA length 4
##    k = len(states[0])
##    if k == 1:
##        pass
##    elif k == 2:
##        pass
##    else:
##        num_states = len(states)
##        tran_probs = np.zeros([num_states,num_states])
##
##        # make prefix-suffix dict -- overlaps of k-1 and k-2
##        prefix = defaultdict(list)
##        for i in range(num_states):
##            pref = states[i][:k-1]
##            prefix[pref].append(i)
##            pref = states[i][:k-2]
##            prefix[pref].append(i)
##
##        ## create transition probs -- can soft code the sampling parameters later if want
##        for i in range(num_states):
##            ## overlap by k-1 (move = 1)
##            current_suffix = states[i][1:k]
##            if unif:
##                trans = np.array([365,365,365,365])
##            else:
##                trans = np.random.poisson(lam=365.0,size=4)
##            t = 0 
##            for j in prefix[current_suffix]:
##                tran_probs[i,j] = trans[t]
##                t += 1
##            if allow_gaps:
##                ## overlap by k-2 (move = 2) -- add additional counts
##                current_suffix = states[i][2:]
##                if unif:
##                    trans = np.array([1]*16)
##                else:
##                    trans = np.random.poisson(lam=4.0, size=16)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] = tran_probs[i,j] + trans[t]
##                    t += 1 
##                ## stay in place: add additional probability to staying in place (move = 0)
##                current_suffix = states[i]
##                if unif:
##                    trans = np.array([3])
##                else:
##                    trans = np.random.poisson(lam=20.0, size=1)
##                tran_probs[i,i] = tran_probs[i,i] + trans
##            
##            ## normalize all counts by sum to create probs that sum to 1
##            tran_probs[i,:] = tran_probs[i,:]/sum(tran_probs[i,:])
##
##    return tran_probs


    
##def generate_random_kmer_transition_probs(states, unif=False):
##    ## if allow_gaps = False, assumes each k-mer has another kmer overlapped by k-1
##    ## can set nonzero.trans to any vector -- for DNA length 4
##    k = len(states[0])
##    if k == 1:
##        pass
##    elif k == 2:
##        pass
##    else:
##        num_states = len(states)
##        tran_probs = np.zeros([num_states,num_states])
##
##        # make prefix-suffix dict -- overlaps of k-1 and k-2
##        prefix = defaultdict(list)
##        for i in range(num_states):
##            pref = states[i][:k-1]
##            prefix[pref].append(i)
##            pref = states[i][:k-2]
##            prefix[pref].append(i)
##            pref = states[i][:k-3]
##            prefix[pref].append(i)
##            pref = states[i][:k-4]
##            prefix[pref].append(i)
##
##        ## create transition probs -- can soft code the sampling parameters later if want
##        for i in range(num_states):
##            ## overlap by k-1 (move = 1)
##            current_suffix = states[i][1:k]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 365
##            else:
##                trans = np.random.poisson(lam=365.0,size=4)
##                t = 0 
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += trans[t]
##                    t += 1
##
##            ## overlap by k-2 (move = 2) -- add additional counts
##            current_suffix = states[i][2:]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 1
##            else:
##                trans = np.random.poisson(lam=4.0, size=16)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += tran_probs[i,j] + trans[t]
##                    t += 1
##
##            ## overlap by k-3 (move = 3)
##            current_suffix = states[i][3:]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 0.5
##            else:
##                trans = np.random.poisson(lam=2.0, size=64)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += tran_probs[i,j] + trans[t]
##                    t += 1
##
##            ## overlap by k-4 (move = 3)
##            current_suffix = states[i][4:]
##            if unif:
##                tran_probs[i,prefix[current_suffix]] += 0.25
##            else:
##                trans = np.random.poisson(lam=4.0, size=256)
##                t = 0
##                for j in prefix[current_suffix]:
##                    tran_probs[i,j] += tran_probs[i,j] + trans[t]
##                    t += 1
##
##            ## no overlap (move = 5)
##            tran_probs[i] += 0.1
##            
##            ## stay in place: add additional probability to staying in place (move = 0)
##            current_suffix = states[i]
##            if unif:
##                tran_probs[i,i] += 3
##            else:
##                tran_probs[i,i] = tran_probs[i,i] + np.random.poisson(lam=20.0, size=1)
##            
##            ## normalize all counts by sum to create probs that sum to 1
##            tran_probs[i,:] = tran_probs[i,:]/sum(tran_probs[i,:])
##
##    return tran_probs











##def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
##    ## states are some type of kmer
##    ## statepath is vector of numbers (indexes)
##    path_length = len(statepath)
##    moves = [0]*path_length ## first move is 0
##    k = len(states[0])
##    end = k-1
##    if k == 1 or k == 2:
##        return "This currently only works with 3-mers as smallest kmer."
##    else:
##        #init
##        seq = states[statepath[0]]
##        moves[0] = 0
##        #iter
##        for i in range(1,path_length):
##            lastSuffix = states[statepath[i-1]][1:]
##            currentPrefix = states[statepath[i]][:k-1]
##            if lastSuffix == currentPrefix:
##                seq += states[statepath[i]][end]
##                moves[i] = 1
##            else:
##                lastSuffix = states[statepath[i-1]][2:]
##                currentPrefix = states[statepath[i]][:k-2]
##                if lastSuffix == currentPrefix:
##                    seq += states[statepath[i]][end-1:]
##                    moves[i] = 2
##                elif statepath[i-1] == statepath[i]:
##                    ## by checking same state last, only heteropolymers affected
##                    ## homopolymers would be caught in first condition
##                    moves[i] = 0
##                    ## nothing is added to sequence
##                    ## could make another fxn that just spits out events and states line by line like 'template events' in f5
##                ## ELSE::: do what? ... in other one just added centroid seq regardless...
##                elif posterior_decoded:
##                    seq += states[statepath[i]][end]
##                    moves[i] = -1 
##                    ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
##                    ## it turns out that adding the base from the illegal move does not hurt the seq overall much
##    return seq, moves





### ALLOW higher gaps 3, 4, 5
def get_sequence_withgaps(states, statepath, checkoverlap=True, posterior_decoded=False):
    ## states are some type of kmer
    ## statepath is vector of numbers (indexes)
    path_length = len(statepath)
    moves = [0]*path_length ## first move is 0
    k = len(states[0])
    end = k-1
    if k == 1 or k == 2:
        return "This currently only works with 3-mers as smallest kmer."
    else:
        #init
        seq = states[statepath[0]]
        moves[0] = 0
        #iter
        for i in range(1,path_length):
            lastSuffix = states[statepath[i-1]][1:]
            currentPrefix = states[statepath[i]][:k-1]
            if lastSuffix == currentPrefix:
                seq += states[statepath[i]][end]
                moves[i] = 1
            elif statepath[i-1] == statepath[i]:
                ## by checking same state last, only heteropolymers affected
                ## homopolymers would be caught in first condition
                moves[i] = 0
                ## nothing is added to sequence
                ## could make another fxn that just spits out events and states line by line like 'template events' in f5
     
            else:
                lastSuffix = states[statepath[i-1]][2:]
                currentPrefix = states[statepath[i]][:k-2]
                if lastSuffix == currentPrefix:
                    seq += states[statepath[i]][end-1:]
                    moves[i] = 2
                else:
                    lastSuffix = states[statepath[i-1]][3:]
                    currentPrefix = states[statepath[i]][:k-3]
                    if lastSuffix == currentPrefix:
                        seq += states[statepath[i]][end-2:]
                        moves[i] = 3
                    else:
                        lastSuffix = states[statepath[i-1]][4:]
                        currentPrefix = states[statepath[i]][:k-4]
                        if lastSuffix == currentPrefix:
                            seq += states[statepath[i]][end-3:]
                            moves[i] = 4
                        else:
                            ## skip 5
                            seq += states[statepath[i]][end-4:]
                            moves[i] = 5
                   ## ELSE::: do what? ... in other one just added centroid seq regardless...
##                elif posterior_decoded:
##                    seq += states[statepath[i]][end]
##                    moves[i] = -1 
##                    ## -1 means it was an "illegal" move (move to a kmer that does not overlap by k-1 or k-2)
                    ## it turns out that adding the base from the illegal move does not hurt the seq overall much
    return seq, moves






## TESTS AND SIMULATIONS DURING EARLY DEV -- CAN PROB DELETE

##def test_viterbi(states=STATES_TRIMERS, length=10):
##    emis=generate_kmer_emission_probs(states)
##    tran=generate_random_kmer_transition_probs(states)
##    init=generate_kmer_initial_probs(states)
##    sp=generate_statepath(tran,init,states,length=length)
##    em=generate_emissions_from_statepath(emis,sp)
##    return viterbi(emis,tran,init,states,em)
##
##def simulate(states=STATES_FIVEMERS, length=10):
##    emis=generate_kmer_emission_probs(states)
##    tran=generate_random_kmer_transition_probs(states)
##    init=generate_kmer_initial_probs(states)
##    sp=generate_statepath(tran,init,states,length=length)
##    em=generate_emissions_from_statepath(emis,sp)
##    print "forward..."
##    start=time.time()
##    f,fs=forward(emis,tran,init,states,em)
##    end=time.time()
##    f1=end-start
##    print "...operation took ", f1, " seconds...."
##    print "backward..."
##    start=time.time()
##    b,bs=backward(emis,tran,init,states,em)
##    end=time.time()
##    b1=end-start
##    print "...operation took ", b1, " seconds...."
##    print "post..."
##    start=time.time()
##    postpath=posterior_decoding(f,fs,b,bs,states)
##    end=time.time()
##    print "...operation took ", end-start, " seconds...."
##    print "viterbi..."
##    start=time.time()
##    vpath,vprob=viterbi(emis,tran,init,states,em)
##    end=time.time()
##    v1=end-start
##    print "...operation took ", v1, " seconds...."
##    print ("").join([str(e) for e in ["...viterbi is ", v1/f1, "x and ", v1/b1, "x slower than F and B respectively"]])
##    print "posterior path vs known:", compare_statepath(sp,postpath)
##    print "viterbi path vs known:", compare_statepath(sp,vpath)
##    print "posterior vs viterbi:", compare_statepath(postpath,vpath)
##
##
##def simulate_delete_me(states=STATES_TRIMERS, length=10):
##    emis=generate_kmer_emission_probs(states)
##    tran=generate_random_kmer_transition_probs(states)
##    init=generate_kmer_initial_probs(states)
##    sp=generate_statepath(tran,init,states,length=length)
##    em=generate_emissions_from_statepath(emis,sp)
##    print "forward..."
##    start=time.time()
##    f,fs=forward(emis,tran,init,states,em)
##    end=time.time()
##    f1=end-start
##    print "...operation took ", f1, " seconds...."
##    print "backward..."
##    start=time.time()
##    b,bs=backward(emis,tran,init,states,em)
##    end=time.time()
##    b1=end-start
##    print "...operation took ", b1, " seconds...."
##    print "post..."
##    start=time.time()
##    postpath=posterior_decoding(f,fs,b,bs,states)
##    end=time.time()
##    print "...operation took ", end-start, " seconds...."
##    print "viterbi..."
##    start=time.time()
##    vpath,vprob=viterbi(emis,tran,init,states,em)
##    end=time.time()
##    v1=end-start
##    print "...operation took ", v1, " seconds...."
##    print "viterbi_fast..."
##    start=time.time()
##    v2path,v2pr=viterbi_fast(emis,tran,init,states,em)
##    end=time.time()
##    v2=end-start
##    print "...operation took ", v2, " seconds...."
##    print "...new viterbi ", v1/v2, "x faster than old one..."
##    print "...new viterbi is ", v2/f1, " x and ", v2/b1, "x slower than F and B respectively"
##    print "posterior path vs known:", compare_statepath(sp,postpath)
##    print "viterbi path vs known:", compare_statepath(sp,vpath)
##    print "viterbi_fast path vs known:", compare_statepath(sp,v2path)
##    print "posterior vs viterbi:", compare_statepath(postpath,vpath)
##    print "viterbi path vs viterbi fast path:", compare_statepath(vpath,v2path)
##    print "posterior vs viterbi fast:", compare_statepath(postpath,v2path)
##


if __name__ == "__main__":
    hmm = HMM()
    print "forward..."
    start = time.time()
    f, fscale = hmm.forward()
    end = time.time()
    print "...operation took ", end-start," seconds...."
    print "backward..."
    start = time.time()
    b, bscale = hmm.backward()
    end = time.time()
    print "...operation took ", end-start, " seconds...."
    print "post..."
    start = time.time()
    postpath = hmm.posterior_decoding(f, fscale, b, bscale)
    end = time.time()
    print "...operation took ", end-start, " seconds...."
    print "viterbi..."
    start = time.time()
    vpath, vprob = hmm.viterbi()
    end = time.time()
    print "...operation took ", end-start, " seconds...."
    print "posterior path vs known:", hmm.compare_statepath(postpath)
    print "viterbi path vs known:", hmm.compare_statepath(vpath)
    print "posterior vs viterbi:", hmm.compare_statepath(postpath, vpath)


