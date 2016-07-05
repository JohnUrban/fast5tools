
from fast5tools.hmm_class import *
from fast5tools.tools import *



class ParsedEvents(object):
    def __init__(self, events=None, f5=None, lead_size=50, hp_half_size=40, max_cutoff=90, max_second_cutoff=70, min_num_lg_events=3, end_trim_size=0, verbose=False):
        ## events is events dict captured from fast5class
        ## if f5class object provided, it can capture the events automatically
        assert (events is not None) or (f5 is not None)
        if events:
            self.events = events
        elif f5:
            self.events = f5.get_events_dict(readtype="input")
        self.nevents = len(self.events['mean'])

        ## find lead
        self.template_start, self.lead_state_path = lead_hmm(self.events['mean'][:lead_size])

        ## find max emission and its index (searching for HP candidate region)
        ## remove Lead Adapter on 5' of template and 3' of complement which can also have high values
        stop_looking = -20 ## cuts events of 3' end that may be high
        self.maxindex_wrt_tempstart, self.maxevent = max_and_index(self.events['mean'][self.template_start:stop_looking])

        ## Adjust maxindex relative to starting event (rather than template_start)
        self.maxindex = self.maxindex_wrt_tempstart + self.template_start

        ## Use heuristic to find a hairpin candidate region
        ## in future, grab multiple candidate regions and evaluate which has highest prob, and if that is higher than some cutoff (if so then hairpin detected, otherwise not).
        if self.maxevent > max_cutoff and sum(self.events['mean'][(self.maxindex-hp_half_size):(self.maxindex+hp_half_size)] > max_second_cutoff) >= min_num_lg_events:
            if verbose:
                sys.stderr.write("HP found\n")
            self.hairpin_detected = True
                
            ## find the starting and ending coordinates of the hairpin in this selected region
            self.hpstart_in_region, self.hpend_in_region, self.hp_state_path = hp_hmm(self.events['mean'][(self.maxindex-hp_half_size):(self.maxindex+hp_half_size)])
            
            ## get hp coord adjusted to starting event
            self.hpstart =  self.hpstart_in_region + (self.maxindex-hp_half_size)
            self.hpend = self.hpend_in_region + (self.maxindex-hp_half_size)

            ## end of template is simply hpstart (not included when used for indexing)
            self.template_end = self.hpstart

            ## beginning of complement is simply hpend+1 
            self.complement_start = self.hpend+1

            ## for now, the lead adapter at 3' end of complement is only trimmed (or not if 0) - in future, there will be comp_lead_hmm to determine if it is present and what to trim
            self.complement_end = -end_trim_size
        
        else: ## No HP candidate regions found
            self.hairpin_detected = False
            if verbose:
                sys.stderr.write("HP not found\n")
            ## template ends at provided trim size (can be 0)
            self.template_end = len(self.events['mean'])-end_trim_size
            ## complement does not exist
            self.complement_start = None
            self.complement_end = None
        

    def get_template_means(self):
        return self.events['mean'][self.template_start:self.template_end] ##up to but not including hpstart
    
    def get_complement_events(self):
        if self.hairpin_detected:
            return self.events['mean'][self.complement_start:self.complement_end]
        else:
            return np.array([])

    def get_event_tuples(self, start, end):
        for i in xrange(start, end):
            yield self.events['mean'][i], self.events['stdev'][i], self.events['start'][i], self.events['length'][i]
                            


def lead_hmm(first_N_events):
    ## 15 states: 0:14, states 0:13 part of lead profile, state14 is end/template state
    emit_probs = np.zeros([2,15]) 
    ## means
    emit_probs[0,:] = [43.93368, 51.82074, 66.3531, 76.30256, 84.15992, 89.97542, 96.22626, 100.97302, 107.33552, 100.54961, 75.71837, 46.63833, 57.33411, 43.53527, 60.0]
    ## stdevs
    emit_probs[1,:] = [2.097209, 3.526526, 2.809502, 1.954605, 1.857928, 1.793586, 1.163202, 1.120078, 2.364349, 2.866541, 13.945599, 1.991525, 16.866727, 2.678975, 5.0]
    ## initial probs - can start anywhere in profile, but mostly first 3 states
    init_probs = np.array([0.4,0.3,0.2,0,0,0,0,0,0,0,0,0,0,0,0])+0.001
    init_probs = init_probs/sum(init_probs)
    ## trans probs -- mostly trans to next state, but can skip states, also somewhat likely to stay in same state
    tran_probs = np.zeros([15,15])
    tran_probs[14,14] = 1.0
    for i in range(14): tran_probs[i,i] = 0.3
    for i in range(13): tran_probs[i,i+1] = 0.35
    for i in range(12): tran_probs[i,i+2] = 0.2
    for i in range(11): tran_probs[i,i+3] = 0.1
    for i in range(10): tran_probs[i,i+4] = 0.001
    for i in range(9): tran_probs[i,i+5] = 0.001
    for i in range(8): tran_probs[i,i+6] = 0.001
    for i in range(7): tran_probs[i,i+7] = 0.001
    ## for now only last 3 states transition to end state
    tran_probs[11,14] = 0.05
    tran_probs[12,14] = 0.1
    tran_probs[13,14] = 0.2
    ## normalize all rows to 1
    for i in range(14): tran_probs[i,:] = tran_probs[i,:]/sum(tran_probs[i,:])

    ## get viterbi path for lead adapter coordinates
    hmm = HMM(states=range(len(init_probs)), emissions=first_N_events)
    hmm.add_initial_probs(init_probs)
    hmm.add_transition_probs(tran_probs)
    hmm.add_emission_probs(emit_probs)
    vpath, vprob = hmm.viterbi()
    template_start = 0
    try:
        while vpath[template_start] != 14:
            template_start += 1
    except IndexError: ## if profile HMM does not find template start in 1st 50, then assume start is at 50
        template_start = len(first_N_events)
    return template_start, vpath
    

def hp_hmm(events,trim=0):
    ## state B,1,2,3,4,5,E = 7 states
    emit_probs = np.zeros([2,7])
    emit_probs[0,] = [65.0, 93.78638, 117.49618, 100.67429, 60.19801, 46.50402, 65.0]
    emit_probs[1,] = [6.0, 6.787453, 8.665963, 4.354063, 6.305904, 1.931336, 6.0]
    init_probs = np.array([0.7,0.2,0.1,0,0,0,0])
    init_probs = init_probs/sum(init_probs)
    tran_probs = np.zeros([7,7])
    tran_probs[6,6] = 1.0
    for i in range(7): tran_probs[i,i] = 0.3
    for i in range(6): tran_probs[i,i+1] = 0.35
    for i in range(5): tran_probs[i,i+2] = 0.2
    for i in range(4): tran_probs[i,i+3] = 0.1
    for i in range(3): tran_probs[i,i+4] = 0.001
    for i in range(2): tran_probs[i,i+5] = 0.001
    tran_probs[3,6] = 0.05
    tran_probs[4,6] = 0.1
    tran_probs[4,5] = 0.7 ## state 4 usually goes directly to 5 (occasionally 2 events, but have not seen more -- all other states tend to stay in-state longer)
    for i in range(7): tran_probs[i,] = tran_probs[i,:]/sum(tran_probs[i,:])

    ## get viterbi path for hairpin adapter coordinates
    hmm = HMM(states=range(len(init_probs)), emissions=events)
    hmm.add_initial_probs(init_probs)
    hmm.add_transition_probs(tran_probs)
    hmm.add_emission_probs(emit_probs)
    vpath, vprob = hmm.viterbi()
    hpstart = 0
    while vpath[hpstart] < 1:
        hpstart += 1
    hpend = len(vpath)-1
    while vpath[hpend] > 5:
        hpend -= 1
    return hpstart-trim, hpend+trim, vpath


