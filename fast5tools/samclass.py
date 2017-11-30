import os, sys, re, pybedtools, pandas
from collections import defaultdict
import numpy as np

class Sam(object):
    ''' Reproduce Same file with:
            sam = Sam(args.sam)
            for record in sam:
                print record
    '''
    def __init__(self, samfile):
        self.sam = self.open(samfile) ## also establishes self.filename, self.filebasename, self.abspath
        self.record = None
        self.header_parsed = False
        self.first_aln_used = False
        self.n_records = 0
        self.HD = {}
        self.SQ = {}
        self.PG = {}
        self.header_string = ''
        self.parse_header()
        
        
    def open(self, samfile):
        if samfile in ("-","stdin"):
            self.filename = "stdin"
            self.filebasename = "stdin"
            self.abspath = "stdin"
            return sys.stdin
        else:
            self.filename = samfile
            self.filebasename = (".").join(samfile.split("/")[-1].split(".")[:-1]) ## takes entire basename except ".sam"
            self.abspath = os.path.abspath(samfile)
            return open(samfile, 'r')


    def updateHeader(self, target, headerline, i):
        ans = headerline[i].split(":")  
        k = ans[0]
        v = ':'.join(ans[1:])
        try:
            target[k].append(v)
        except:
            target[k] = [v]
        
        
    def updateHD(self,headerline, i):
        self.updateHeader(self.HD, headerline, i)

    def updateSQ(self,headerline, i):
        self.updateHeader(self.SQ, headerline, i)


    def updatePG(self,headerline, i):
        self.updateHeader(self.PG, headerline, i)
            
    def parse_header(self):
        if not self.header_parsed: #self.alignment is None:
            headerline = self.sam.readline()
            while headerline[0] == '@':
                self.header_string += headerline
                headerline = headerline.strip().split('\t')
                tag = headerline[0][1:3]
                if tag == 'HD':
                    self.updateHD(headerline, 1)
                    self.updateHD(headerline, 2)
                elif tag == 'SQ':
                    self.updateSQ(headerline, 1)
                    self.updateSQ(headerline, 2)
                elif tag == 'PG':
                    self.updatePG(headerline, 1)
                    self.updatePG(headerline, 2)
                    self.updatePG(headerline, 3)
                    self.updatePG(headerline, 4)
                headerline = self.sam.readline()
            self.record = SamRecord(headerline)
            self.n_records += 1 ## assumes there is an alignment....
            self.header_parsed = True
            self.header_string = self.header_string.strip() ## takes last \n off

    def get_header_string(self):
        return self.header_string
                    
            
    def get_next_sam_record(self):
        '''Obtains next alignment (or first alignment if immediately after parsing header).'''
        if not self.header_parsed: #self.alignment is None:
            self.parse_header()
        elif self.header_parsed and not self.first_aln_used:
            self.first_aln_used = True ## NOT SURE THIS APPLIES ANYMORE... may be able to take condition out
        elif self.header_parsed and self.first_aln_used:
            self.record = SamRecord( self.sam.readline() )
            self.n_records += 1
        return self.record


    def __iter__(self):
        return self

    def next(self):
        try:
            self.get_next_sam_record()
            return self.record
        except Exception as e:
            raise StopIteration





class SamRecord(object):
    def __init__(self, sam_record):
        self.cigar = None
        self.parsed_aln = {} 
        self.seq_len = None 
        self.read_len = None 
        self.edit_dist = None 
        self.clipping_dist = None
        self.fast5info = None
        self.alignment = sam_record
        self.parse_alignment()

    def __str__(self):
        ## Return string_from_parsed b/c it updates may have been made
        return self.string_from_parsed_alignment()

    def parse_alignment(self):
        aln = self.alignment.strip().split('\t')
        self.parsed_aln['QNAME'] = aln[0]
        self.parsed_aln['FLAG'] = aln[1]
        self.parsed_aln['RNAME'] = aln[2]
        self.parsed_aln['POS'] = int(aln[3])
        self.parsed_aln['MAPQ'] = int(aln[4])
        self.parsed_aln['CIGAR'] = aln[5]
        self.parsed_aln['RNEXT'] = aln[6]
        self.parsed_aln['PNEXT'] = aln[7]
        self.parsed_aln['TLEN'] = aln[8]
        self.parsed_aln['SEQ'] = aln[9]
        self.parsed_aln['QUAL'] = aln[10]
        self.parsed_aln['EXTRA'] = ('\t').join( aln[11:] )

    def is_aligned(self):
        ## TODO - better definition?
        ''' This definition is useful for single-end/long reads - but not for PE reads.'''
        return self.get_rname_field() != '*'

    def is_unaligned(self):
        ## TODO - better definition?
        ''' This definition is useful for single-end/long reads - but not for PE reads.'''
        return self.get_rname_field() == '*'


    def get_sam_fields(self):
        return ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'EXTRA']

    def get_qname_field(self):
        return self.parsed_aln['QNAME']

    def get_flag_field(self):
        return self.parsed_aln['FLAG']

    def get_rname_field(self):
        return self.parsed_aln['RNAME']

    def get_pos_field(self):
        return self.parsed_aln['POS']

    def get_mapq_field(self):
        return self.parsed_aln['MAPQ']

    def get_cigar_field(self):
        return self.parsed_aln['CIGAR']

    def get_rnext_field(self):
        return self.parsed_aln['RNEXT']

    def get_pnext_field(self):
        return self.parsed_aln['PNEXT']

    def get_tlen_field(self):
        return self.parsed_aln['TLEN']

    def get_seq_field(self):
        return self.parsed_aln['SEQ']

    def get_qual_field(self):
        return self.parsed_aln['QUAL']

    def get_extra_fields(self):
        return self.parsed_aln['EXTRA']
    
    def string_from_parsed_alignment(self):
        return ('\t').join( [ str(self.parsed_aln[k]) for k in self.get_sam_fields() ] )


    def get_cigar_counts(self):
        # if present as single '*', then CIGAR is unavailable
        # cigar chars are in 'MIDNSHP=X'
        # some ways to get seq and read lengths from CIGAR if not storing cigar dict
        #   self.seq_len = sum([int(e) for e in re.findall('(\d+)[MIS=X]', self.parsed_aln['CIGAR'])])
        #   self.read_len = sum([int(e) for e in re.findall('(\d+)[HMIS=X]', self.parsed_aln['CIGAR'])])
        self.cigar = {"M":0, "I":0, "D":0, "N":0, "S":0, "H":0, "P":0, "=":0, "X":0}
        for count,char in re.findall('(\d+)([MIDNSHP=X])', self.parsed_aln['CIGAR']):
            self.cigar[char] += int(count)
        return self.cigar

    def get_SEQ_len(self):
        '''SEQ as in the 10th SAM field -- or what should be there as in case of secondary alignments that often leave it as "*"'''
        if self.cigar is None: 
            self.get_cigar_counts()
        if self.seq_len is None:
            self.seq_len = sum([self.cigar[e] for e in 'MIS=X']) 
        return self.seq_len

    def get_read_len(self):
        if self.cigar is None: 
            self.get_cigar_counts()
        if self.read_len is None:
            self.read_len = sum([self.cigar[e] for e in 'HMIS=X']) 
        return self.read_len

    def get_SEQ_len_without_clipped_regions(self):
        ''' Length of actually aligned portion.'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return sum([self.cigar[e] for e in 'MI=X']) 


    def update_pos_field(self, add=0, replace=None): ## example use: doubled ecoli genome for mapping over breakpoint; now want to subtract G_size from anything greater than G_size to put back on first copy
        assert add == 0 or replace is None
        if add != 0:
            self.parsed_aln['POS'] += add
        elif replace is not None:
            self.parsed_aln['POS'] = replace


    def find_given_extra_field(self, regex, output_type=None):
        '''Returns only first match. There should only be one match.'''
        for e in self.parsed_aln['EXTRA'].split('\t'):
            found = re.search(regex, e)
            if found is not None:
                break
        try:
            if output_type is None:
                return found.group(1)
            else:
                return output_type(found.group(1))
        except:
            return None
            


    def get_edit_dist_field(self):
        '''Assumes NM:i: field is present for now...'''
        if self.edit_dist is None:
            self.edit_dist = self.find_given_extra_field(regex='NM:i:(\d+)', output_type=int)
        return self.edit_dist


    def get_fast5_field(self):
        '''Only works with SAM files that have the F5:Z: tag from fast5tools'''
        if self.fast5info is None: 
            self.fast5info = self.find_given_extra_field(regex='F5:Z:(.*)', output_type=str)
        return self.fast5info
            
    def get_AS_field(self):
        '''Assumes AS:i: field is present for now...'''
        return self.find_given_extra_field(regex='AS:i:(.*)', output_type=int)

    def get_XS_field(self):
        '''Assumes XS:i: field is present for now...'''
        return self.find_given_extra_field(regex='XS:i:(.*)', output_type=int)

    def get_MD_field(self):
        '''Assumes NM:i: field is present for now...'''
        return self.find_given_extra_field(regex='MD:Z:(.*)', output_type=str)

    def get_XA_field(self):
        '''Assumes NM:i: field is present for now...'''
        return self.find_given_extra_field(regex='XA:Z:(.*)', output_type=str)


    def get_clipping_dist(self):
        if self.cigar is None: 
            self.get_cigar_counts()
        if self.clipping_dist is None: 
            self.clipping_dist = sum([self.cigar[e] for e in 'HS'])
        return self.clipping_dist
            

    def get_cigar_list(self):
        return [(int(count),char) for count,char in re.findall('(\d+)([MIDNSHP=X])', self.get_cigar_field())]

##    def convert_cigar_element_

    def get_edit_dist_with_clipping(self):
        if self.is_aligned():
            return int(self.get_edit_dist_field()) + int(self.get_clipping_dist())
        else:
            return None

    def get_reference_aln_len(self):
        ''' This gives length of alignment from POV of reference according to pos field and cigar.
            Add up matches to ref and deletions from ref to get ref len in aln.
            The reference end coord is Pos + MD - 1'''
        return sum([self.cigar[e] for e in 'MD=X']) ## used to be 'MD', which I validated -- however, I then reasoned '=X' should be added for aligners that use those instead of M...

    def get_reference_end_pos(self):
        ''' This gives tuple of 1-based (as SAM is) start/end on reference according to pos field and cigar.
            The reference end coord is Pos + MD - 1 (i.e. add up matches to ref and deletions to ref to get ref len in aln)'''
        ## Note: obtain "start pos" with self.get_pos_field()
        return self.get_pos_field() + self.get_reference_aln_len() - 1


    def get_clipping_length(self, cigaridx1, cigaridx2, clip=0):
        ''' Do not necessarily see use for this beyond 5' and 3' clip lengths - but left as a fxn to simplify'''
        ''' 5' clip len: cigaridx1=0, cigaridx2=1'''
        ''' 3' clip len: cigaridx1=-1, cigaridx2=-2'''
        if self.get_clipping_dist() > 0:
            cigar_list = self.get_cigar_list()
            if cigar_list[cigaridx1][1] in ('H', 'S'):
                clip += cigar_list[cigaridx1][0]
                if cigar_list[cigaridx2][1] == 'S': #possible if H was first
                    clip += cigar_list[cigaridx2][0]
        return clip

    def get_5prime_clip_len(self):
        return self.get_clipping_length(cigaridx1=0, cigaridx2=1, clip=0)

    def get_3prime_clip_len(self):
        return self.get_clipping_length(cigaridx1=-1, cigaridx2=-2, clip=0)

    def estimate_reference_to_read_length_ratio(self):
        return self.get_reference_aln_len() / float( self.get_SEQ_len_without_clipped_regions() )

    def estimate_reference_length_under_5prime_clip(self):
        ## could us SCALE = sum(MD=X)/sum(MI=X) to get scale factor for clip_lens
        ## Do int(round(Clip_Len * SCALE)) <---self.estimate_reference_to_read_length_ratio()
        ## This will shrink clip_len if reference part is smaller -- make it bigger if it is bigger
        ## Pitfall: as this calculation would estimate the scale_factor potentially from small portions of a long read
        ##          it might not reflect the scale_factor that would be obtained over all aligned portions of read
        return self.get_5prime_clip_len() * self.estimate_reference_to_read_length_ratio()


    def estimate_reference_length_under_3prime_clip(self):
        ## see NOTES in estimate_reference_len_under_5prime_clip
        return self.get_3prime_clip_len() * self.estimate_reference_to_read_length_ratio()

    def estimate_reference_length_for_full_read_alignment(self, use_real_clip_sizes=False):
        if use_real_clip_sizes:
            return self.get_5prime_clip_len() + self.get_reference_aln_len() + self.get_3prime_clip_len()
        else:
            return self.estimate_reference_length_under_5prime_clip() + self.get_reference_aln_len() + self.estimate_reference_length_under_3prime_clip()

    def get_adjusted_pos(self, initialpos, clip=0, extra=0, direction=0):
        '''This used to be a more involved function, but has been broken into smaller pieces that basically
            render it as a simple arithmetic fxn.
            initialpos is just some integer (typically the start POS or end POS of alignment
            clip is just some integer (typically the 5' or 3' clipping len)
            extra is just some integer (typically a fudge factor such as 10% of read length)
            direction is just 1 or -1 (though can be used as a scale factor).
                if direction=1, clip and extra will be added to initial pos
                    (if abs(direction) > 1, then the scaled factor will be added)
                if direction=-1, clip and extra will be subtracted from initial pos
                    (if abs(direction) > 1, then the scaled factor will be subtracted)
                if direction=0 (default), initialpos is returned as is
            '''
        return initialpos + direction*clip + direction*extra

    def get_clipping_adjusted_start(self, extra=0, use_real_clip_sizes=True):
        ''' Adjusts start pos to be the number of HS-clipped bases 5' to reported POS'''
        initialpos = self.get_pos_field()
        if use_real_clip_sizes:
            clip = self.get_5prime_clip_len()
        else:
            clip = self.estimate_reference_length_under_5prime_clip()
        direction = -1
        return self.get_adjusted_pos_field(initialpos, clip, extra, direction)

    def get_clipping_adjusted_end(self, extra=0):
        initialpos = self.get_reference_end_pos()
        if use_real_clip_sizes:
            clip = self.get_3prime_clip_len()
        else:
            clip = self.estimate_reference_length_under_3prime_clip()
        direction = 1
        return self.get_adjusted_pos(initialpos, clip, extra, direction)
        

    def get_proportion_of_read_len(self, proportion=0.1):
        '''Returns int'''
        return int( self.get_read_len() * proportion )

    def get_genomic_window_around_alignment(self, flank=0, adjust_for_clipping=True, identifier=False):
        '''Get genomic window surrounding an alignment with some buffer/flanks -- e.g. for local re-alignment.'''
        ''' Default: W/o buffer/flanking sequence, it gives the clipping-adjusted guesses for start and end positions.'''
        ''' Add buffer/flanks two ways:'''
        '''     (1) int > 1 adds/subtracts that int.
                (2) float [0,1] adds/subtracts that proportion of read length
                    NOTE: 1.0 gives proportion while 1 gives 1 bp.'''
        '''1-based, closed.'''
        '''In splitread class - can check for overlap among various genomic windows from split alns of same read.'''
        assert flank >= 0 ## no negative numbers are meaningful here
        if flank == 0:
            extra = 0
        elif flank > 1 or (flank == 1 and isinstance(flank, int)):
            extra = int(flank) 
        elif flank >= 0 and flank < 1 or (flank == 1 and isintance(flank, float)):
            extra = int(self.get_proportion_of_read_len(proportion))
        chrom = self.get_record(i).get_rname_field()
        if adjust_for_clipping:
            start = self.get_clipping_adjusted_start(extra=extra)
            end = self.get_clipping_adjusted_end(extra=extra)
        else:
            start = self.get_pos_field() - extra
            end = self.get_reference_end_pos() + extra
        if identifier:
            return (chrom, start, end, identifier)
        else:
            return (chrom, start, end)










### SPLIT ALN ANALYSES
### E.g. BWA can align long reads across multiple records
### The SamSplitAlnAggregator class inherits most of the Sam class
###     However, it takes in SAMs sorted by name
###     And yields >= 1 Sam Record in a list corresponding to all records with same QNAME
### The SplitReadSamRecord class takes in the list of SamRecord objects to handle desired/sterotyped analyses on those records

class SamSplitAlnAggregator(Sam):
    '''Assumes SAM file is sorted by name'''
    def __init__(self, samfile):
        self.endfile = False
        self.records = []
        self.n_splits = 0
        Sam.__init__(self, samfile)

    def parse_header(self):
        if not self.header_parsed: #self.alignment is None:
            headerline = self.sam.readline()
            while headerline[0] == '@':
                self.header_string += headerline
                headerline = headerline.strip().split('\t')
                tag = headerline[0][1:3]
                if tag == 'HD':
                    self.updateHD(headerline, 1)
                    self.updateHD(headerline, 2)
                elif tag == 'SQ':
                    self.updateSQ(headerline, 1)
                    self.updateSQ(headerline, 2)
                elif tag == 'PG':
                    self.updatePG(headerline, 1)
                    self.updatePG(headerline, 2)
                    self.updatePG(headerline, 3)
                    self.updatePG(headerline, 4)
                headerline = self.sam.readline()
            self.record = SamRecord(headerline)
            self.current_record_name = self.record.get_qname_field()
            self.n_records += 1 ## assumes there is an alignment....
            self.header_parsed = True
            self.header_string = self.header_string.strip() ## takes last \n off

    def get_all_records_with_shared_named(self):
        ''' Starting with list of records given (empty list by default),
            continue to get next record until name doesn't match'''
        while self.record.get_qname_field() == self.current_record_name:
            self.records.append( self.record )
            try:
                self.record = SamRecord( self.sam.readline() )
                self.n_records += 1
            except: ## no more records to add?
                self.record = None
                self.endfile = True
                break
        if not self.endfile:
            self.current_record_name = self.record.get_qname_field()
##        return self.records
        return SplitReadSamRecord( self.records )

    def get_next_set_of_sam_records(self):
        '''Obtains next alignment (or first alignment if immediately after parsing header).'''
        if not self.header_parsed: ## first record is established in parsing header (done by default so should always be present)
            self.parse_header()
        else:
            self.records = [] ## reset records
        return self.get_all_records_with_shared_named()

    def next(self):
        try:
            self.get_next_set_of_sam_records()
            return self.records
        except Exception as e:
            raise StopIteration






class SplitReadSamRecord(object):
    '''Current objectives for this class:
        Obj1: Use [samrecords] to identify most likely genomic region a read came from.
            If one alignment, then that region.
            If two or more alignments, then some type of majority rule.
        Obj2: Use all parts of split alignment to get better idea of pct identity.
            - pct identity of longest aln
    '''
    def __init__(self, samrecords):
        self.records = samrecords
        self.num_aln = None

    def get_num_aln(self):
        if self.num_aln is None:
            self.num_aln = len(self.records)
        return self.num_aln

    def get_record(self,index):
        return self.records[index]
    
    def genomic_window_ovlps(self):
        ## Use pybedtools to find overlaps
        pass


    def ovlp(a,b,d=0):
        ''' a and b are 3-tuples of (chr,start,end)'''
        ''' assumes sorted such that, if ovlp, a is before b along number line.'''
        ''' d is distance - allows a gap up to d between intervals to still be an overlap - default 0'''
	if a[0] == b[0] and a[1] <= b[1] and b[1] <= a[2]:
            return True
	return False

    def merge(l,d=0):
        ''' Takes in list of 3-tuples of (chr.start,end)'''
        ''' Returns updated list where overlapping intervals are merged'''
        ''' d is distance - allows a gap up to d between intervals to still be an overlap - default 0'''
        a = l[0]
        update = []
        for i in range(1,len(l)):
            b = l[i]
            if self.ovlp(a,b,d):
                a = merge_tuples(a,b) ### new merge_tuples is identifier-aware; old = (a[0],a[1],b[2])
            else:
                update.append(a)
                a = b
        ## at end, if last thng ovlpd, then it merged interval in the making (a) needs to added to 'update'
        ##         if last thing didn't ovlp, then final interval, named 'a' before exiting, needs to be added to 'update'
        update.append(a)
        return update

    def merge_tuples(a,b):
        if len(a) == 3:
            return (a[0],a[1],b[2])
        elif len(a) == 4: ##has identifier
            identifier = str(a[3]) + "," + str(b[3])
            return (a[0], a[1], b[2], identifier)


    def determine_window_length(self, genomic_window):
        ''' Meant to take in single genomic_window tuple and report length back'''
        return genomic_window[2] - genomic_window[1] + 1

    def determine_window_lengths(self, genomic_windows):
        ''' Meant to take in list of genomic_window tuples, and return list with window lengths.'''
        return np.array([self.determine_window_length( gw ) for gw in genomic_windows])


    def determine_longest_window(self, genomic_windows):
        ''' Meant to take list of genomic_window tuples and report the indexes of longest window(s) in list.
        More than one index is returned only when there is a tie.'''
        a = self.determine_window_lengths(genomic_windows)
        return np.array(range(len(a)))[a == a.max()]

    def determine_window_proportions(self, genomic_windows, scale_factor=False):
        ''' Default scaling is to summed window length - giving proportion each window makes up of all windows..
            Example of another scale_factor one might want to provide, is read length to determine if some window
            makes up a desired proportion of the read length...'''
        a = self.determine_window_lengths(genomic_windows)
        if not scale_factor: ## Then scale_factor is the sum of window lengths (default)
            scale_factor = sum(a)
        return a.astype(np.float)/scale_factor

    def determine_majority_window(self, genomic_windows, majority=0.5, scale_factor=False):
        ''' Meant to look at list of genomic_window tuples and report whether one window
            makes up the majority percentage of the summed length of all windows. (reports index)
            This does not mean the largest proportion (which would be same index as longest window).
            This needs the max window to be larger than a given cutoff.
            Majority is not simply 51% - it should be some number that makes sense.
            For example: 0.7 is saying I'd like at least 70% of the total genomic window to be encompassed by one.
            ...So long as majority >0.5, there should be only one gw index returned.
            ...Nonetheless, it is coded anticipating other uses that may allow more than one window index returned.
            If a single window = 0.5 and there are >2 windows, then it is majority...'''
        gw_props = self.determine_window_proportions(genomic_windows, scale_factor)
        gw_max = gw.props.max()
        if gw_max >= majority:
            return np.array(range(len(gw_props)))[gw_props == gw_max]
        else:
            return False


    def determine_longest_window_ratio(self, genomic_windows):
        ''' Meant to look at list of genomic_window tuples and report whether the ratio between
            the longest and second longest window.'''
        gw_props = self.determine_window_proportions(genomic_windows, scale_factor)
        if len(gw_props) == 1:
            return float('inf')
        else:
            gw_props.sort() #smallest to largest gw_props[-1]/gw_props[-2]
            gw_props = gw_props[-1::-1] #largest to smallest
            return gw_props[0] / float(gw_props[1])

    def alignments_ordered_like_read(self, indexes=[]):
        ''' indexes = list of indexes of the records you'd like to check.
                e.g. you might have 3 split alignments but only want to check 2 and 3.
            This will ensure that the records actually occur on same chrom (just to be safe),
                then checks to see if the 5' clipping in reads agrees with POS field ordering.
            BWA gives partitioned split alignments. So if 3 split alignments occur consecutively along
            a genomic region, then their 5' clipping lengths should each be longer than the last.'''
        if not indexes:
            indexes = range(self.get_num_aln())
        d = {'pos':[], 'clip':[]}
        for i in indexes:
            d['pos'].append( self.get_record(i).get_pos_field() )
            d['clip'].append( self.get_record(i).get_5prime_clip_len() )
        df = pandas.DataFrame(d)
        ans = df.sort('pos')['pos'] == c.sort('clip')['pos']
        return sum(ans) == len(indexes)
        
                

    def get_genomic_window(self, flank=0.1, merge_dist=0, majority=0.7):
        '''Returns 3-tuple'''
        '''flank = as in get_genomic_window_around_alignment() described below'''
        ''' merge_dist = d from self.merge(): allows a gap up to d between intervals to still be an overlap - default 0'''
        ''' From the individual get_genomic_window_around_alignment():'''
        '''Get genomic window surrounding an alignment with some buffer/flanks -- e.g. for local re-alignment.'''
        ''' Default: W/o buffer/flanking sequence, it gives the clipping-adjusted guesses for start and end positions.'''
        ''' Add buffer/flanks two ways:'''
        '''     (1) int > 1 adds/subtracts that int.
                (2) float [0,1] adds/subtracts that proportion of read length
                    NOTE: 1.0 gives proportion while 1 gives 1 bp.'''
        '''1-based, closed.'''
        '''In splitread class - can check for overlap among various genomic windows from split alns of same read.'''
        if self.get_num_aln() == 1:
            return self.get_record(0).get_genomic_window_around_alignment(flank=flank, adjust_for_clipping=True, identifier="single_alignment")
        elif self.get_num_aln() > 1:
            genomic_windows = []
##            bedstr = ''
            for i in range(self.get_num_aln):
                gw = self.get_record(i).get_genomic_window_around_alignment(flank=flank, adjust_for_clipping=True, identifier=i) 
                genomic_windows.append( gw )
##                bedstr = chrom + "\t" + str(gw[0]-1) + "\t" + str(gw[1])
            #bedtools approach
##            bedtool = pybedtools.BedTool( bedstr, from_string=True )
##            merged = bedtool.merge(d=merge_dist)
            #python-only approach
            sorted_genomic_windows = sorted(genomic_windows) ## use sorted not .sort(), so can access genomic_windows below
            gw_merge = self.merge(l=sorted_genomic_windows, d=merge_dist)
            num_after_merge = len(gw_merge)
            if num_after_merge == 1:
                #return the composite genomic window formed by overlapping genomic windows
                identifier = "single_merged_window"
                single_merged_window = (gw_merge[0], gw_merge[1], gw_merge[2], identifier)
                return single_merged_window
            elif num_after_merge > 1:
                # find majority
                # if num_merge == self.get_num_aln(): then can use SamRecord objects for more info about majority alignment
                # else: some merges occurred -- can find majority based on aln len -- can also dig into SamRecords if nec
                # can look first to see if all merges are on the same chr -- then look at start of elemnt1 and end of the last element
                #   to determine total length/genomic window.... if it is "reasonable", just give the entire window...
                #       "reasonable" might be 2X read length or something...
                #       a problem with using the bp readlen is that the base-caller might seriously shorten the reads when it is dealing with lots of analog...
                #       could do 10x readlen or min(X, 10xreadlen) .... but X needs to be > readlen + fudgefactor...
                #       could also make it proportional to the number of events or number of raw data points...
                #           Just checked - N=1 - same DNA sequence for BrdU, EdU, IdU, and T -- all produce similar base-seq lens but analogs segment into 1.3x more events than T (or ~25% more than T)
                ####
                #### Note b/c the genomic windows were extended to at least the length of reads by adjusting for clipping,
                ####    they should not be used as is for determining the majority alignment
                ####    Nor necessarily should it be a % of read length, since large portions of reads might not align
                ####        Moreover - since these reads can be highly error prone (w/ normal bases), 50% of read might not map to 50% of alignment
                ####                -  and w/ the analogs... read length might be really off... only would want to trust aligned regions...
                ####    Best might be as a % of total alignment length MD=X.
                if num_after_merge == self.get_num_aln():
                    ## nothing merged in above conditions -- meaning alignments are relatively far away from each other
                    ## is there a majority? Aln that makes up > X% of summed alignment length?
                    genomic_alignments = []
                    for i in range(self.get_num_aln):
                        ref_coords = self.get_record(i).get_genomic_window_around_alignment(flank=0, adjust_for_clipping=False, identifier=i) 
                        genomic_alignments.append( ref_coords )
                    majority = self.determine_majority_window(genomic_alignments, majority)
                    if majority:
                        ## need to get genomic window around alignment for majority alignment
                        ## these types of genomic windows were obtained above in genomic_windows -- not genomic_alignments
                        majority_alignment = genomic_windows[majority]
                        identifier = 'majority_alignment:' #+ proportion_aligned_bases + proportion_estimated_alignment_len + num_alignments + 2nfbest...
                        majority_window = (majority_alignment[0], majority_alignment[1], majority_alignment[2], identifier)
                        return majority_window
                    else:
                        ## return longest? return multiple? use alignment scores?
                        identifier = 'longest_alignment:' ## plus info....
                        pass
                else: # compare merged and single alns all at same time?
##                    if_majority here
##                    if_merged_majority
##                    longest_btwn_merged_and_single

                    pass
                    # there was some merging, meaning there is a site in the genome with more than one alignment within reasonable distance (in above conditions)
                    # check if the merged alignments make up a majority of all alignments...
                    #  ... if not -- see if a single alignment makes up a majority...
                    #       or should this logic be reversed...? look for majority single... then look for majority from merged...
                    #       Well if a single makes up 0.5, then the merge can only sum to 0.5 anyway...
                    #       Does one trust one contiguous aln more than >=2 possibly smaller disjointed alignments?
                    #           Id say ..sure...
                    #       If no single majority... if no merged majority, then what?
                    #       which is longer... the single alignment or the sum of merged alignments? go with that one...
                    #           note proportion of alignments taken by final... note proportion of readlength... or estimated ref len...
                    # Do I check len of merged_aln_1_start to merged_aln_N_end?
                    # Or do I just sum their aln lengths?
                    # if doing merged_ref_len, then shouldnt I check to see if that is longer than longest single aln in beginng
                    # i.e. majority calculation would have different denominators if using sum_aln_lens vs front_to_end_merge_len
                

                
            
            






#Functions on SAM line
def sam_seq_len(x):
    ''' Length of SEQ (part of READ that is soft-clipped and/or aligned) as defined by SAM specifications.'''
    return sum([int(e) for e in re.findall('(\d+)[MIS=X]', x)])

def sam_read_len(x):
    ''' Length of entire READ as defined by SAM specifications.'''
    return sum([int(e) for e in re.findall('(\d+)[HMIS=X]', x)])

def sam_alnseq_len(x):
    ''' Length of READ that is aligned only (no soft clipping) as defined by SAM specifications.'''
    return sum([int(e) for e in re.findall('(\d+)[MI=X]', x)])

def sam_refalnseq_len(x):
    '''Length of REFERENCE sequence underlying the aligned portion of read.'''
    return sum([int(e) for e in re.findall('(\d+)[MD=X]', x)])

#### STUFF TO WORK IN TO MAYBE THE SAM CLASS
##import re, sys
##from collections import defaultdict
##import numpy as np
##
##readlengths = defaultdict(list)
##alnlengths = defaultdict(list)
##editdists = defaultdict(list)
##alnscores = defaultdict(list)
##mapqs = defaultdict(list)
##contigs = defaultdict(set)
##UN_alnscores = defaultdict(list)
##UN_mapqs = defaultdict(list)
##names = set([])
##totalaln = 0
##totalunaln = 0
##for line in sys.stdin:
##    line = line.strip().split()
##    name = line[0]
##    names.add(name)
##    if line[1] != '4':
##        ##rlen = sum([int(e) for e in re.findall('(\d+)[MIS=X]', line[5])])
##        rlen = sum([int(e) for e in re.findall('(\d+)[MISH=X]', line[5])])
##        readlengths[name].append( rlen  )
##        totalaln += 1
##        alnlengths[name].append( sum([int(e) for e in re.findall('(\d+)[MI=X]', line[5])]) )
##        editdists[name].append( int(line[11].split(':')[2])  )
##        alnscores[name].append( int(line[13].split(':')[2])  )
##        mapqs[name].append( int(line[4])  )
##        contigs[name].add( line[2]  )
##    else:
##        readlengths[name].append( len(line[9]) ) 
##        totalunaln += 1
##        UN_alnscores[name].append( 0 )
##        UN_mapqs[name].append( 0 )
##
###initate stats vars
##summapq_old = 0
##num_0_ctg = 0
##num_1_ctg = 0
##num_multi_ctg = 0
##readlensum = 0
##alnlensum = 0
##alnscoresum = 0
##alnscoresum2 = 0
##mapqsum = 0
##mapqsum2 = 0
##editsum = 0
##editsumwithunmapped = 0
##numaln = 0
##numaln_un = 0
##numunaln = 0
##unmappedreadlensum = 0
##matchsum = 0
##matchsum2 = 0
##totaleditsum = 0
##unalnportionsum = 0
##sumalnscore_per_aln = 0
###initiate per-read outfile
##out = open('$PRE-per-read.txt','w')
##
##for name in list(names):
##    ALN = name in alnlengths
##    UN = name in UN_alnscores
##    readlen = readlengths[name][0]
##    if ALN:
##        nctg = len( contigs[name] )
##        aln_n = len( alnlengths[name] )
##        alnlen = sum( alnlengths[name] )
##        edit = sum( editdists[name] )
##        alnscore = np.mean( alnscores[name] )
##        mapq = np.mean( mapqs[name] )
##        unmappedreadlen = 0
##        if UN:
##            editun = 0
##            aln_n2 = aln_n + len( UN_alnscores[name] )
##            alnscore2 = np.mean( alnscores[name]+UN_alnscores[name] )
##            mapq2 = np.mean( mapqs[name]+UN_mapqs[name] )
##        else:
##            editun = 0
##            aln_n2 = aln_n
##            alnscore2 = alnscore
##            mapq2 = mapq
##    elif UN:
##        nctg = 0
##        alnlen = 0
##        aln_n = 0
##        edit = 0
##        editun = readlen
##        alnscore = 0
##        mapq = 0
##        aln_n2 = 0
##        alnscore2 = 0
##        mapq2 = 0
##        unmappedreadlen = readlen
##    unalnportion = readlen - alnlen
##    totaledit = edit + unalnportion
##    ##match = readlen - unmappedreadlen - edit ## error b/c it counts unaligned portions as matches
##    match = alnlen - edit
##    match2 = readlen - totaledit ## same as alnlen-edit --> readlen-totaledit = readlen-(edit+unalnportion) = readlen-(edit+readlean-alnlen) = readlen-edit-readlen+alnlen = alnlen+(readlen-readlen-edit) = alnlen-edit
##    out.write( ('\t').join([str(e) for e in [name, nctg, aln_n, readlen, alnlen, unalnportion, edit, totaledit, editun, alnscore, mapq, aln_n2, alnscore2, mapq2, match, match2]]) + '\n' )
##
##    #additional analysis
##    summapq_old += sum( mapqs[name]+UN_mapqs[name] )
##    sumalnscore_per_aln += sum( alnscores[name]  )
##    if nctg == 0: 
##        num_0_ctg += 1
##    elif nctg == 1: 
##        num_1_ctg += 1
##    elif nctg > 1: 
##        num_multi_ctg += 1
##    readlensum += readlen
##    alnlensum += alnlen
##    mapqsum += mapq
##    editsum += edit
##    editsumwithunmapped += edit + editun
##    alnscoresum += alnscore
##    mapqsum2 += mapq2
##    alnscoresum2 += alnscore2
##    unmappedreadlensum += unmappedreadlen
##    matchsum += match
##    matchsum2 += match2
##    totaleditsum += totaledit
##    unalnportionsum += unalnportion
##    if ALN and UN:
##        numaln += 1
##        numaln_un += 1
##    elif ALN:
##        numaln += 1
##    elif UN:
##        numunaln += 1
###close per-read outfile
##out.close()
##
##
###final analysis - open stats out file and write
###totalaln = 0
###totalunaln = 0
##numentries = totalaln + totalunaln
##numuniqentries = numaln + numunaln
##totalnumaln = totalaln
##totalnumuniqaln = numaln
####summapq_old = summapq_old
##avgmapq_old = summapq_old / float( totalaln )
##pctaln_old = totalnumuniqaln / float( numuniqentries )
##pctunaln_old = (numuniqentries - totalnumuniqaln) / float( numuniqentries )
##alnratio = totalnumaln / float( totalnumuniqaln )
##
##avgalnlen_allaln = float(alnlensum) / totalnumaln
##avgalnlen_allalnread = float(alnlensum) / totalnumuniqaln
##avgalnlen_allread = float(alnlensum) / numuniqentries
##
##avgmapq = mapqsum/float(numaln)
##avgmapq2 = mapqsum2/float(numaln)
####avgmapq_with_unmap = ((mapqsum+mapqsum2)/2.0)/float(numaln+numunaln)
##avgmapq_with_unmap = mapqsum/float(numaln+numunaln)
##avgalnscore = alnscoresum/float(numaln)
##avgalnscore2 = alnscoresum2/float(numaln)
####avgalnscore_with_unmap = ((alnscoresum+alnscoresum2)/2.0)/float(numaln+numunaln)
##avgalnscore_with_unmap = alnscoresum/float(numaln+numunaln)
##avgalnscore_per_aln = sumalnscore_per_aln / float( totalaln )
##
##
##avgedit = editsum/float(alnlensum)
###editalnratio = editsumwithunmapped / float(alnlensum)
###editmatchratio = editsumwithunmapped / (float(alnlensum) - editsum)
##
####matches = readlensum - unmappedreadlensum - editsum ## erroneously counts unaligned portions of aligned reads as matches.... (i.e. readlen - alnlen gives unaln portion of read)
#### All below agree -- only need 1 for printing -- e.g. matches and pctmatches
##matches = alnlensum - editsum   
##pctmatches = float(matches)/readlensum
##matches2 = readlensum - totaleditsum
##pctmatches2 = float(matches2)/readlensum
##pctmatchsum = float(matchsum)/readlensum
##pctmatchsum2 = float(matchsum2)/readlensum
##
##pctalnlen = float(alnlensum)/readlensum
##
###totaleditsum-unmappedreadlensum
###readlensum-unmappedreadlensum
##pctedit_allreads = float(totaleditsum)/(readlensum)
##pctmatch_allreads = float(matches)/(readlensum)
##pctedit_alnreads = float(totaleditsum-unmappedreadlensum)/(readlensum-unmappedreadlensum)
##pctmatch_alnreads = float(matches)/(readlensum-unmappedreadlensum)
##pctedit_alns = editsum/float(alnlensum)
##pctmatch_alns = float(matches)/float(alnlensum)
##
##
##pct_0_ctg = num_0_ctg/float(num_0_ctg + num_1_ctg + num_multi_ctg)
##pct_1_ctg = num_1_ctg/float(num_0_ctg + num_1_ctg + num_multi_ctg)
##pct_multi_ctg = num_multi_ctg/float(num_0_ctg + num_1_ctg + num_multi_ctg)
##


