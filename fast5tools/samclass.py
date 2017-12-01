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
        self._uses_M = None

    def __str__(self):
        ## Return string_from_parsed b/c it updates may have been made
        return self.string_from_parsed_alignment()

    def parse_alignment(self):
        aln = self.alignment.strip().split('\t')
        self.parsed_aln['QNAME'] = aln[0]
        self.parsed_aln['FLAG'] = int(aln[1])
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

    def uses_M(self):
        '''Assumes alignment present...'''
        ## Note... this fails silently in that if there are no MX= in cigar string, self._uses_M stays as None
        if self._uses_M is None:
            if 'M' in self.get_cigar_field():
                self._uses_M = True
            elif 'X' in self.get_cigar_field() or '=' in self.get_cigar_field():
                self._uses_M = False
        return self._uses_M

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

    def get_num_cigar_matches_and_mismatches(self):
        assert self.uses_M() == True
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['M']

    def get_num_cigar_insertions_to_ref(self):
        ''' i.e. deletions from read'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['I']

    def get_num_cigar_deletions_from_ref(self):
        ''' i.e. insertions to read'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['D']

    def get_num_cigar_hard_clips(self):
        ''' H can only be present as the first and/or last operation.'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['H']

    def get_num_cigar_soft_clips(self):
        ''' S may only have H operations between them and the ends of the CIGAR string.'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['S']

    def get_num_cigar_skipped_region_from_ref(self):
        ''' '''
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['N']

    def get_num_cigar_padded_silent_deletions_from_ref(self):
        ''' The SAM format is typically used to describe alignments against an unpadded reference sequence,
            but it is also able to describe alignments against a padded reference. In the latter case, we say we are using
            a padded SAM. A padded SAM is a valid SAM, but with the difference that the reference and positions in use are padded.
            There may be more than one way to describe the padded representation. We recommend the following.

            In a padded SAM, alignments and coordinates are described with respect to the padded reference sequence. 
            Unlike traditional padded representations like the ACE file format where pads/gaps are recorded in reads using *'s, 
            we do not write *'s in the SEQ field of the SAM format.12 Instead, we describe pads in the query sequences as deletions 
            from the padded reference using the CIGAR 'D' operation. In a padded SAM, the insertion and padding CIGAR operations 
            ('I' and 'P') are not used because the padded reference already considers all the insertions.

            Above copies from: https://samtools.github.io/hts-specs/SAMv1.pdf
            More info there.
            '''
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['P']


    def get_num_cigar_matches(self):
        assert self.uses_M() == False
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['=']

    def get_num_cigar_mismatches(self):
        assert self.uses_M() == False
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.cigar['X']

    

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
        return self.get_adjusted_pos(initialpos, clip, extra, direction)

    def get_clipping_adjusted_end(self, extra=0, use_real_clip_sizes=True):
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

    def get_genomic_window_around_alignment(self, flank=0, adjust_for_clipping=True, identifier=None):
        '''Get genomic window surrounding an alignment with some buffer/flanks -- e.g. for local re-alignment.'''
        ''' Default: W/o buffer/flanking sequence, it gives the clipping-adjusted guesses for start and end positions.'''
        ''' Add buffer/flank lengths to each side in two ways:'''
        '''     (1) int > 1 adds/subtracts that int.
                (2) float [0,1] adds/subtracts that proportion of read length
                    NOTE: 1.0 gives proportion while 1 gives 1 bp.'''
        '''1-based, closed.'''
        ''' turn adjust_for_clipping off to get genomic window around aligned portion only -- if off and flank=0, just gives aligned portion (MD=x)'''
        '''In splitread class - can check for overlap among various genomic windows from split alns of same read.'''
        assert flank >= 0 ## no negative numbers are meaningful here
        if flank == 0:
            extra = 0
        elif flank > 1 or (flank == 1 and isinstance(flank, int)):
            extra = int(flank) 
        elif flank >= 0 and flank < 1 or (flank == 1 and isintance(flank, float)):
            extra = int(self.get_proportion_of_read_len(flank))
        chrom = self.get_rname_field()
        if adjust_for_clipping:
            start = self.get_clipping_adjusted_start(extra=extra)
            end = self.get_clipping_adjusted_end(extra=extra)
        else:
            start = self.get_pos_field() - extra
            end = self.get_reference_end_pos() + extra
        if identifier is not None:
            return (chrom, start, end, identifier)
        else:
            return (chrom, start, end)


    def reference_strand(self):
        # With help understanding the bitwise flags from: http://blog.nextgenetics.net/?e=18
        if self.get_flag_field() & 16:
            return -1
        else:
            return 1

    def on_positive_strand(self):
        return self.reference_strand() > 0

    def on_negative_strand(self):
        return self.reference_strand() < 0







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
            return self.get_next_set_of_sam_records()
##            self.get_next_set_of_sam_records()
##            return self.records
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
        ''' samrecords is a list of SAM records describing one source read.
            This list is typically made by parsing a sorted-by-readname SAM file with SamSplitAlnAggregator() class.'''
        self.records = samrecords
        self.num_aln = None
        self.readname = None
        self.fast5info = None
        self.readlength = None

    def get_read_name(self):
        if self.readname is None:
            self.readname = self.records[0].get_qname_field()
            for record in self.records:
                assert record.get_qname_field() == self.readname
        return self.readname

    def get_read_length(self):
        if self.readlength is None:
            self.readlength = self.records[0].get_read_len()
            for record in self.records:
                assert record.get_read_len() == self.readlength
        return self.readlength

    def get_fast5_info(self):
        if self.fast5info is None:
            self.fast5info = self.records[0].get_fast5_field()
            for record in self.records:
                assert record.get_fast5_field() == self.fast5info
        return self.fast5info

    def get_mapq_fields(self):
        return [record.get_mapq_field() for record in self.records]
    
    def get_AS_fields(self):
        return [record.get_AS_field() for record in self.records]

    def get_edit_dist_fields(self):
        return [record.get_edit_dist_field() for record in self.records]

    def get_pct_identity(self):
        '''As far as I can tell, each NM:i: editdist (in BWA) only applies to the part of read that maps to underlying reference POS to POS+sum(MDX=).
            If one looks at number of soft clips or hard clips, both often exceed NM -- therefore, they must not be included in it.
            So pct identity of an individual aligned region is 100*NM/read.get_SEQ_len_without_clipped_regions().
            And the pct identity of all aligned regions is 100*sum(NM)/sum(read.get_SEQ_len_without_clipped_regions())

            Wait:
            EDIT DIST = NM = n_mismatch + n_del + n_ins = number of changes needed to transform string1 to string2
            But qlen = n_match + n_mismatch + n_ins
            And rlen = n_match + n_mismatch + n_del
            So if you have a cigar of 5M5I5M and another of 5M5D5M, and NM=6 for each,
                then NM/qlen != NM/rlen b/c qlen != rlen (not necessarily anyway)
            Few things:
                n_mismatch = NM - n_del - n_ins
                So we can know that.
                Therefore, we can know n_match as M-n_mismatch (where M is M from cigar string).
                n_match = M - (NM - n_del - n_ins)
                so pct_match could be n_match/(n_match+n_mismatch) = (M-(NM-n_del-n_ins))/M
                But pct_match ignores indels...
                One could get pct of bases in read that form matches by:
                pct_read = n_match/read_len = (M-(NM-ndel-nins))/read_len

            Blast pct id?
            If you give blast: Query v Subject (below are fake used to convey small point)
            If you give blast: AAAAA v AAAAA: Identities = 5/5 (100%)
            If you give blast: AAGAA v AAAAA: (mm) Identities = 4/5 (80%)
            If you give blast: AAGAAA v AAAAA: (ins) Identities = 5/6 (83.3%)
            If you give blast: AAAA   v AAGAA: (del) Identities = 4/5 (80%)
            In general it is: number_matches / local_alignment_length
            ...where local_aln_length = n_match + n_mismatch + n_ins + n_del = M+I+D in cigar-speak
            So technically one cannot 'exactly' get pct identity across whole read if there are H and S present... one would need to know how many MDI it would form.
            But I do show a way to approximate it below b/c while pcr_identity as described above is great to talk about the
                aligned portions of a read, it is not fantastic for talking about the quality of the entire read.
            
            To ponder:
            One can get a normalized EDIT by:
                EDIT/CIGAR_LEN = (n_mm + n_ins + n_del) / (n_mm + n_ins + n_del + n_match)
            One can get normalized edit for read across all split alignments by:
                Sum(EDITS)/Sum(CIGAR_LENS)
            However, this of course might ignore the fact that large portions of the read did not align -- i.e. clipping.
            Thus, a 100kb read with 1 kb aligned and Edit=100 would get same score as a fully aligned 1 kb read with Edit=100.
            As far as the aligned sections go, we should consider them equivalent.
            As far as the entirety of the reads, there clearly is more wrong with the 100kb read than the 1kb read.
            So if one wants more info about the full read need to consider its full length and how much is aligned vs not.
            
            One can get a read_length normalized EDIT by:
                sum(NM)/readlen

            However, this could is overly optimistic for reads with LOTS of clipping.
            The normalized Edit here would seem smaller for that 100 kb read than the 1 kb read described above.
            This is b/c it is not incorporating unaligned info.
            
            Another way would be to do:
                total_distance = (sum(NM)+sum(unused_clipping))
                n_EDIT = (sum(NM)+sum(unused_clipping))/readlen   -- where unused clipping is only amount of clipped bases that do not appear as aligned bases in another one of mutliple split alignments
                       = total_distance/readlen
            This way would clearly judge the 100kb read as worse... which is more in line with what we'd want to know.
            This one might be overly pessimistic -- since some of the clipped stuff might have been mathces -- but I wouldn't think by too much.
                -> pessimistic b/c:  some of the clipped stuff might have been mathces
                -> pessimistic b/c:  the cigar for the clipped regions if they were to align as well as across all the aligned regions is likely longer than the readlen
                    ..therefore the normalized edit appears higher than it would be given the fully aligned known cigar
                -> one could potentially remedy this by using cigar lengths
                    adjusted_read_len = cigar_of_all_aligned_regions+sum(unused_clipping)
                    adjusted_read_len = bst guess at total alignment len
                best way to estimate normalized EDIT:
                    n_EDIT = (sum(NM)+sum(unused_clipping))/(sum(cigar_of_all_aligned_regions)+sum(unused_clipping))
                           = (sum(NM)+sum(unused_clipping))/
                           = total_distance/total_alignment

            Turn this around:
                num_matches = total_alignment - total_distance
                            = (sum(cigar_of_all_aligned_regions) + sum(unused_clipping)) - (sum(NM)+sum(unused_clipping))
                            = (sum(cigar_of_all_aligned_regions) - sum(NM)
                            = (n_match + n_mismatch + n_ins + n_del) - (n_mismatch + n_ins + n_del)
                            = n_match
                pct_identity_including_clipped_regions = (total_alignment - total_distance)/total_alignment
            
            
            As a proxy, I used to use:
                (read_len-NM)/readlen
            However, I guess I could do better as above.
            
        '''
        

    def get_num_aln(self):
        ''' This is more accurately "number of sam records".
            However, I left it as is, since the functions that use it will
            typically be using it when 1+ sam records is an alignment.
            This will report unaligned SAM records as 1 alignment.... so be careful to use with "has_alignments()"'''
        if self.num_aln is None:
            self.num_aln = len(self.records)
        return self.num_aln

    def has_alignments(self):
        ''' Looks for at least one record with an alignment.'''
        for i in range(self.get_num_aln()):
            if self.get_record(i).is_aligned():
                return True
        return False

    def has_only_records_that_are_aligned(self):
        ''' This ensures all SAM records are for alignments - not unaligned reads.
            This is essentially redundant with has_alignments() since it appears the possibilities for a SAM record are:
                1. unaligned read
                2. aligned read - only alignment
                3. aligned read - 1 of N alignments'''
        count = 0
        for i in range(self.get_num_aln()):
            if self.get_record(i).is_aligned():
                count+=1
        return self.get_num_aln() == count


    def get_record(self,index):
        return self.records[index]
    
    def genomic_window_ovlps(self):
        ## Use pybedtools to find overlaps
        pass


    def determine_longest_alignment(self):
        if self.has_alignments():
            maxalnlen = 0
            for i in range(self.get_num_aln()):
                if self.get_reference_aln_len() > maxalnlen:
                    maxalnlen = self.get_reference_aln_len()
                    maxalnidx = i
            return (maxalnlen, maxalnidx)
        else:
            return None

    def determine_highest_alignment_score(self):
        if self.has_alignments():
            maxalnscore = float('-inf')
            for i in range(self.get_num_aln()):
                AS = self.get_record(i).get_AS_field()
                if AS > maxalnscore:
                    maxalnscore = AS 
                    maxalnidx = i
            return (maxalnscore, maxalnidx)
        else:
            return None

    def determine_highest_alignment_score_ratio(self):
        ''' How much bigger is highest AS than next biggest AS?'''
        if self.has_alignments():
            maxalnscore = float('-inf')
            secondplace = float('-inf')
            for i in range(self.get_num_aln()):
                AS = self.get_record(i).get_AS_field()
                if AS > maxalnscore:
                    secondplace = maxalnscore
                    maxalnscore = AS
                elif AS > secondplace:
                    secondplace = AS
            if secondplace is not float('-inf'):
                return float(maxalnscore)/secondplace
            else:
                return float('inf')
        else:
            return None

    def ovlp(self, a,b,d=0):
        ''' a and b are 3-tuples of (chr,start,end)'''
        ''' assumes sorted such that, if ovlp, a is before b along number line.'''
        ''' d is distance - allows a gap up to d between intervals to still be an overlap - default 0'''
	if a[0] == b[0] and a[1] <= b[1] and b[1] <= a[2]:
            return True
	return False

    def merge(self, l, d=0):
        ''' Takes in list of 3-tuples of (chr.start,end)'''
        ''' Returns updated list where overlapping intervals are merged'''
        ''' d is distance - allows a gap up to d between intervals to still be an overlap - default 0'''
        a = l[0]
##        print 1, a
        update = []
        for i in range(1,len(l)):
            b = l[i]
            if self.ovlp(a,b,d):
                a = self.merge_tuples(a,b) ### new merge_tuples is identifier-aware; old = (a[0],a[1],b[2])
            else:
                update.append(a)
                a = b
        ## at end, if last thng ovlpd, then it merged interval in the making (a) needs to added to 'update'
        ##         if last thing didn't ovlp, then final interval, named 'a' before exiting, needs to be added to 'update'
        update.append(a)
##        print l,  update
        return update

    def merge_tuples(self, a,b):
        if len(a) == 3:
            return (a[0],a[1],b[2])
        elif len(a) == 4: ##has identifier
##            print a
            identifier = str(a[3]) + "," + str(b[3])
            return (a[0], a[1], b[2], identifier)


    def get_indexes_from_merged_tuples(self, genomic_windows):
        ''' After taking genomic_windows with indexes corresponding to which SAM record each is from,
            followed by merging these tuples, one gets a list of (possibly) merged tuples with comma-separated identifiers (in 4th position).
            This returns two answers.
            1. List of indexes for tuples with a merge.
            2. For each merged tuple, list of indexes of the SAM records that were merged.
            These answers are the keys and values of a single dictionary:
            {merge_idx:[sam_Record_idxs]}
            If there are no merges, an empty dict is returned.'''
        idx = {}
        for i in range(len(genomic_windows)):
##            print genomic_windows
            sam_idx = [int(e) for e in str(genomic_windows[i][3]).split(",")]
            if len(sam_idx) > 1:
                idx[i] = sam_idx
        return idx

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

    def determine_majority_window(self, genomic_windows, majority=0.5, scale_factor=False, proportions_provided=False):
        ''' Meant to look at list of genomic_window tuples and report whether one window
            makes up the majority percentage of the summed length of all windows. (reports index)
            This does not mean the largest proportion (which would be same index as longest window).
            This needs the max window to be larger than a given cutoff.
            Majority is not simply 51% - it should be some number that makes sense.
            For example: 0.7 is saying I'd like at least 70% of the total genomic window to be encompassed by one.
            ...So long as majority >0.5, there should be only one gw index returned.
            ...Nonetheless, it is coded anticipating other uses that may allow more than one window index returned.
            If a single window = 0.5 and there are >2 windows, then it is majority...'''
        gw_props = genomic_windows
        if not proportions_provided: #default is that proportions are obtained herein, but it makes more sense for some use cases to create/store the proportions array for other purposes in which case it can be provided.
            gw_props = self.determine_window_proportions(genomic_windows, scale_factor)
        gw_max = gw_props.max()
        if gw_max >= majority:
            return np.array(range(len(gw_props)))[gw_props == gw_max]
        else:
            return False


    def determine_longest_window_ratio(self, genomic_windows, proportions_provided=False):
        ''' Meant to look at list of genomic_window tuples and report whether the ratio between
            the longest and second longest window.
            Don't need proportions here.... same answer when not scaled to sum... but I wanted the np.array'''
        gw_props = genomic_windows
        if not proportions_provided: #default is that proportions are obtained herein, but it makes more sense for some use cases to create/store the proportions array for other purposes in which case it can be provided.
            gw_props = self.determine_window_proportions(genomic_windows, scale_factor)
        if len(gw_props) == 1:
            return float('inf')
        else:
            gw_props = sorted(gw_props)[-1::-1] #largest to smallest
            return gw_props[0] / float(gw_props[1])

    def alignments_ordered_like_read(self, indexes=[]):
        ''' ASSUMES all are on the same chrom... though violating this will likely return False
                In future: can have it ensure that the records actually occur on same chrom (just to be safe)...
            indexes = list of indexes of the records you'd like to check.
                e.g. you might have 3 split alignments but only want to check 2 and 3.
            Checks to see if the 5' clipping in reads agrees with POS field ordering.
            BWA gives partitioned split alignments. So if 3 split alignments occur consecutively along
            a genomic region, then their 5' clipping lengths should each be longer than the last.'''
        if not indexes:
            indexes = range(self.get_num_aln())
        d = {'pos':[], 'clip':[]}
        for i in indexes:
            d['pos'].append( self.get_record(i).get_pos_field() )
            d['clip'].append( self.get_record(i).get_5prime_clip_len() )
        df = pandas.DataFrame(d)
        ans = df.sort('pos')['pos'] == df.sort('clip')['pos']
        return sum(ans) == len(indexes)


    def strands_of_alignments(self, indexes=[]):
        if not indexes:
            indexes = range(self.get_num_aln())
        strands = np.zeros(len(indexes))
        for i in range(len(indexes)):
            strands[i] = self.get_record(indexes[i]).reference_strand() 
        return strands
    
    def alignments_on_same_strand(self, indexes=[]):
        ''' ASSUMES all are on the same chrom... though violating this will likely return False
                In future: can have it ensure that the records actually occur on same chrom (just to be safe)...
            indexes = list of indexes of the records you'd like to check.
                e.g. you might have 3 split alignments but only want to check 2 and 3.
            Checks to see if the 5' clipping in reads agrees with POS field ordering.
            BWA gives partitioned split alignments. So if 3 split alignments occur consecutively along
            a genomic region, then their 5' clipping lengths should each be longer than the last.'''
        strands = self.strands_of_alignments(indexes)
        positive = sum(strands == 1)
        negative = sum(strands == -1)
        allsame = (positive == 0 and negative == len(indexes)) or (negative == 0 and positive == len(indexes))
        return allsame, positive, negative

    def get_spanning_alignment_coordinates(self, indexes=[]):
        ''' Given list of indexes specifying a (sub)set of SAM records,
            return (1-based, closed) 3-tuple of:
                (chrom, 5'-most start position, 3'-most end position).
            This constitutes the reference stretch spanned by these alignments.
            If all records do not share the same RNAME field, None is returned.
            In context of alignments captured in a genomic window similar to the read length (+/-)
                that are shown to be ordered along the read correctly and the alignments are on the same strand,
                it is plausible (though not guaranteed) that the DNA molecule that gave rise to the read
                spanned at least this stretch of the reference --- with segments of bad base-calling splitting up the alignment.
                NOTE: spurious alignments can add noise to the order and strand info so be careful in rejecting real things...'''
        if not indexes:
            indexes = range(self.get_num_aln())
        d = {'start':[], 'end':[]}
        # initialize
        d['start'].append( self.get_record(indexes[0]).get_pos_field() )
        d['end'].append( self.get_record(indexes[0]).get_reference_end_pos() )
        chrom = self.get_record(indexes[0]).get_rname_field()
        for i in indexes[1:]:
            if self.get_record(i).get_rname_field() != chrom:
                return None
            d['start'].append( self.get_record(i).get_pos_field() )
            d['end'].append( self.get_record(i).get_reference_end_pos() )
        start = sorted(d['start'])[0]
        end = sorted(d['end'])[-1]
        return (chrom, start, end)

    def get_merge_with_max_spanning_alignment(self, idx, require_order=False, require_strand=False):
        ''' idx = doctionary output of get_index_from_merged_tuples().
            require_order = when True, output from alignments_ordered_like_read() must be True to be considered a valid merge.
            require_strand = when True, output from alignments_on_same_strand() must be True to be considered a valid merge.
            Only 'valid merges' are considered for max spanning alignment:
                all merges by default,
                filtered with require_order and require_strand as described above.
            Returns index of max_merge and the max span len.
            max_span is set to False when there are no valid merges.'''
        max_span_len = 0
        same_order_as_read = True ## Set to True when require_order=False -- although may not be True. Can be detected as False if require_order=True.
        same_strand = True ## Set to True when require_strand=False -- although may not be True. Can be detected as False if require_strand=True.
        max_span_i = None
        for merge_i in idx.keys():
            sam_idx_list = idx[merge_i]
            span_coords = self.get_spanning_alignment_coordinates(sam_idx_list)
            merge_span_len = span_coords[2] - span_coords[1] + 1
            if require_order:
                same_order_as_read = self.alignments_ordered_like_read(sam_idx_list)
            if require_strand:
                same_strand = self.alignments_on_same_strand(sam_idx_list)
            if merge_span_len > max_span_len and same_order_as_read and same_strand:
                max_span_len = merge_span_len
                max_span_i = merge_i
##        print 111, max_span_i
##        print 111, max_span_len
        return (max_span_i, max_span_len)
                

                

    def get_genomic_window(self, flank=0.1, merge_dist=0, majority=0.5, require_order=False, require_strand=False):
        '''Returns 3-tuple'''
        '''flank = as in get_genomic_window_around_alignment() described below'''
        ''' merge_dist = d from self.merge(): allows a gap up to d between intervals to still be an overlap - default 0
            majority = threshold to exceed to be consider a majority.
            require_order = as described in get_merge_with_max_spanning_alignment()
            require_strand = as described in get_merge_with_max_spanning_alignment()'''
        ''' Passage from the individual SamRecord function get_genomic_window_around_alignment():'''
        '''Get genomic window surrounding an alignment with some buffer/flanks -- e.g. for local re-alignment.'''
        ''' Default: W/o buffer/flanking sequence, it gives the clipping-adjusted guesses for start and end positions.'''
        ''' Add buffer/flanks two ways:'''
        '''     (1) int > 1 adds/subtracts that int.
                (2) float [0,1] adds/subtracts that proportion of read length
                    NOTE: 1.0 gives proportion while 1 gives 1 bp.'''
        '''1-based, closed.'''
        '''In splitread class - can check for overlap among various genomic windows from split alns of same read.'''
        if not self.has_alignments():
            print -1
            return None ## No alignments, therefore no genomic window.
        elif self.get_num_aln() == 1: ## although get_num_aln() will report 1 for unaligned SAM records, those are taken care of above.
            print 0
            return self.get_record(0).get_genomic_window_around_alignment(flank=flank, adjust_for_clipping=True, identifier="single_alignment")
        elif self.get_num_aln() > 1:
            print 1
            genomic_windows = []
##            bedstr = ''
            for i in range(self.get_num_aln()):
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
                print 2,1
                #return the composite genomic window formed by overlapping genomic windows
                identifier = "single_merged_window"
                single_merged_window = (gw_merge[0][0], gw_merge[0][1], gw_merge[0][2], identifier)
                return single_merged_window
            elif num_after_merge > 1:
                print 2,2
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
                split_alignments = []
                for i in range(self.get_num_aln()):
                    ref_coords = self.get_record(i).get_genomic_window_around_alignment(flank=0, adjust_for_clipping=False, identifier=i) 
                    split_alignments.append( ref_coords )
                split_aln_lengths = self.determine_window_lengths(split_alignments)
                split_alignment_proportions = self.determine_window_proportions(split_alignments)
                split_alignment_majority = self.determine_majority_window(split_alignment_proportions, majority, proportions_provided=True)
                if split_alignment_majority:
                    print 3,1
                    chosen_alignment_index = split_alignment_majority[0] ## FOR NOW JUST TAKE FIRST IF MORE THAN ONE (should be rare)
                else:
                    print 3,2
                    chosen_alignment_index = self.determine_longest_window(split_alignments)[0] ## FOR NOW JUST TAKE FIRST IF MORE THAN ONE (should be rare)
                AS = self.get_record(chosen_alignment_index).get_AS_field()    
                split_alignment_top2_ratio = self.determine_longest_window_ratio(split_alignment_proportions, proportions_provided=True)
                split_alignment_highest_AS = self.determine_highest_alignment_score()[0]
                split_alignment_top2_AS_ratio = self.determine_highest_alignment_score_ratio()
                longest_aln = split_aln_lengths[chosen_alignment_index]
                split_alignment_label = 'numaln:' + str(self.get_num_aln()) + "|longest_aln_len:" + str(longest_aln) + "|longest_aln_proportion:" + str(split_alignment_proportions[chosen_alignment_index]) + "|Top2AlnLengthRatio:" + str(split_alignment_top2_ratio) + "|AS:" + str(AS) + "|LongestAStoHighestASratio:" + str( AS/float(split_alignment_highest_AS) ) + "|top2ASratio:" + str(split_alignment_top2_AS_ratio)
                chosen_alignment = genomic_windows[chosen_alignment_index]
                if num_after_merge == self.get_num_aln():
                    print 4,1
                    ## nothing merged in above conditions -- meaning alignments are relatively far away from each other
                    ## is there a majority? Aln that makes up > X% of summed alignment length?
                    if split_alignment_majority:
                        print 5,1
                        ## need to get genomic window around alignment for majority alignment
                        ## these types of genomic windows were obtained above in genomic_windows -- not genomic_alignments

                        identifier = 'majority_alignment|' +  split_alignment_label
                    else:
                        print 5,2
                        ## return longest ((Other options would be to: return multiple, use alignment scores to pick))
                        identifier = 'longest_alignment|' + split_alignment_label
                    chosen_alignment_genomic_window = (chosen_alignment[0], chosen_alignment[1], chosen_alignment[2], identifier)
                    return chosen_alignment_genomic_window
                else:
                    print 6
                    # compare merged and single alns all at same time?
                    #is there a merge that exceeds the max align length?
                    #   gather information:
                    #     does that merge contain alignments in same order as found along read?
                    #     are the alignments inside merge from the same strand?
                    #     is the sum of the alignments in that span also longer than longest single aln?
                    #     is the sum of the aln scores higher than the single aln score?
                    merge_idx_dict = self.get_indexes_from_merged_tuples(gw_merge)
                    # determine if a given merge is interesting:
                    # -- note require_strand and require_order should be false for now. Noise and other artifacts can cause the ignoring of an actually valid merge.
                    #   potentially better to just add that information to identifier-string        
                    merge_with_max_spanning_aln = self.get_merge_with_max_spanning_alignment(idx=merge_idx_dict, require_order=require_order, require_strand=require_strand)

                    if merge_with_max_spanning_aln[0] is not None: ## (i,j): i = False when no valid merges, merge_idx when there is; j = span_len
                        print 7,1
                        sam_idx_list = merge_idx_dict[merge_with_max_spanning_aln[0]]
                        same_order_as_read = self.alignments_ordered_like_read(sam_idx_list)
                        same_strand = (',').join( [str(e) for e in self.alignments_on_same_strand(sam_idx_list)] )
                        sum_aln = sum([split_aln_lengths[i] for i in sam_idx_list])
                        sum_aln_gt_longest = sum_aln > longest_aln
                        longest_in_merge = chosen_alignment_index in sam_idx_list
                        sum_AS = sum([self.get_record(i).get_AS_field() for i in sam_idx_list])
                        sum_AS_gt_longest_AS = sum_AS > AS
                        sum_AS_gt_highest_AS = sum_AS > split_alignment_highest_AS
                        merge_label = 'num_aln_in_merge:' + str(len(sam_idx_list)) + '|same_order_as_read:' + str(same_order_as_read) + '|same_strand:' + str(same_strand) + '|span_len:' + str(merge_with_max_spanning_aln[1]) + '|sum_mergealn_len:' + str(sum_aln) + '|sum_mergealn_len_gt_longest_aln:' + str(sum_aln_gt_longest) + '|longest_aln_in_merge:' + str(longest_in_merge) + '|sum_merge_AS:' + str(sum_AS) + '|sum_merge_AS_gt_AS_of_longest_aln:' + str(sum_AS_gt_longest_AS) + '|sum_merge_AS_gt_highest_AS:' + str(sum_AS_gt_highest_AS) 
                        if merge_with_max_spanning_aln[1] > longest_aln:
                            print 8,1
                            chosen_genomic_window = gw_merge[merge_with_max_spanning_aln[0]]
                            if split_alignment_majority:
                                print 9,1
                                identifier = 'longest_merge_span_gt_majority_aln|' + merge_label + '|' + split_alignment_label
                            else:
                                print 9,2, merge_with_max_spanning_aln
                                identifier = 'longest_merge_span_gt_longest_aln|' + merge_label + '|' + split_alignment_label
                        else:
                            print 8,2
                            chosen_genomic_window = chosen_alignment
                            if split_alignment_majority:
                                print 9,3
                                identifier = 'longest_merge_span_lt_majority_aln|' + merge_label + '|' + split_alignment_label
                            else:
                                print 9,4
                                identifier = 'longest_merge_span_lt_longest_aln|' + merge_label + '|' + split_alignment_label
                        chosen_genomic_window = (chosen_genomic_window[0], chosen_genomic_window[1], chosen_genomic_window[2], identifier)
                        return chosen_genomic_window
                    else:
                        print 7,2, 'same as', 4,1
                        ## REPEATED CODE BELOW....(repeated from above single_aln area)... can maybe have function?
                        ## No valid merge - choose from alignments
                        if split_alignment_majority:
                            ## need to get genomic window around alignment for majority alignment
                            identifier = 'majority_alignment|' +  split_alignment_label
                        else:
                            ## return longest ((Other options would be to: return multiple, use alignment scores to pick))
                            identifier = 'longest_alignment|' + split_alignment_label
                        chosen_alignment_genomic_window = (chosen_alignment[0], chosen_alignment[1], chosen_alignment[2], identifier)
                        return chosen_alignment_genomic_window

                
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
    ''' Length of SEQ (part of READ that is soft-clipped and/or aligned) as defined by SAM specifications. x = cigar string.'''
    return sum([int(e) for e in re.findall('(\d+)[MIS=X]', x)])

def sam_read_len(x):
    ''' Length of entire READ as defined by SAM specifications. x = cigar string.'''
    return sum([int(e) for e in re.findall('(\d+)[HMIS=X]', x)])

def sam_alnseq_len(x):
    ''' Length of READ that is aligned only (no soft clipping) as defined by SAM specifications. x = cigar string.'''
    return sum([int(e) for e in re.findall('(\d+)[MI=X]', x)])

def sam_refalnseq_len(x):
    '''Length of REFERENCE sequence underlying the aligned portion of read. x = cigar string.'''
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

##out = open('$PRE-per-read-stats.txt','w')
##out.write( 'numentries\t'+str(numentries)+'\n' )
##out.write( 'numuniqentries\t'+str(numuniqentries)+'\n' )
##out.write( 'totalnumaln\t'+str(totalnumaln)+'\n' )
##out.write( 'totalnumuniqaln\t'+str(totalnumuniqaln)+'\n' )

##out.write( 'numaln_with_un_should_be_0\t'+str(numaln_un)+'\n' )
##out.write( 'num_unaln\t'+str(numunaln)+'\n' )
##out.write( 'num_0_ctg\t'+str(num_0_ctg)+'\n' ) ## should be same as num_unaln
##out.write( 'num_1_ctg\t'+str(num_1_ctg)+'\n' ) ## should be same as totalnumuniqaln
##out.write( 'num_multi_ctg\t'+str(num_multi_ctg)+'\n' )
##out.write( 'readlensum\t'+str(readlensum)+'\n' )
##out.write( 'alnlensum\t'+str(alnlensum)+'\n' )
##out.write( 'unalnportionsum\t'+str(unalnportionsum)+'\n' )
##out.write( 'unmappedreadlensum\t'+str(unmappedreadlensum)+'\n' )
##out.write( 'editsum\t'+str(editsum)+'\n' )

##out.write( 'totaleditsum\t'+str(totaleditsum)+'\n' )
##out.write( 'sum_matches\t'+str(matches)+'\n' )
##
##
##
##out.write( 'pctalnlen\t'+str(pctalnlen)+'\n' )
##out.write( 'pctedit_allreads\t'+str(pctedit_allreads)+'\n' )
##out.write( 'pctmatch_allreads\t'+str(pctmatch_allreads)+'\n' )
##out.write( 'pctedit_alnreads\t'+str(pctedit_alnreads)+'\n' )
##out.write( 'pctmatch_alnreads\t'+str(pctmatch_alnreads)+'\n' )
##out.write( 'pctedit_alns\t'+str(pctedit_alns)+'\n' )
##out.write( 'pctmatch_alns\t'+str(pctmatch_alns)+'\n' )
##
##
##out.write( 'pct_reads_that_aln\t'+str(pctaln_old)+'\n' )
##out.write( 'alnratio\t'+str(alnratio)+'\n' )
##out.write( 'pct_0_ctg\t'+str(pct_0_ctg)+'\n' ) ## should be same as pct_reads_that_dont_align
##out.write( 'pct_1_ctg\t'+str(pct_1_ctg)+'\n' ) ## should be <= pct_reads_that_aln
##out.write( 'pct_multi_ctg\t'+str(pct_multi_ctg)+'\n' ) ## should be <= pct_reads_that_aln
##
##out.write( 'avgalnlen_allaln\t'+str(avgalnlen_allaln)+'\n' )
##out.write( 'avgalnlen_allalnread\t'+str(avgalnlen_allalnread)+'\n' )
##out.write( 'avgalnlen_allread\t'+str(avgalnlen_allread)+'\n' )
##
##out.write( 'summapq_per_aln\t'+str(summapq_old)+'\n' )
##out.write( 'avgmapq_per_aln\t'+str(avgmapq_old)+'\n' )
##out.write( 'mapqsum_per_read\t'+str(mapqsum)+'\n' )
##out.write( 'avgmapq_per_read\t'+str(avgmapq)+'\n' )
##out.write( 'avgmapq_with_unmap\t'+str(avgmapq_with_unmap)+'\n' )
##
##
##out.write( 'alnscoresum_per_aln\t'+str(sumalnscore_per_aln)+'\n' )
##out.write( 'avgalnscore_per_aln\t'+str(avgalnscore_per_aln)+'\n' )
##out.write( 'alnscoresum_aln_reads\t'+str(alnscoresum)+'\n' )
##out.write( 'avgalnscore_aln_reads\t'+str(avgalnscore)+'\n' )
##out.write( 'avgalnscore_allreads\t'+str(avgalnscore_with_unmap)+'\n' )
##
##
##out.write( 'mapqsum2\t'+str(mapqsum2)+'\n' )
##out.write( 'avgmapq2\t'+str(avgmapq2)+'\n' )
##out.write( 'alnscoresum2\t'+str(alnscoresum2)+'\n' )
##out.write( 'avgalnscore2\t'+str(avgalnscore2)+'\n' )
##
##
##out.write( '\t'+str()+'\n' )
##out.close()
##
##
#####################################
##
##out = open('$PRE-per-read-stats.simple.txt','w')
##
##out.write( 'alnlensum\t'+str(alnlensum)+'\n' )
##out.write( 'sum_matches\t'+str(matches)+'\n' )
##
##out.write( 'pctalnlen\t'+str(pctalnlen)+'\n' )
##out.write( 'pctmatch_allreads\t'+str(pctmatch_allreads)+'\n' )
##out.write( 'pctmatch_alnreads\t'+str(pctmatch_alnreads)+'\n' )
##out.write( 'pctmatch_alns\t'+str(pctmatch_alns)+'\n' )
##
##out.write( 'pct_reads_that_aln\t'+str(pctaln_old)+'\n' )
##out.write( 'alnratio\t'+str(alnratio)+'\n' )
##out.write( 'pct_multi_ctg\t'+str(pct_multi_ctg)+'\n' ) ## should be <= pct_reads_that_aln
##
##out.write( 'avgalnlen_allaln\t'+str(avgalnlen_allaln)+'\n' )
##out.write( 'avgalnlen_allalnread\t'+str(avgalnlen_allalnread)+'\n' )
##out.write( 'avgalnlen_allread\t'+str(avgalnlen_allread)+'\n' )
##
##out.write( 'avgmapq_per_aln\t'+str(avgmapq_old)+'\n' )
##out.write( 'avgmapq_per_read\t'+str(avgmapq)+'\n' )
##out.write( 'avgmapq_with_unmap\t'+str(avgmapq_with_unmap)+'\n' )
##
##out.write( 'avgalnscore_per_aln\t'+str(avgalnscore_per_aln)+'\n' )
##out.write( 'avgalnscore_aln_reads\t'+str(avgalnscore)+'\n' )
##out.write( 'avgalnscore_allreads\t'+str(avgalnscore_with_unmap)+'\n' )
##
##out.close()

