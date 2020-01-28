import os, sys, re, pybedtools, pandas, operator
from collections import defaultdict
import numpy as np
from string import maketrans

isthis = {'gt':operator.gt, 'ge':operator.ge, 'lt':operator.lt, 'le':operator.le, 'eq':operator.eq, 'ne':operator.ne} 

def Is(a,op,b):
    return isthis[op](a,b)



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

    def string_from_parsed_alignment(self):
        return ('\t').join( [ str(self.parsed_aln[k]) for k in self.get_sam_fields() ] )
    
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
    


    def get_cigar_counts(self, get_item=False):
        # if present as single '*', then CIGAR is unavailable
        # cigar chars are in 'MIDNSHP=X'
        # some ways to get seq and read lengths from CIGAR if not storing cigar dict
        #   self.seq_len = sum([int(e) for e in re.findall('(\d+)[MIS=X]', self.parsed_aln['CIGAR'])])
        #   self.read_len = sum([int(e) for e in re.findall('(\d+)[HMIS=X]', self.parsed_aln['CIGAR'])])
        if self.cigar is None:
            self.cigar = {"M":0, "I":0, "D":0, "N":0, "S":0, "H":0, "P":0, "=":0, "X":0}
            for count,char in re.findall('(\d+)([MIDNSHP=X])', self.parsed_aln['CIGAR']):
                self.cigar[char] += int(count)
        if not get_item:
            return self.cigar
        else:
            # expects item to be in MIDNSHP=X
            return self.cigar[get_item]

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

    def get_alignment_length(self):
        ''' Sum of all CIGAR operations'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return sum([self.cigar[e] for e in 'HSMID=XNP'])

    def get_alignment_length_MDI(self):
        ''' Sum of M+D+I = n_match + n_mismatch + n_ins + n_del
            No clipping included.
            No N or P operations included.
            =X included in case aligner used those instead of M'''
        if self.cigar is None: 
            self.get_cigar_counts()
        return sum([self.cigar[e] for e in 'MDI=X'])

    def get_num_matches_MDI_NM(self):
        if self.cigar is None: 
            self.get_cigar_counts()
        return self.get_alignment_length_MDI() - self.get_edit_dist_field()

    def get_num_mismatches_MDI_NM(self):
        if self.cigar is None: 
            self.get_cigar_counts()
        # RETURN: M - (MDI - NM) = M + NM - MDI = NM - DI
        # Below: M - n_match = n_mis+n_match - n_match
        return self.cigar['M'] - self.get_num_matches_MDI_NM()


##    def get_pct_identity_given_alignment(self):
##        ''' Returns pct id of local alignment as BLAST:
##                pct id  = 100.0 * n_match / (n_match + n_mismatch + n_ins + n_del)
##                        = 100.0 * (((n_match + n_mismatch + n_ins + n_del) - (n_mismatch + n_ins + n_del)) / (n_match + n_mismatch + n_ins + n_del)
##                        = 100.0 * (MDI - NM)/(MDI)
##                '''
##        if self.cigar is None: 
##            self.get_cigar_counts()
##        if self.uses_M():
##            MDI = pass
##        else:
##            pass

    def pos_field_is(self, op, b):
        return Is(self.get_pos_field(), op, b)

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

    def get_BC_field(self):
        '''Only works with SAM files that have the BC:Z: tag from fast5tools'''
        return self.find_given_extra_field(regex='BC:Z:(.*)', output_type=str)
 
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

    def get_start_pos_on_read(self):
        ''' Using 1-based, closed indexing.
            Where does alignment start on the read?
            Count up 5' clipping.
            Start = 5primeClipLen + 1'''
        return self.get_5prime_clip_len() + 1


    def get_end_pos_on_read(self):
        ''' Using 1-based, closed indexing.
            Where does alignment END on the read?
                End = 5primeClipLen + SEQ_len_without_clipped_regions  ## no + 1 b/c closed indexing
                    = start_pos_on_read + SEQ_len_without_clipped_regions - 1
            Assertion:
                This should give same answer as:
                End = ReadLen - 3primeClipLen
            '''
        end = self.get_5prime_clip_len()  + self.get_SEQ_len_without_clipped_regions()
        assert end == self.get_read_len() - self.get_3prime_clip_len()
        return end
    

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
        self.bcinfo = None
        self.readlength = None
        self.pctid = None
##        self.readcoords = None
        self.alnstring = None
        self.pct_aln = None
        self.pct_unaln = None
        self.alncounts = None
        self.pct_id = None

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

    def get_BC_info(self):
        ''' BC is fast5tools bardecoder BarCode field - BC:Z:_____'''
        if self.bcinfo is None:
            self.bcinfo = self.records[0].get_BC_field()
            for record in self.records:
                assert record.get_BC_field() == self.bcinfo
        return self.bcinfo

    def get_mapq_fields(self):
        return [record.get_mapq_field() for record in self.records]
    
    def get_AS_fields(self):
        return [record.get_AS_field() for record in self.records]

    def get_edit_dist_fields(self):
        return [record.get_edit_dist_field() for record in self.records]

    def get_SEQ_fields(self):
        return [record.get_SEQ_field() for record in self.records]

    def get_read_SEQ(self):
        ''' For 1 alignment, the read seq will typicall be there.
            For multiple alignments, sometimes the SEQ is clipped (with H) or not shown in secondary alignments.
            This attempts to find the longest SEQ that is ALSO the expected length.'''
        seqs = self.get_SEQ_fields()
        maxlen = 0
        keepseq = ''
        for seq in seqs:
            seqlen = len(seq)
            if seqlen > maxlen:
                maxlen = seqlen
                keepseq = seq
        assert maxlen == self.get_read_length()
        return keepseq

                
    
    def get_MDI_alignment_lengths(self):
        return [record.get_alignment_length_MDI() for record in self.records]

    def get_list_of_num_cigar_deletions_from_ref(self):
        return [record.get_num_cigar_deletions_from_ref() for record in self.records]


    
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

    def get_pct_identity(self, aligned_regions_only=False):
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
                    n_mismatch  = NM - n_del - n_ins
                                = NM - D - I
                So we can know that.
                Therefore, we can know n_match as M-n_mismatch (where M is M from cigar string).
                    n_match     = M - (NM - n_del - n_ins)
                                = M - (NM - D -I)
                                = M+D+I - NM
                                = MDI - NM
                note: you can use MDI to represent the sum of MDI from all alignments and NM to represnt sum NM of all alignments (( or both from one aln ))
                so pct_match could be:
                    pct_match   = 100 * n_match/(n_match+n_mismatch)
                                = 100 * (MDI - NM)/M
                But pct_match ignores indels by subtracting them out... it is telling you something about M
                
                One could get pct of bases in read that form matches by:
                    pct_matches_in_read    = 100 * n_match/read_len
                                = 100 * (M-(NM-ndel-nins))/read_len
                                = 100 * (MDI - NM)/read_len
                                ...where MDI and NM are sums from all split alignments if they exist
                This has some issues:
                    1. if there is overlap between split alignments, it can over-estimate the
                        number of bases in the read that are involved in match operations.
                    2. the denominator of read length instead of alignment length treats all non-matches as if they were mismatches.
                Can fix #2 by:
                    pct_matches_in_read_aln = pct matches in alignment
                                            = 100 * (MDI - NM)/MDI
                This clearly makes more sense since it translates to:
                                            = 100 * ( (n_match + n_mismatch + n_ins + n_del) - (n_mismatch + n_ins + n_del)/ (n_match + n_mismatch + n_ins + n_del)
                                            = 100 * n_match / (n_match + n_mismatch + n_ins + n_del)
                However, now parts of the read not involved in an alignment are ignored. These are clipping regions in all alignments.
                So we need to adjust this to:
                            pct_matches     = 100 * n_match / (n_match + n_mismatch + n_ins + n_del + n_unaligned)
                                            = 100 * (MDI - NM) / (MDI + n_unaligned)
                This has a few issues, but is best we can do for now (as also discussed below, coming at this from different angle):
                    1. Treats all unaligned bases as mismatches rather than collection of matches, mismatches, insertions, and deletions.
                        This likely results in a smaller denominator and therefore slightly larger %match.
                        One could also say that this treats the unaligned bases as insertions though, or mismatches and insertions (to reference).
                        Since they were involved only in clipping and were not aligned though, we miss out on the number of deletions that might be there...
                        An alternative could be to multiply the unaligned base length by 2.
                        This would treat the unaligned bases as an alternating set of insertions and deletions in the reference.
                        That operation would therefore likely over-estimate the alignment length if forced.
                        Probably thinking of them as mismatches or insertions makes most sense for now.
                    2. Overlapping sections between split alignments are given more weight in determining the pct_matches than
                        other non-overlapping regions of read...
                        This effect -- in pilot analyses -- seems small though. Determining average from NM and CIGAR alone is not possible.
                        Would need MD string - but not using it for now.

            I also thought of a scaling alternative to treating unaligned bases as mismatches:
                pct matches across read = pct matches in alignments scaled by pct of read that is aligned.
                                        =  100 * [(MDI - NM) / (MDI)] * [n_alignned/(n_aligned+n_unaligned)]
                                        = 100 * [(MDI-NM)*(n_aligned)] / [(MDI)*(n_aligned+n_unaligned)]
                                        = 100 * [n_match * n_aligned] / [MDI * readlen]
                                        = 100 * n_match * n_aligned / (n_match + n_mismatch + n_ins + n_del) * readlen
                                        = 100 * n_aligned/readlen * n_match/(n_match + n_mismatch + n_ins + n_del)
                                        = scalng pct_match to aligned proportion of read
                                        = penalizing as aligned proportion gets smaller
                                        = rewarding as aligned proportion gets bigger
                                        = returning same score when 100% of read is aligned.
                        This seems to typically return a lower (but not by much) estimate than treating unaligned bases as all mismatches (or insertions).
                        Scaling is typically >99% the same number. Nonetheless, that would suggest scaling is ever so slightly more conservative.
                        Both methods can converge at 100% match.
                        Some examples:
                            Case1:
                                100 bp read, 10 bp align with 10 matches, 90 bp do not align. EDIT=0 for alignment.
                                Local pct identity = n_match/MDI = 10/10 = 100%
                                Global pct identity for read
                                        Using MDI+unaligned_bases method = n_match/(MDI+90) = 10/(10+90) = 10/100 = 0.1 (10%)
                                        Using scaling method = (n_match/MDI) * aln/(aln+unaln) = 10/10 * 10/100 = 1*0.1 = 0.1 (10%)
                            Case 2:
                                100 bp read, 10 bp align with 5 matches, 5 mismatches, 90 bp do not align. EDIT=5 for alignment.
                                Local pct identity = n_match/MDI = 5/10 = 50%
                                Global pct identity for read
                                        Using MDI+unaligned_bases method = n_match/(MDI+90) = 5/(5+5+90) = 5/100 = 0.05 (5%)
                                        Using scaling method = (n_match/MDI) * aln/(aln+unaln) = 5/10 * 10/100 = 0.5*0.1 = 0.05 (5%)
                            Case3:
                                100 bp read, 15 bp align w/ 5 matches, 5 mismatches, 5 insertions, 5 deletions.
                                85 bp do not align. EDIT = 15. MDI = 20
                                Local pct identity = n_match/MDI = 5/20 = 25%
                                Global pct identity for read
                                        Using MDI+unaligned_bases method = n_match/(MDI+85) = 5/(20+85) = 5/105 = 0.04762 (4.762%)
                                        Using scaling method = (n_match/MDI) * aln/(aln+unaln) = 5/20 * 15/100 = 0.25*0.15 = 0.0375 (3.75%)
                                Scaling is more conservative b/c it is propagating information about deletions whereas MDI+unaln treats all unaligned
                                    bases as mismatches or insertions.
                                If one increases the number of deletions, MDI gets bigger while the proportion of read aligned stays the same.
                                    10 deletions (MDI=5+5+5+10=25, a/total=15/100 = M+I/total = 0.15):
                                        method1 = 5/(25+85) = 5/110 = 0.455
                                        method2 = 5/25 * 15/100 = 0.2 * 0.15 = 0.03
                                    15 deletions (MDI=5+5+5+15=30, a/total = M+I/total = 5+5+5 / 100 = 15/100):
                                        method1 = 5/(30+85) = 5/115 = 0.435
                                        method2 = 5/30 * 15/100 = 5/30 * 0.15 = 0.025 = 2.5%
                                Meanwhile, if one increases number of matches, mismatches, or insertions, MDI gets bigger AND proportion of read aligned gets bigger.
                                    10 insertions (MDI=5+5+10+5=25, a/total = M+I/readlen = 5+5+10/100 = 20/100)
                                        method1 = 5/(25+80) = 5/105 = 0.4762 (4.76%) ## same as above
                                        method2 = 5/25 * 20/100 = 0.2 * 0.2 = 0.04
                                    15 insertions (MDI=5+5+15+5=30, a/total = M+I/readlen = 5+5+15 = 25/100)
                                        method1 = 5/(30+75) = 5/105 = 0.4762 (4.76%) ## same as above since unaligned are treated as insertions (or mismatches) meaning this just switches where the Is are counted
                                        method2 = 5/30 * 25/100 = 5/30 * 0.25 = 0.04167 (4.167%)
                                Note that for DELETIONS, having 5, 10, or 15 deletions affects both methods by bringing pct_id down.
                                    However, for INSERTIONS, having 5, 10, or 15 insertions has no effect on resulting pct_id for METHOD_1.
                                        This is b/c the unaligned bases are treated as insertions anyway, so they could be added to the I in MDI.
                                    For METHOD_2, having more insertions in the aligned portion lowers the pct_id locally, and this is generalized out the read when scaled to aligned bases.
                                        On some level, that feels more "Bayesian" -- we believe the pct_id is lower after seeing 15 insertions than we do after seeing 5....
                                        So the scaling method seems to incorporate information about extent of deletions and insertions in the local alignment.
                                        Another way to interpret it is in a 100 bp read with 25 aligned bases, the best score you can get 25/25 matches in local and is 25/100 = 0.25 for global...
                                            If the local alignment has 100% identity, then it becomes 25/25 * 25/100 = 1*0.25 = 0.25
                                                                                                    = pct_id_seen * max_pct_id_possible
                                            It also converges on method1 here b/c having seen no indels, it is treating all unaligned bases as mismatches (or insertions)...
                                        Wait -- let me look at this again:
                                            Method2 insertions 5,10,15 gave pctids of: 3.75, 4.0, and 4.167.
                                                Somehow the pctid is getting better with more insertions...
                                            Do these two methods converge with higher levels of insertions? Likely so since all unaligned in method1 are treated as an insertion.
                                            So if method 2 is converging up to method1 despite more error... is method1 better?
                                            
                            Case 4: Do these methods converge when there are NO DELETIONS?
                                100 bp read, 15 bp align w/ 5 matches, 5 mismatches, 5 insertions, 0 deletions.
                                85 bp do not align. EDIT = 15. MDI = 15
                                Local pct identity = n_match/MDI = 5/15 = 33%
                                Global pct identity for read
                                        Using MDI+unaligned_bases method = n_match/(MDI+85) = 5/(15+85) = 5/100 = 0.05 (5%)
                                        Using scaling method = (n_match/MDI) * aln/(aln+unaln) = 5/15 * 15/100 = 0.05 (5%)
                                YES THEY DO....
                                Toggle mismatches and insertions
                                5 matches, 7 mismatches, 3 insertions, 0 deletions. MDI=15. aligned=15.
                                    Local = 5/15 = 33%
                                    Global1 = 5/(15+85) = 5% (does not change b/c insertions and mismatches treated same)
                                    Global2 = 5/15 * 15/100 = 5% (does not change)
                        A difference is definitely deletions....
                        How to make them agree with deletions?
                        Method 1 can be modified to be:
                            n_match/(MDI + unaligned + unaligned*D/MDI)
                            = n_match/(MDI + U*M/MDI + U*I/MDI + 2*U*D/MDI)
                        In 55/55 cases, this brought the two methods closer together w/ avg distances going from 0.31 to 0.042.
                        So -- overall, the methods differ in:
                            1. Whether expected deletion rate is modeled into final answer.
                            2. Whether observing more insertions (or mismatches) affects final answer.
                                




            Below comes to the pct_matches above through a different thought process I had - but both thought processes led to that answer.
            Blast pct id?
            If you give blast: Query v Subject (below are fake used to convey small point)
            If you give blast: AAAAA v AAAAA: Identities = 5/5 (100%)
            If you give blast: AAGAA v AAAAA: (mm) Identities = 4/5 (80%)
            If you give blast: AAGAAA v AAAAA: (ins) Identities = 5/6 (83.3%)
            If you give blast: AAAA   v AAGAA: (del) Identities = 4/5 (80%)
            In general it is: number_matches / local_alignment_length
            ...where local_aln_length = n_match + n_mismatch + n_ins + n_del = M+I+D in cigar-speak
            So technically one cannot 'exactly' get pct identity across whole read if there are H and S present in CIGAR... one would need to know how many MDI it would form.
            But I do show a way to approximate it below b/c while pct_identity as described above is great to talk about the
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
                pct matches in the alignments (with unused-clipped regions treated like mismatches)
            
                NOTE:   This assumes split alignments are non-overlapping, which I have found they aren't (i.e. they do overlap some times, and not just at clipped parts).
                        It is not super clear how to handle this scenario, but there are a few ways.
                        1. Keep calculation as is.
                            Keeping it as is somewhat gives the average over the overlapping sections - since it will normalize the matches in each to the length of each
                            - however, I guess it technically gives more weight/influence to the overlapping sections to the overall sum(matches)/sum(alnlens).
                        2. Identify the part of the read involved in both alignments - add the pct_match from each subregion and divide by 2 to get average.
                        3. Identify the part of read involved - and take the one with lower/higher pct matches.

                        #2 is clearly the ideal way... however, there is no quick way of knowing from NM and CIGAR, which Ms are matches and which are mismatches...
                        Therefore #3 will have to do.
                        Only a small % of read is involved in overlapping alignments anyway. Ones I checked <= 3.6% - though this could be a biased sample.


            
            As a proxy, I used to use:
                (read_len-NM)/read_len
                or
                (read_len-total_distance)/read_len
                I like this a little less than the complement (NM/read_len) b/c at least you can interpret NM/readlen as the number of
                    problems in the alignment per base in the read.
                It is unclear how to interpret (read_len-NM)/read_len to me,
                although since they are complements, they are equivalently useful for ranking and such.
                Perhaps it could be interpreted as the absolute minimum number of good bases in the read per bases in the read,
                    whereas (total_alignment - total_distance)/total_alignment 
                So it is like a conservative guess at how good the read is.

            Correlation of different methods (using 55 aligned reads):
            V1 = global method adding in unaligned bases
            V2 = global method of scaling local sum to pct aligned bases
            V3 = just the sum of locals (not accounting for unaligned portion of read)
            V4 = proxy above
                          V1        V2        V3        V4
                V1 1.0000000 0.9994857 0.9051552 0.7084775
                V2 0.9994857 1.0000000 0.8988643 0.7007634
                V3 0.9051552 0.8988643 1.0000000 0.9352874
                V4 0.7084775 0.7007634 0.9352874 1.0000000
                
            
        '''
        # RETURN:
        # GLOBAL: 100 * (MDI - NM) / (MDI + n_unaligned)
        # SUMLOCAL: 100 * (MDI - NM) / (MDI)
        # SCALED_LOCAL: pct_aligned_bases * (MDI - NM) / (MDI) = 100 * A/(U+A) * (MDI-NM) / NM 
        if self.pct_id is not None:
            return self.pct_id
        MDI = self.get_sum_MDI_alignment_lengths()
        n_match = MDI - self.get_total_edit_dist()
        U = self.get_number_bases_in_read_not_aligned()
        MDIU = MDI + U
        D = self.get_sum_cigar_deletions_from_ref()
        M = sum([record.get_cigar_counts('M') for record in self.records])
        I = sum([record.get_cigar_counts('I') for record in self.records])
        n_mismatch = MDI - D - I - n_match
##        Ud = U - U*M/float(MDI) - U*I/float(MDI)
        MDIUD = MDIU + U*D/float(MDI)
        PROP_U = U*M/float(MDI) + U*I/float(MDI) + 2*U*D/float(MDI)
        MDI_PROPU = MDI + PROP_U
        aln_only = 100.0 * n_match / MDI
        with_unaln = 100.0 * n_match / MDIU
        scaled_aln_only = self.get_pct_of_read_aligned() * n_match / MDI
        with_unaln_accounting_for_dels = 100.0 * n_match / MDIUD
        alt = 100.0 * n_match / MDI_PROPU
        self.pct_id = {'pctid':[with_unaln, scaled_aln_only, aln_only, with_unaln_accounting_for_dels, alt], 'match':n_match, 'mismatch':n_mismatch, 'del':D, 'ins':I, 'MDI':MDI, 'unaligned':U, 'MDIU':MDIU}
        #return with_unaln, scaled_aln_only, aln_only, with_unaln_accounting_for_dels, alt
        return self.pct_id

    def get_pct_identity_proxy(self):
        ''' Proxy is "(read_len-NM)/read_len" as described above in pct_identity blurb.
            This definitely seems to correlate with other pct_id estimates, but is less precise as explained in the pct_id blurb.
            Around 50% of the time it seems to be global_pct_id** <= proxy <= local_pct_id
                **obtained either by scaling to pct aligned bases or including unaligned bases in denom.
            When I looked at 55 reads that had alignments, 49.1% of the time it was the relationshp above, 29.1% it was hgher than the local, 21.8% of the time it was lower than the globals.
            As shown above, it has a correlation of 0.7 with the globals and 0.93 with local sum (using these 55 alignments).
            Overall - I'd just use one of the more formalized and interpretable pct_alns from now on...
                '''
        return 100.0 * (self.get_read_length() - self.get_total_edit_dist()) / self.get_read_length()

    def get_sum_MDI_alignment_lengths(self):
        return sum( self.get_MDI_alignment_lengths() )

    def get_total_edit_dist(self):
        return sum(self.get_edit_dist_fields())
    
    def get_total_dist(self):
        ''' as defined in get_pct_identity() blurb:
                total dist = total edit dist + total unaligned (not used in any of the alignments)
            Unaligned bases are effectively treated as mismatches.
            Also - alignments that have overlaps (typically not a lot of overlap) will'''
        return self.get_total_edit_dist() + self.get_number_bases_in_read_not_aligned()

    def get_sum_cigar_deletions_from_ref(self):
        return sum(self.get_list_of_num_cigar_deletions_from_ref())

    def get_coords_along_read(self, indexes=[]):
        ''' Assumes there alignments.
            Gets list of start and end positions of alignments along the read.
            Can return for any index or set of indexes given.
            Default is to return for all.'''
        assert self.has_alignments()
        if not indexes:
            indexes = range(self.get_num_aln())
        l = []
        for i in indexes:
            start =  self.get_record(i).get_start_pos_on_read() 
            end =  self.get_record(i).get_end_pos_on_read()
            l.append( (start,end,i) )
        return l

    def get_per_base_align_status_for_read(self, indexes=[]):
        ''' Reports CIGAR-like string pertaining to bases in the read regarding whether they are aligned or not, given a set of alignments.
            For example, if there is a 1000 bp read where the inner 800 bp is aligned
            and outer 100 bp on each side is not, return:
                100U800A100U
            Can also use on subset of indexes to obtain status of each bp given subset of alignments.

            Current implementation assumes that reads are nearly perfectly partitioned in the splits.
            That means it allows neighboring alignments to have overlap, but does not anticipate overlap
            between further away alignments. This seems to hold for BWA (bwa mem).

            Currently asserts that the sum of A+U should be the read length.
            This should be True - the assertion will fire up if the code is wrong (or an unanticipated case is found).'''

        ## Can be used on subsets of the alignments, but empty indexes does all and saves into a variable that is returned here if it already exists
        if not indexes and self.alnstring is not None:
            return self.alnstring
        

        #initialize
        alnstatus = ''
        lengths = {'Total':0, 'Un':0, 'Aln':0}

        coords = sorted( self.get_coords_along_read(indexes) )
        if coords[0][0] != 1: ## Then there is clipping/unaligned stuff at beginning
            unalnlen = coords[0][0] - 1
            lengths['Total'] += unalnlen
            lengths['Un'] += unalnlen
            alnstatus += str(unalnlen)+'U'
            unalnlen = 0
        alnlen = coords[0][1]-coords[0][0] + 1
        
        #iterate
        for i in range( 1, len(coords) ):
            if coords[i][0] > coords[i-1][1]:
                alnstatus += str(alnlen)+'A'
                lengths['Total'] += alnlen
                lengths['Aln'] += alnlen
                alnlen = 0
                unalnlen = (coords[i][0] - 1) - (coords[i-1][1] + 1) + 1
                lengths['Total'] += unalnlen
                lengths['Un'] += unalnlen
                alnstatus += str(unalnlen)+'U'
                unalnlen = 0
                alnlen += coords[i][1]-coords[i][0] + 1
                
            elif coords[i][0] <= coords[i-1][1]:
                alnlen += coords[i][1]-coords[i][0] + 1 - (coords[i-1][1]-coords[i][0]+1)
                
        #finalize
        if alnlen > 0:
            alnstatus += str(alnlen)+'A'
            lengths['Total'] += alnlen
            lengths['Aln'] += alnlen
            alnlen = 0
        rl = self.get_read_length()
        if coords[-1][1] < rl:
            unalnlen = rl - coords[-1][1]
            lengths['Total'] += unalnlen
            lengths['Un'] += unalnlen
            alnstatus += str(unalnlen)+'U'

        ## ASSERT: require the U+A len (T for total U+A) to equal the read length!
        assert rl == lengths['Total']
        assert lengths['Total'] == lengths['Aln']+lengths['Un']

        ## Can be used on subsets of the alignments, but empty indexes does all and saves into a variable that is returned here if it already exists
        if not indexes and self.alnstring is None:
            self.alnstring = alnstatus
            self.alncounts = lengths
        ## RETURN
        return alnstatus

    def get_number_bases_in_read_aligned(self):
        if self.alncounts is None:
            self.get_per_base_align_status_for_read()
        return self.alncounts['Aln']

    def get_number_bases_in_read_not_aligned(self):
        if self.alncounts is None:
            self.get_per_base_align_status_for_read()
        return self.alncounts['Un']
        

    def get_pct_of_read_aligned(self):
        ''' This is to offer a different perspective than the pct_id and stuff normalized to alignment lengths.
            This tells one what percent of a read is involved in an alignment (single or multiple splits) regardless
            of how many of those bases are in matches.'''
        if self.pct_aln is None:
             self.pct_aln = 100.0 * self.get_number_bases_in_read_aligned() / (self.get_number_bases_in_read_aligned() + self.get_number_bases_in_read_not_aligned() )
        return self.pct_aln
    
    def get_pct_of_read_not_aligned(self):
        ''' This is to offer a different perspective than the pct_id and stuff normalized to alignment lengths.
            This tells one what percent of a read is NOT involved in an alignment at all
            -- e.g. only involved in hard/soft clipping.'''
        if self.pct_unaln is None:
             self.pct_unaln = 100.0 * self.get_number_bases_in_read_not_aligned() / (self.get_number_bases_in_read_aligned() + self.get_number_bases_in_read_not_aligned() )
        return self.pct_unaln
    

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
            if secondplace not in (float('-inf'), 0):
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

##        print df.sort_values('pos')['pos']
##        print
##        print df.sort_values('clip')['pos']
##        print
        def getans(x,y):
            try:
                return x == y
            except:
                return [0]
        try:
            df.sort('pos')['pos']
            ans = getans(df.sort('pos')['pos'], df.sort('clip')['pos'])
        except: ## .sort deprecated
            ans = getans(df.sort_values('pos')['pos'], df.sort_values('clip')['pos'])
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

    def get_merge_with_max_spanning_alignment(self, idx, require_order=False, require_strand=False, return_max_span_aln_tuple=False):
        ''' idx = dictionary output of get_index_from_merged_tuples().
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
                max_span_aln_tuple = span_coords
##        print 111, max_span_i
##        print 111, max_span_len
        if not return_max_span_aln_tuple:
            return (max_span_i, max_span_len)
        else:
            return (max_span_i, max_span_len, max_span_aln_tuple)
                

                

    def get_genomic_window(self, flank=0.1, merge_dist=0, majority=0.5, require_order=False, require_strand=False, adjust_for_clipping_in_output=True):
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


        if self.has_alignments():
            alnstring = self.get_per_base_align_status_for_read()
            pctaln = self.get_pct_of_read_aligned()
            pctid = self.get_pct_identity()['pctid'][0]
            pctstr = 'pctid:' + str(pctid) + '|pctaln:' + str(pctaln) + '|alnstring:' + alnstring
        if not self.has_alignments():
##            print -1
            return None ## No alignments, therefore no genomic window.
        elif self.get_num_aln() == 1: ## although get_num_aln() will report 1 for unaligned SAM records, those are taken care of above.
##            print 0
            AS = self.get_record(0).get_AS_field()
            ## adjust_for_clipping needs to be "False" for "length" b/c determining alignment coords/length on reference
            length = self.determine_window_length( self.get_record(0).get_genomic_window_around_alignment(flank=0, adjust_for_clipping=False, identifier=0) )
            identifier = 'single_alignment|numaln:' + str(self.get_num_aln()) + '|longest_aln:' + str(length) + '|AS:' + str(AS) + '|' + pctstr
            ## adjust_for_clipping should be "True" below (and there ought to be some flank) if you are trying to extract a genomic window that definitely subsumes the read.
            ##      Therefore, it uses the variable that defaults to "True". Setting to "False" will just return the alignment coords.
            ##      return self.get_record(0).get_genomic_window_around_alignment(flank=flank, adjust_for_clipping=True, identifier=identifier)
            return self.get_record(0).get_genomic_window_around_alignment(flank=flank, adjust_for_clipping=adjust_for_clipping_in_output, identifier=identifier)
        elif self.get_num_aln() > 1:
##            print 1
            genomic_windows = []
##            bedstr = ''
            for i in range(self.get_num_aln()):
                ## adjust_for_clipping should be "False" below in order to see if the flanks inferred by clip lengths results in grouping more than 1 alignment.
                gw = self.get_record(i).get_genomic_window_around_alignment(flank=flank, adjust_for_clipping=True, identifier=i) 
                genomic_windows.append( gw )
##                bedstr = chrom + "\t" + str(gw[0]-1) + "\t" + str(gw[1])
            #bedtools approach
##            bedtool = pybedtools.BedTool( bedstr, from_string=True )
##            merged = bedtool.merge(d=merge_dist)
            #python-only approach
            sorted_genomic_windows = sorted(genomic_windows) ## use sorted not .sort(), so can access genomic_windows below
            gw_merge = self.merge(l=sorted_genomic_windows, d=merge_dist) ## merge uses overlap detection that requires sorted coords
            num_after_merge = len(gw_merge)
            if num_after_merge == 1:
##                print 2,1
                #return the composite genomic window formed by overlapping genomic windows
                
                #get aln information
                split_alignments = []
                for i in range(self.get_num_aln()):
                    ## adjust_for_clipping needs to be "False" for "length" b/c determining alignment coords/length on reference
                    ref_coords = self.get_record(i).get_genomic_window_around_alignment(flank=0, adjust_for_clipping=False, identifier=i) 
                    split_alignments.append( ref_coords )
                split_aln_lengths = self.determine_window_lengths(split_alignments)
                split_alignment_proportions = self.determine_window_proportions(split_alignments)
                chosen_alignment_index = self.determine_longest_window(split_alignments)[0] ## FOR NOW JUST TAKE FIRST IF MORE THAN ONE (should be rare)
                AS = self.get_record(chosen_alignment_index).get_AS_field()    
                split_alignment_top2_ratio = self.determine_longest_window_ratio(split_alignment_proportions, proportions_provided=True)
                split_alignment_highest_AS = self.determine_highest_alignment_score()[0]
                split_alignment_top2_AS_ratio = self.determine_highest_alignment_score_ratio()
                longest_aln = split_aln_lengths[chosen_alignment_index]
                LtoH_ASR = 1.0 if AS == float(split_alignment_highest_AS) else (float('inf') if split_alignment_highest_AS == 0 else AS/float(split_alignment_highest_AS))
                split_alignment_label = 'numaln:' + str(self.get_num_aln()) + "|longest_aln_len:" + str(longest_aln) + "|longest_aln_proportion:" + str(split_alignment_proportions[chosen_alignment_index]) + "|Top2AlnLengthRatio:" + str(split_alignment_top2_ratio) + "|AS:" + str(AS) + "|LongestAStoHighestASratio:" + str( LtoH_ASR ) + "|top2ASratio:" + str(split_alignment_top2_AS_ratio)
                chosen_alignment = genomic_windows[chosen_alignment_index]
                
                #get merge information
                merge_idx_dict = self.get_indexes_from_merged_tuples(gw_merge)
                #merge_with_max_spanning_aln = self.get_merge_with_max_spanning_alignment(idx=merge_idx_dict, require_order=require_order, require_strand=require_strand)
                merge_with_max_spanning_aln = self.get_merge_with_max_spanning_alignment(idx=merge_idx_dict, require_order=require_order, require_strand=require_strand, return_max_span_aln_tuple=True)
                sam_idx_list = merge_idx_dict[merge_with_max_spanning_aln[0]]
                same_order_as_read = self.alignments_ordered_like_read(sam_idx_list)
                same_strand = (',').join( [str(e) for e in self.alignments_on_same_strand(sam_idx_list)] )
                sum_aln = sum([split_aln_lengths[i] for i in sam_idx_list])
                sum_aln_gt_longest = sum_aln > longest_aln
                longest_in_merge = chosen_alignment_index in sam_idx_list
                sum_AS = sum([self.get_record(i).get_AS_field() for i in sam_idx_list])
                sum_AS_gt_longest_AS = sum_AS > AS
                sum_AS_gt_highest_AS = sum_AS > split_alignment_highest_AS
                merge_alnstring = self.get_per_base_align_status_for_read(indexes=sam_idx_list)
                mergealnstr_sam_as_alnstr = merge_alnstring == alnstring
                merge_label = 'num_aln_in_merge:' + str(len(sam_idx_list)) + '|same_order_as_read:' + str(same_order_as_read) + '|same_strand:' + str(same_strand) + '|span_len:' + str(merge_with_max_spanning_aln[1]) + '|sum_mergealn_len:' + str(sum_aln) + '|sum_mergealn_len_gt_longest_aln:' + str(sum_aln_gt_longest) + '|longest_aln_in_merge:' + str(longest_in_merge) + '|sum_merge_AS:' + str(sum_AS) + '|sum_merge_AS_gt_AS_of_longest_aln:' + str(sum_AS_gt_longest_AS) + '|sum_merge_AS_gt_highest_AS:' + str(sum_AS_gt_highest_AS) + '|merge_alnstring:' + merge_alnstring + '|merge_alnstring_sam_as_alnstring:' +  str(mergealnstr_sam_as_alnstr)

                identifier = 'single_merged_window|numaln:' + str(self.get_num_aln()) + '|' +  merge_label + '|' + split_alignment_label + '|' + pctstr
                
                if not adjust_for_clipping_in_output:
                    single_merged_window = (merge_with_max_spanning_aln[2][0], merge_with_max_spanning_aln[2][1], merge_with_max_spanning_aln[2][2], identifier)
                else:
                    single_merged_window = (gw_merge[0][0], gw_merge[0][1], gw_merge[0][2], identifier)
                return single_merged_window
            elif num_after_merge > 1:
##                print 2,2
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
                    ## adjust_for_clipping needs to be "False" for "length" b/c determining alignment coords/length on reference
                    ref_coords = self.get_record(i).get_genomic_window_around_alignment(flank=0, adjust_for_clipping=False, identifier=i) 
                    split_alignments.append( ref_coords )
                split_aln_lengths = self.determine_window_lengths(split_alignments)
                split_alignment_proportions = self.determine_window_proportions(split_alignments)
                split_alignment_majority = self.determine_majority_window(split_alignment_proportions, majority, proportions_provided=True)
                if split_alignment_majority is not False:
##                    print 3,1
                    chosen_alignment_index = split_alignment_majority[0] ## FOR NOW JUST TAKE FIRST IF MORE THAN ONE (should be rare)
                else:
##                    print 3,2
                    chosen_alignment_index = self.determine_longest_window(split_alignments)[0] ## FOR NOW JUST TAKE FIRST IF MORE THAN ONE (should be rare)
                AS = self.get_record(chosen_alignment_index).get_AS_field()    
                split_alignment_top2_ratio = self.determine_longest_window_ratio(split_alignment_proportions, proportions_provided=True)
                split_alignment_highest_AS = self.determine_highest_alignment_score()[0]
                split_alignment_top2_AS_ratio = self.determine_highest_alignment_score_ratio()
                longest_aln = split_aln_lengths[chosen_alignment_index]
                LtoH_ASR = 1.0 if AS == float(split_alignment_highest_AS) else (float('inf') if split_alignment_highest_AS == 0 else AS/float(split_alignment_highest_AS)) #AS/float(split_alignment_highest_AS)
                split_alignment_label = 'numaln:' + str(self.get_num_aln()) + "|longest_aln_len:" + str(longest_aln) + "|longest_aln_proportion:" + str(split_alignment_proportions[chosen_alignment_index]) + "|Top2AlnLengthRatio:" + str(split_alignment_top2_ratio) + "|AS:" + str(AS) + "|LongestAStoHighestASratio:" + str( LtoH_ASR ) + "|top2ASratio:" + str(split_alignment_top2_AS_ratio)
                chosen_alignment = genomic_windows[chosen_alignment_index]
                if num_after_merge == self.get_num_aln():
##                    print 4,1
                    ## nothing merged in above conditions -- meaning alignments are relatively far away from each other
                    ## is there a majority? Aln that makes up > X% of summed alignment length?
                    if split_alignment_majority is not False:
##                        print 5,1
                        ## need to get genomic window around alignment for majority alignment
                        ## these types of genomic windows were obtained above in genomic_windows -- not genomic_alignments

                        identifier = 'majority_alignment|' +  split_alignment_label + '|' + pctstr
                    else:
##                        print 5,2
                        ## return longest ((Other options would be to: return multiple, use alignment scores to pick))
                        identifier = 'longest_alignment|' + split_alignment_label + '|' + pctstr
                    ## chosen alignment is already adjusted for clipping OR NOT b/c it comes from "genomic_windows" which is
                    chosen_alignment_genomic_window = (chosen_alignment[0], chosen_alignment[1], chosen_alignment[2], identifier)
                    return chosen_alignment_genomic_window
                else:
##                    print 6
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
                    #merge_with_max_spanning_aln = self.get_merge_with_max_spanning_alignment(idx=merge_idx_dict, require_order=require_order, require_strand=require_strand)
                    merge_with_max_spanning_aln = self.get_merge_with_max_spanning_alignment(idx=merge_idx_dict, require_order=require_order, require_strand=require_strand, return_max_span_aln_tuple=True)

                    if merge_with_max_spanning_aln[0] is not None: ## (i,j): i = False when no valid merges, merge_idx when there is; j = span_len
##                        print 7,1
                        sam_idx_list = merge_idx_dict[merge_with_max_spanning_aln[0]]
                        same_order_as_read = self.alignments_ordered_like_read(sam_idx_list)
                        same_strand = (',').join( [str(e) for e in self.alignments_on_same_strand(sam_idx_list)] )
                        sum_aln = sum([split_aln_lengths[i] for i in sam_idx_list])
                        sum_aln_gt_longest = sum_aln > longest_aln
                        longest_in_merge = chosen_alignment_index in sam_idx_list
                        sum_AS = sum([self.get_record(i).get_AS_field() for i in sam_idx_list])
                        sum_AS_gt_longest_AS = sum_AS > AS
                        sum_AS_gt_highest_AS = sum_AS > split_alignment_highest_AS
                        merge_alnstring = self.get_per_base_align_status_for_read(indexes=sam_idx_list)
                        mergealnstr_sam_as_alnstr = merge_alnstring == alnstring
                        merge_label = 'num_aln_in_merge:' + str(len(sam_idx_list)) + '|same_order_as_read:' + str(same_order_as_read) + '|same_strand:' + str(same_strand) + '|span_len:' + str(merge_with_max_spanning_aln[1]) + '|sum_mergealn_len:' + str(sum_aln) + '|sum_mergealn_len_gt_longest_aln:' + str(sum_aln_gt_longest) + '|longest_aln_in_merge:' + str(longest_in_merge) + '|sum_merge_AS:' + str(sum_AS) + '|sum_merge_AS_gt_AS_of_longest_aln:' + str(sum_AS_gt_longest_AS) + '|sum_merge_AS_gt_highest_AS:' + str(sum_AS_gt_highest_AS) + '|merge_alnstring:' + merge_alnstring + '|merge_alnstring_sam_as_alnstring:' +  str(mergealnstr_sam_as_alnstr)
                        if merge_with_max_spanning_aln[1] > longest_aln:
##                            print 8,1
                            if not adjust_for_clipping_in_output:
                                chosen_genomic_window = (merge_with_max_spanning_aln[2][0], merge_with_max_spanning_aln[2][1], merge_with_max_spanning_aln[2][2])
                            else:
                                chosen_genomic_window = gw_merge[merge_with_max_spanning_aln[0]]
                            if split_alignment_majority is not False:
##                                print 9,1
                                identifier = 'longest_merge_span_gt_majority_aln|' + merge_label + '|' + split_alignment_label + '|' + pctstr
                            else:
##                                print 9,2, merge_with_max_spanning_aln
                                identifier = 'longest_merge_span_gt_longest_aln|' + merge_label + '|' + split_alignment_label + '|' + pctstr
                        else:
##                            print 8,2
                            chosen_genomic_window = chosen_alignment
                            if split_alignment_majority is not False:
##                                print 9,3
                                identifier = 'longest_merge_span_lt_majority_aln|' + merge_label + '|' + split_alignment_label + '|' + pctstr
                            else:
##                                print 9,4
                                identifier = 'longest_merge_span_lt_longest_aln|' + merge_label + '|' + split_alignment_label + '|' + pctstr
                        ## chosen_genomic_window is already adjusted for clipping OR NOT b/c it comes from "genomic_windows" which is
                        chosen_genomic_window = (chosen_genomic_window[0], chosen_genomic_window[1], chosen_genomic_window[2], identifier)
                        return chosen_genomic_window
                    else:
##                        print 7,2, 'same as', 4,1
                        ## REPEATED CODE BELOW....(repeated from above single_aln area)... can maybe have function?
                        ## No valid merge - choose from alignments
                        if split_alignment_majority is not False:
                            ## need to get genomic window around alignment for majority alignment
                            identifier = 'majority_alignment|' +  split_alignment_label
                        else:
                            ## return longest ((Other options would be to: return multiple, use alignment scores to pick))
                            identifier = 'longest_alignment|' + split_alignment_label
                        ## chosen_genomic_window is already adjusted for clipping OR NOT b/c it comes from "genomic_windows" which is
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




### MPILEUP FROM SAMTOOLS ANALYSIS OF BAM
class MpileupRecord(object):
    def __init__(self, rec, trantab, pattern, stranded=False):
        self.stranded=stranded
        self.trantab = trantab
        self.pattern = pattern
        self.seq, self.pos, self.ref, self.depth, self.composition, self.quals = rec.strip().split()
        if self.has_depth():
            ## This needs to catch mpileup records that report 0 depth and show "*" "*" in the reads and quals columns
            ## Those cause errors in _get_composition_counts - however, this is a band-aid over another problem
            ## If there are >= 1 reads but both show deletions for a base, then the marg2_totals will be divided by 0
            ## The answer to %mm or %m of 0mm and 0m is undefined 0/(0+0)
            ##  -- it is not abundantly clear how to handle it at the moment b/c reporting 0% mm for marg2 might make it seem like it has no errors
            ##  -- unless you also looked at %m which would also be 0% which would make it seem like it had nothing right (which would be more correct anyway).
            ## -- So truly it is undefined in that sense.
            ## With high enough coverage - one could always add pseudocounts to avoid this (like small priors)
            ##  Maybe that is a TODO
            ## -- I'd only add a pseudocount to the match though representing a small belief that the reference was assembled correctly
            ## -- In that sense it is just like counting the reference that was aligned to as one of the alignments
            ## For now I will look for the ZeroDivisionError and return '.' or NA
            self.counts = self._get_composition_counts(self.composition, self.ref, self.trantab, self.pattern, self.depth)
           
    def equal_or_above_threshold(self, threshold=1.0, stringency=1, strand=0):
        ''' stringency levels 1,2,3:
            3 means p(match) = 1 given X=DN (there can only be matches) - no mismatch, no del, no N
            2 means marg_over_N_p(match) = 1 -- no mismatches, no dels (Ns are ignored)
            1 means marg_over_N_D_p(match) = 1 -- no mismatches (Ns and dels are ignored)
            strand: 0=nonstranded, 1=plus,2=minus'''
##        if threshold > 1.0: #intended to print all
##            self.testvalue = -1.0
##            return False #Nothing above 1
        if stringency >= 3:
            self.testvalue = self.counts[strand]['p=']
        elif stringency < 3 and stringency >= 2:
            self.testvalue = self.counts[strand]['marg1_p=']
        else: #<=1
            self.testvalue = self.counts[strand]['marg2_p=']
        return self.testvalue >= threshold

    def perfect(self, stringency=1, strand=0):
        return self.equal_or_above_threshold(threshold=1.0, stringency=stringency, strand=strand)

    def passes_threshold(self, threshold=1.0, stringency=1, strand=0, above=False, reject_undefined=True):
        if self.equal_or_above_threshold(threshold, stringency, strand):
            passes = True if above else False
        else:
            passes = True if not above else False
        if self.testvalue is 'NA':
            passes = False if reject_undefined else True
        return passes

    def has_depth(self):
        return int(self.depth) > 0

    def has_min_depth(self, mindepth=0):
        return int(self.depth) > mindepth
    
    def __str__(self):
        if self.has_depth():
            probs = ['p=', 'pX', 'pD', 'pN', 'marg1_p=', 'marg1_pX', 'marg1_pD', 'marg2_p=', 'marg2_pX', 'pA', 'pC', 'pG', 'pT', 'marg1_pA',  'marg1_pC',  'marg1_pG',  'marg1_pT', 'marg2_pA', 'marg2_pC', 'marg2_pG', 'marg2_pT']
            counts = ['=', 'X', 'D', 'I', 'A', 'C', 'G', 'T', 'N']
            ends = ['rI', 'bI', '^', '$' ]
            symbols = counts + probs + ends
            nonstranded = [self.counts[0][k] for k in symbols]
            pos = [self.seq, self.pos, self.ref, self.depth]
            out = pos + nonstranded
            if self.stranded:
                plus = [self.counts[1][k] for k in symbols[:-4]]
                minus = plus = [self.counts[2][k] for k in symbols[:-4]]
                out += plus
                out += minus
        else:
            out = []
        return '\t'.join([str(e) for e in out])

    def _catch_zerodivisionerror(self, numerator, denominator):
        try:
            return numerator/denominator
        except ZeroDivisionError:
            return 'NA'
    def _get_probs(self, d):
        if sum(d.values()) > 0:
            #['p=', 'pX', 'pD', 'pN', 'marg1_p=', 'marg1_pX', 'marg1_pD', 'marg2_p=', 'marg2_pX', 'pA', 'pC', 'pG', 'pT', 'marg1_pA',  'marg1_pC',  'marg1_pG',  'marg1_pT', 'marg2_pA', 'marg2_pC', 'marg2_pG', 'marg2_pT']  
            ## Prob that a read has Match, Mismatch, or Deletion wrt reference (No need to include insertions here)
            marg2_total = float(d['='] + d['X']) #Only % match or mismatch
            marg1_total = float(marg2_total + d['D']) # % match or mismatch or del
            total = float(marg1_total + d['N']) # % match or mismatch or ambig
            d['p='] = self._catch_zerodivisionerror(float(d['=']), total)
            d['pX'] = self._catch_zerodivisionerror(float(d['X']), total)
            d['pD'] = self._catch_zerodivisionerror(float(d['D']), total)
            d['pN'] = self._catch_zerodivisionerror(float(d['N']), total)
            d['marg1_p='] = self._catch_zerodivisionerror(float(d['=']), marg1_total)
            d['marg1_pX'] = self._catch_zerodivisionerror(float(d['X']), marg1_total)
            d['marg1_pD'] = self._catch_zerodivisionerror(float(d['D']), marg1_total)
            d['marg2_p='] = self._catch_zerodivisionerror(float(d['=']), marg2_total)
            d['marg2_pX'] = self._catch_zerodivisionerror(float(d['X']), marg2_total)
            ## Prob that a read has A, C, G, T, N, or deletion over reference position
            ## Deletion prob already defined.
            new_marg2_total = float(d['A'] + d['C'] + d['G'] + d['T']) # Only % bases
            new_marg1_total = float(new_marg2_total + d['D']) # % bases or dels
            new_total = float(new_marg1_total + d['N']) # % bases or dels or ambig
            assert new_marg2_total == marg2_total and new_marg1_total == marg1_total and new_total == total
            d['pA'] = self._catch_zerodivisionerror(float(d['A']), new_total)
            d['pC'] = self._catch_zerodivisionerror(float(d['C']), new_total)
            d['pG'] = self._catch_zerodivisionerror(float(d['G']), new_total)
            d['pT'] = self._catch_zerodivisionerror(float(d['T']), new_total)
            ## pD should be same as above
            ## pN should be same as above
            d['marg1_pA'] = self._catch_zerodivisionerror(float(d['A']), new_marg1_total)
            d['marg1_pC'] = self._catch_zerodivisionerror(float(d['C']), new_marg1_total)
            d['marg1_pG'] = self._catch_zerodivisionerror(float(d['G']), new_marg1_total)
            d['marg1_pT'] = self._catch_zerodivisionerror(float(d['T']), new_marg1_total)
            ## marg1_pD should be same as above
            d['marg2_pA'] = self._catch_zerodivisionerror(float(d['A']), new_marg2_total)
            d['marg2_pC'] = self._catch_zerodivisionerror(float(d['C']), new_marg2_total)
            d['marg2_pG'] = self._catch_zerodivisionerror(float(d['G']), new_marg2_total)
            d['marg2_pT'] = self._catch_zerodivisionerror(float(d['T']), new_marg2_total)
        else:
            for e in ['p=', 'pX', 'pD', 'pN', 'marg1_p=', 'marg1_pX', 'marg1_pD', 'marg2_p=', 'marg2_pX', 'pA', 'pC', 'pG', 'pT', 'marg1_pA',  'marg1_pC',  'marg1_pG',  'marg1_pT', 'marg2_pA', 'marg2_pC', 'marg2_pG', 'marg2_pT']:
                d[e] = '.'

        return d

        
    def _get_composition_counts(self, composition, ref, trantab, pattern, depth):
        ''' Provide composition string.'''
        ref=ref.upper()
        symbols = '=XDIACGTN^$0' ## Right now '0' represents the 'deletion boundary' -[0-9]+[ACGTNacgtn]+ that occur before deletions (not on them).
        plus = {k:0 for k in symbols[:-2]}
        minus = {k:0 for k in symbols[:-2]}
        nonstranded = {k:0 for k in symbols}
        refrc = ref.translate(trantab)
        for pattern_match in re.findall(pattern, composition):
            #print pattern_match
            n = len(pattern_match)
            if pattern_match.startswith('.'): #base match Plus
                plus[ref] += n
                nonstranded[ref] += n
                plus['='] += n
                nonstranded['='] += n
            elif pattern_match.startswith(','): #base match Minus
                minus[refrc] += n
                nonstranded[ref] += n
                minus['='] += n
                nonstranded['='] += n
            elif pattern_match[0] in 'ACGTN': ## 1 or more mismatches Plus
    ##            plus['X'] += n
    ##            nonstranded['X'] += n
                for b in pattern_match:
                    if b != 'N': ## N in read will not be considered a mismatch, it will just be considered ambiguous
                        plus['X'] += 1
                        nonstranded['X'] += 1
                    plus[b] += 1
                    nonstranded[b] += 1
            elif pattern_match[0] in 'acgtn': ## 1 or more mismatches Minus -- but the lower-case base reflects the plus strand
    ##            minus['X'] += n
    ##            nonstranded['X'] += n
                for b in pattern_match.upper():
                    if b != 'N': ## N in read will not be considered a mismatch, it will just be considered ambiguous
                        minus['X'] += 1
                        nonstranded['X'] += 1
                    minus[b.upper().translate(trantab)] += 1 ## reflecting composition on reverse strd
                    nonstranded[b.upper()] += 1
            elif pattern_match.startswith('*'): ## deletions over position
                nonstranded['D'] += n
            elif pattern_match.startswith('-'): ## Deletion boundary -- analogous to how insertions, but not a del over pos
                nonstranded['0'] += 1
                if pattern_match[1] in 'ACGTN':
                    plus['0'] += 1
                elif pattern_match[1] in 'acgtn':
                    minus['0'] += 1
            elif pattern_match.startswith('+'):
                nonstranded['I'] += 1
                if pattern_match[1] in 'ACGTN':
                    plus['I'] += 1
                elif pattern_match[1] in 'acgtn':
                    minus['I'] += 1
            elif pattern_match.startswith('^'): ## read start
                nonstranded['^'] += 1 ## Only thing I did not have the regex grab more than 1 of... so just add 1
            elif pattern_match.startswith('$'): ## ends of read segments
                nonstranded['$'] += n

        ## Calculate proportions/probs from counts
        nonstranded = self._get_probs(nonstranded)
        plus = self._get_probs(plus)
        minus = self._get_probs(minus)

        ## Calculate insertion rate given depth and plus:minus insertion bias ratio
        nonstranded['rI'] = nonstranded['I']/float(depth)
        try:
            nonstranded['bI'] = float((plus['I']-minus['I']))/(plus['I']+minus['I'])
        except ZeroDivisionError:
            nonstranded['bI'] = 0.0
        return (nonstranded, plus, minus)
    
class Mpileup(object):
    def __init__(self, mpileup, stranded=False, header=False, ):
        self.stranded=stranded
        self.trantab = maketrans('ACGTacgtNn', 'TGCAtgcaNn')
        ##pattern = re.compile('\.|,|\^.|\$|[\+|-][0-9]+[A|C|G|T|N|a|c|g|t|n]+')
        ##pattern = re.compile('\.|,|\^.|\$|[\+|-]?[0-9]*[A|C|G|T|N|a|c|g|t|n]+')
        #self.pattern = re.compile('\.+|,+|\^.|\$+|[\+|-]?[0-9]*[A|C|G|T|N]+|[\+|-]?[0-9]*[a|c|g|t|n]+')
        ## Here are important patterns to look for -- 
        fwd_matches_on_pos = '\.+'
        rev_matches_on_pos = ',+'
        readstart_w_mapq_on_pos = '\^.'
        readends_on_pos = '\$+'
        dels_on_pos = '\*+'
        fwd_ins_after_pos = '\+[0-9]*[A|C|G|T|N]+'
        rev_ins_after_pos = '\+[0-9]*[a|c|g|t|n]+'
        fwd_del_after_pos = '\-[0-9]*[A|C|G|T|N]+' ## Not doing much with this for now
        rev_del_after_pos = '\-[0-9]*[a|c|g|t|n]+'
        fwd_mismatches_on_pos = '[A|C|G|T|N]+' ## patterns are evaluated left to right - tests show that it doesnt mater if these come before or after indels (possible b/c it prefers the greediest one), but it seems safer to keep them after (or combined as below)
        rev_mismatches_on_pos = '[a|c|g|t|n]+'
        fwd_insafter_delafter_or_mms = '[\+|-]?[0-9]*[A|C|G|T|N]+'
        rev_insafter_delafter_or_mms = '[\+|-]?[0-9]*[a|c|g|t|n]+'
        self.patterns = '|'.join( [fwd_matches_on_pos, rev_matches_on_pos, dels_on_pos, readstart_w_mapq_on_pos, readends_on_pos, fwd_insafter_delafter_or_mms, rev_insafter_delafter_or_mms] )
        self.pattern = re.compile(self.patterns)
        if mpileup in ('-', 'stdin'):
            self.mpileup = sys.stdin
            self.need_to_close_file = False
        else:
            self.mpileup = open(mpileup)
            self.need_to_close_file = True
    def __iter__(self):
        return self
    def next(self):
        nextrec = self.mpileup.readline()
##        try:
##            return MpileupRecord(nextrec, self.trantab, self.pattern, self.stranded)
##        except ValueError, e:
##            raise StopIteration
        if nextrec:
            return MpileupRecord(nextrec, self.trantab, self.pattern, self.stranded)
        else:
            raise StopIteration

    def close(self):
        if self.need_to_close_file:
            self.mpileup.close()
    def header(self):
        pos = ['chr', 'pos', 'ref', 'depth']
        probs = ['p=', 'pX', 'pD', 'pN', 'marg1_p=', 'marg1_pX', 'marg1_pD', 'marg2_p=', 'marg2_pX', 'pA', 'pC', 'pG', 'pT', 'marg1_pA',  'marg1_pC',  'marg1_pG',  'marg1_pT', 'marg2_pA', 'marg2_pC', 'marg2_pG', 'marg2_pT']
        counts = ['=', 'X', 'D', 'I', 'A', 'C', 'G', 'T', 'N']
        ends = ['rI', 'bI', '^', '$' ]
        nonstranded = counts + probs + ends
        preout = pos + nonstranded
        if self.stranded:
            plus = ['plus_'+e for e in nonstranded[:-4]]
            minus = ['minus_'+e for e in nonstranded[:-4]]
            preout += plus + minus
        out = [str(i+1)+'='+preout[i] for i in range(len(preout))]
        return '\t'.join(out)

 
