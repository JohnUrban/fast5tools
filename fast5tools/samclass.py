import os, sys, re
from collections import defaultdict

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


    def update_pos_field(self, add=0, replace=None):
        assert add == 0 or replace is None
        if add != 0:
            self.parsed_aln['POS'] += add
        elif replace is not None:
            self.parsed_aln['POS'] = replace

    def find_given_extra_field(self, regex):
        '''Returns only first match. There should only be one match.'''
        for e in self.parsed_aln['EXTRA'].split('\t'):
            found = re.search(regex, e)
            if found is not None:
                break
        return found.group(1)


    def get_edit_dist_field(self):
        '''Assumes NM:i: field is present for now...'''
        if self.edit_dist is None: 
            self.edit_dist = int(self.find_given_extra_field(regex='NM:i:(\d+)'))
        return self.edit_dist


    def get_fast5_field(self):
        '''Only works with SAM files that have the F5:Z: tag from fast5tools'''
        if self.fast5info is None: 
            self.fast5info = self.find_given_extra_field(regex='F5:Z:(.*)')
        return self.fast5info
            


    def get_clipping_dist(self):
        if self.cigar is None: 
            self.get_cigar_counts()
        if self.clipping_dist is None: 
            self.clipping_dist = sum([self.cigar[e] for e in 'HS'])
        return self.clipping_dist
            

    def get_cigar_list(self):
        return [(int(count),char) for count,char in re.findall('(\d+)([MIDNSHP=X])', self.get_cigar_field())]

    def convert_cigar_element_

    def get_edit_dist_with_clipping(self):
        return int(self.get_edit_dist_field()) + int(self.get_clipping_dist())

    def get_reference_aln_len(self):
        ''' This gives length of alignment from POV of reference according to pos field and cigar.
            Add up matches to ref and deletions from ref to get ref len in aln.
            The reference end coord is Pos + MD - 1'''
        return sum([self.cigar[e] for e in 'MD'])

    def get_reference_end_pos(self):
        ''' This gives tuple of 1-based (as SAM is) start/end on reference according to pos field and cigar.
            The reference end coord is Pos + MD - 1 (i.e. add up matches to ref and deletions to ref to get ref len in aln)'''
        ## Note: obtain "start pos" with self.get_pos_field()
        return self.get_pos_field() + self.get_reference_aln_len() - 1


    def get_adjusted_pos_field(self, adjust_start=True):
        '''If adjust_front is False, it adjusts end pos.'''
        if adjust_start:
            initialpos = self.get_pos_field()
            cigaridx1=0
            cigaridx2=1
            direction=-1
        else:
            initialpos = self.get_reference_end_pos()
            cigaridx1=-1
            cigaridx2=-2
            direction=1
        clip = 0
        if self.get_clipping_dist() > 0:
            cigar_list = self.get_cigar_list()
            if cigar_list[cigaridx1][1] in ('H', 'S'):
                clip += cigar_list[cigaridx1][0]
                if cigar_list[cigaridx2][1] == 'S': #possible if H was first
                    clip += cigar_list[cigaridx2][0]
        return initialpos + direction*clip  

    def get_clipping_adjusted_start(self):
        return self.get_adjusted_pos_field(cigaridx1=0, cigaridx2=1, direction=-1)

    def get_clipping_adjusted_end(self):
        return self.get_adjusted_pos_field(cigaridx1=-1, cigaridx2=-2, direction=1)
        

    def get_proportion_of_read_len(self, proportion=0.1):
        '''Returns int'''
        return int( self.get_read_len() * proportion )









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


