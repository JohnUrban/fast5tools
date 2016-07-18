## JOHN URBAN (2015, 2016)

#info
import h5py, os, sys, tarfile, shutil
##import cStringIO as StringIO
from Bio import SeqIO
from glob import glob
import datetime
from random import randint


##class FastXMolecule(object):
##    ## wraps over fasta/fastq files that may have 1 or many contents
##    ## May want a separate class that records are given to
##    ## records will be BioPython in nature....
##    def __init__(self, filename):
##        self.filename = filename
##        self.filebasename = (".").join(filename.split("/")[-1].split(".")[:-1])
##        self.abspath = os.path.abspath(filename)
##        self.is_open = self.open() #self.f5 is the f5 file object
##        if self.is_open:# and self.is_nonempty():
##            self._has_read["template"] = self.has_template()
##            self._has_read["complement"] = self.has_complement()
##            self._has_read["2d"] = self.has_2d()
##            self._has_read["input"] = True # arbitrarily true for now
##            self.molecule = None
##            self.molequal = None
##            self.log = None
##            self.fq_path = None
##            self.given_name = {"template":None, "complement":None, "2d":None}
##            self.seq = {"template":None, "complement":None, "2d":None}
##            self.fq_sep = {"template":None, "complement":None, "2d":None}
##            self.quals = {"template":None, "complement":None, "2d":None}
##            self.quals_as_int = {"template":None, "complement":None, "2d":None}
##            self.info_name = {"template":None, "complement":None, "2d":None}
##            self.base_info_name = None
##            self.GENERAL_PATH = None
##            self._get_general_path()
##
##    def __str__(self):
##        return self.filename
##    
##    def open(self):
##        """
##        Open an ONT Fast5 file, assuming HDF5 format
##        """
##        try:
####            self.fx = SeqIO.(self.filename, 'r') 
##            return True
##        except Exception, e:
####            sys.stderr.write("Cannot open file: %s \n" % self.filename)
####            logger.warning("Cannot open file: %s" % self.filename)
##            return False
##
##    def close(self):
##        """
##        Close an open an ONT Fast5 file, assuming HDF5 format
##        """
##        if self.is_open:
##            self.f5.close()
##
##    def is_template(self):
##        pass
##
##    def is_complement(self):
##        pass
##
##    def is_2d(self):
##        pass
##
##    def _define_molecule(self):
##        ## if 2d, use 2d as molecule
##        ## if both t and c present, use the longer one
##        ## if only t present, use t
##        ### NOTE: as is - when fast5 is "empty" (lacking even template), this gives answer as template
##        ### i.e. assumes non-empty
##        if self.has_read("2d"):
##            self.molecule = "2d"
##        elif self.has_read("complement"):
##            if int(self.get_seq_len("complement")) > int(self.get_seq_len("template")):
##                self.molecule = "complement"
##            else:
##                self.molecule = "template"
##        else:
##            self.molecule = "template"
##
##    def _define_molequal(self): ## added July 6, 2016 -- copied/modified from define_molecule
##        ## if 2d, use 2d as molecule
##        ## **if both t and c present, use the one of higher quality** <- how it differs from molecule
##        ## if only t present, use t
##        ### NOTE: as is - when fast5 is "empty" (lacking even template), this gives answer as template
##        ### i.e. assumes non-empty
##        if self.has_read("2d"):
##            self.molequal = "2d"
##        elif self.has_read("complement"):
##            if int(self.get_mean_qscore("complement")) > int(self.get_mean_qscore("template")):
##                self.molequal = "complement"
##            else:
##                self.molequal = "template"
##        else:
##            self.molequal = "template"
##
##
##    def use_molecule(self):
##        if self.molecule == None:
##            self._define_molecule()
##        return self.molecule
##
##    def use_molequal(self):
##        if self.molequal == None:
##            self._define_molequal()
##        return self.molequal
##
##    def get_mean_qscore(self, readtype):
##        '''readtpye in 2d, template, complement'''
##        return self._get_attr(path = self._get_attr_path(readtype), attr = "mean_qscore")
##
##
##    def get_seq_len(self, readtype):
##        '''readtpye in 2d, template, complement'''
##        return self._get_attr(path = self._get_attr_path(readtype), attr = "sequence_length")
##
##
##    def get_device_id(self):
##        return self.f5[TRACKING_ID].attrs["device_id"]
##
##    def get_asic_id(self):
##        return self.f5[TRACKING_ID].attrs["asic_id"]
##
##
##    def get_run_id(self):
##        return self.f5[TRACKING_ID].attrs["run_id"]
##
##
##
##    def get_file_version(self):
##        return self.f5['/'].attrs['file_version']
##
##
##    def get_read_number(self): ## June 22, 2016 -- not tested for all yet -- developing for non-basecalled files
##        return self.f5["/Analyses/EventDetection_000/Reads/"].keys()[0]
##
##    def get_model_type(self):
##        return self.f5[self.GENERAL_PATH].attrs['model_type']
##
##
##
##    def _parse_fastq_info(self, readtype):
##        if self.seq[readtype] == None:
##            try:
##                (self.given_name[readtype], self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]) = self.f5[self._get_location_path(readtype,"Fastq")][()].strip().split('\n')
##                self.given_name[readtype] = self.given_name[readtype].lstrip('@')
##            except:
##                (self.given_name[readtype], self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]) = None, None, None, None
##
##
##    def _get_pore_info_name(self, readtype):
##        ## TODO -- this will store most read-pertinent info
##        if self.base_info_name == None:
##            self._get_base_info_name()
##        if self.info_name[readtype] == None:
##            info = []
##            info.append(readtype)
##            info.append("len_"+str(self.get_seq_len(readtype)))
##            info.append("Q_"+str(self.get_mean_qscore(readtype)))
##            self.info_name[readtype] = ("_").join(info)
##        return self.info_name[readtype] + "_" + self.base_info_name
##    
##    def _get_base_info_name(self):
##        info = []
##        info.append("channel_"+self.get_channel_number())
##        info.append(self.get_read_number())
####        info.append("file_"+self.get_file_number())
##        info.append("asic_"+self.get_asic_id())
##        info.append("run_"+self.get_run_id()) ## Though ASIC ID is different for different, if someone re-uses a flowcell, then it will not be. But I believe the RunID will be.... unless it is given the same tag dependent on ASIC ID and/or other constants
##        info.append("device_"+self.get_device_id())
##        info.append("model_"+self.get_model_type())
####            info.append(("-").join(self.get_time_stamp().split()))
##        #could add start time, duration to help distinguish...
##        self.base_info_name = ("_").join(info)
##        
##    def get_base_info_name(self):
##        if self.base_info_name == None:
##            self._get_base_info_name()
##        return self.base_info_name
##
##    
##    def get_fastq(self, readtype):
##        if self.has_read(readtype):
##            self._parse_fastq_info(readtype)
##            return '\n'.join(['@'+self._get_pore_info_name(readtype), self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]])
##
##    def get_fasta(self, readtype):
##        if self.has_read(readtype):
##            self._parse_fastq_info(readtype)
##            return '\n'.join(['>'+self._get_pore_info_name(readtype), self.seq[readtype]])
##
##    def get_quals(self, readtype):
##        if self.has_read(readtype):
##            self._parse_fastq_info(readtype)
##            return '\n'.join(['>'+self._get_pore_info_name(readtype), self.quals[readtype]])
##
##    def get_quals_as_int(self, readtype):
##        if self.has_read(readtype):
##            self._parse_fastq_info(readtype)
##            if self.quals_as_int[readtype] == None:
##                self.quals_as_int[readtype] = StringIO.StringIO()
##                SeqIO.convert(StringIO.StringIO('\n'.join(['@'+self._get_pore_info_name(readtype), self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]])), "fastq", self.quals_as_int[readtype], "qual")
##            return self.quals_as_int[readtype].getvalue().rstrip()
##
##    def get_start_time(self, which):
##        pass
##        ## may want to add start times to fastx reads to filter by relative times
##
##    def get_time_stamp(self):
##        #perhaps in future
##        pass
##
##    def get_channel_number(self):
##        pass
##    
##    def get_heatsink_temp(self):
##        #perhaps in future
##        pass
##    
##    def get_asic_temp(self):
##        #perhaps in future
##        pass
##
##    def get_basecall_version(self):
##        pass
##        ## perhaps in future
##


class FastXMolecule(object):
    def __init__(self, fx):
        self.reads = {}
        self.reads['template'] = None
        self.reads['2d'] = None
        self.reads['complement'] = None
        self._has = {}
        self._has['2d'] = False
        self._has['complement'] = False
        self._has['template'] = False
        self.molecule_name = None
        self.add_fx(fx)

    def is_empty(self):
        return self.reads['template'] == None and self.reads['complement'] == None and self.reads['2d'] == None and self.molecule_name == None

    def is_insertable(self, fx):
        if fx.get_molecule_name() == fx.get_molecule_name():
            if self.reads[fx.get_read_type()] == None:
                return True
        return False

    def add_fx(self, fx):
        ##
        if self.is_empty():
            self.molecule_name = fx.get_molecule_name()
        if self.is_insertable(fx):
            self.reads[fx.get_read_type()] = fx
            self._has[fx.get_read_type()] = True

    def get_molecule_name(self):
        return self.molecule_name

    def get_molecule_name_as_str(self, delim="_"):
        return (delim).join([str(e) for e in self.molecule_name])
    
    def has_template(self):
        return self._has["template"]

    def has_complement(self):
        return self._has["complement"]
    
    def has_2d(self):
        return self._has["2d"]

    def template(self):
        return self.reads["template"]

    def complement(self):
        return self.reads["complement"]

    def twod(self):
        return self.reads["2d"]

    def get_read(self, readtype):
        return self.reads[readtype]

    def molecule(self):
        if self.has_2d():
            return self.reads["2d"]
        elif self.has_complement():
            assert self.has_template()
            if self.template().get_seq_len() > self.complement().get_seq_len():
                return self.template()
            else:
                return self.complement()
        return self.template()

    def molequal(self):
        if self.has_2d():
            return self.reads["2d"]
        elif self.has_complement():
            assert self.has_template()
            if self.template().get_mean_qscore() > self.complement().get_mean_qscore():
                return self.template()
            else:
                return self.complement()
        return self.template()
            
    def passes_filter(self, readtype, minlen=0, maxlen=float('inf'), minq=0, maxq=float('inf'), channel=None, readnum=None, asic=None, runid=None, deviceid=None, modelid=None, minscore=1, filter_rule="and"):
        if self._has[readtype]:
            return self.get_read(readtype).passes_filter(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid)
        return False
        
    def interpret(self, readtype):
        if readtype == "molecule":
            return self.molecule().get_read_type()
        elif readtype == "MoleQual":
            return self.molequal().get_read_type()
        else:
            return readtype
    
    def _get_read_type_list(self, readtype):
        if readtype == "molecule":
            return [self.molecule().get_read_type()]
        elif readtype == "MoleQual":
            return [self.molequal().get_read_type()]
        if readtype == "all":
            return ["template", "complement", "2d"]
        else:
            return [readtype]
        
    def get_fastx_entry(self, readtype, outtype, i):
        ## i is for falconizing name
        for rtype in self._get_read_type_list(readtype):
            if self._has[readtype]:
##                    if self.passes_filter(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid):
                return self.get_read(readtype).get_fastx_entry(outtype, i)
##                    self.get_read(readtype).passes_filter(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid)
        


    def get_molecule_stats(self):
        molecule_stats = []
        molecule_stats.append(self.get_molecule_name_as_str())
        molecule_stats.append(self.molecule().get_seq_len())
        molecule_stats.append(int(self.has_complement()))
        molecule_stats.append(int(self.has_2d()))
        molecule_stats.append(self.twod().get_seq_len()) if self.has_2d() else molecule_stats.append("-")
        molecule_stats.append(self.template().get_seq_len()) if self.has_template() else molecule_stats.append("-")
        molecule_stats.append(self.complement().get_seq_len()) if self.has_complement() else molecule_stats.append("-")
        molecule_stats.append(self.twod().get_mean_qscore()) if self.has_2d() else molecule_stats.append("-")
        molecule_stats.append(self.template().get_mean_qscore()) if self.has_template() else molecule_stats.append("-")
        molecule_stats.append(self.complement().get_mean_qscore()) if self.has_complement() else molecule_stats.append("-")
        return molecule_stats 

    def get_molecule_stats_string(self, delim="\t"):
        return (delim).join([str(e) for e in self.get_molecule_stats()])
                                                                                                                   

##        if filter_rule == "and":
##            return self.is_in_intersect(minlen, maxlen, minq, maxq, readtype, channel, readnum, asic, run, device, model)
##        elif filter_rule == "or":
##            return self.is_in_union(minlen, maxlen, minq, maxq, readtype, channel, readnum, asic, run, device, model)
##        return None
##
##    def is_in_intersect(self, minlen, maxlen, minq, maxq, readtype, channel, readnum, asic, run, device, model):
##        pass
##
##    def is_in_union(self, minlen, maxlen, minq, maxq, readtype, channel, readnum, asic, run, device, model):
##        pass
##
##    def get_reads(self, readtype):
##        pass




class FastXSeqFromFast5(object):
    def __init__(self, fxseqio):
        #fxseqio is seqio object from SeqIO.parse (biopython)
        self._is_2d = None
        self._is_template = None
        self._is_complement = None
        self.fx = fxseqio
        self.parse_fx_name()
        self._len_correct = None
        self.qualschecked = False
        
    def parse_fx_name(self):
        info = self.fx.name.split("|")
        self.info = {'readtype':info[0]}
        for i in range(1,len(info)):
            key, value = info[i].split(":")
            self.info[key] = value

    def determine_read_type(self):
        if self.info['readtype'] in ('2d','2D'):
            self._is_2d = True
        else:
            self._is_2d = False
        if self.info['readtype'].lower() == 'template':
            self._is_template = True
        else:
            self._is_template = False
        if self.info['readtype'] == 'complement':
            self._is_complement = True
        else:
            self._is_complement = False
        assert sum([self._is_2d, self._is_template, self._is_complement]) == 1

    def get_name(self):
        return self.fx.name

    def get_description(self):
        return self.fx.description

    def get_molecule_name(self):
        ## return tuple, string, ..? right now its a tuple.
        return (self.get_asic_id(), self.get_run_id(), self.get_channel(), self.get_read_number())
    
    def get_seq(self):
        return str(self.fx.seq)

    def get_reverse_complement(self):
        return str(self.fx.seq.reverse_complement())

    def get_run_id(self):
        return self.info['run']

    def get_device_id(self):
        return self.info['device']

    def get_asic_id(self):
        return self.info['asic']

    def get_seq_len(self, as_int=True):
        if as_int:
            return int(self.info['len'])
        return self.info['len']

    def get_mean_qscore(self, as_float=True):
        if as_float:
            return float(self.info['Q'])
        return self.info['Q']

    def get_channel(self):
        return self.info['channel']
    
    def get_read_number(self):
        return self.info['Read']

    def get_read_type(self):
        return self.info['readtype']

    def get_model_type(self):
        return self.info['model']

    def is_readtype(self, readtype):
        return self.info['readtype'] == readtype

    def is_length(self, minlen=0, maxlen=float('inf')):
        return self.get_seq_len() >= minlen and self.get_seq_len() <= maxlen

    def is_quality(self, minq=0, maxq=float('inf')):
        return self.get_mean_qscore() >= minq and self.get_mean_qscore() <= maxq

    def is_from_channel(self, channel):
        return self.info['channel'] == channel

    def is_read_number(self, readnum):
        return self.info['Read'] == readnum

    def is_from_asic(self, asic):
        return self.info['asic'] == asic

    def is_from_run(self, runid):
        return self.info['run'] == runid

    def is_from_device(self, deviceid):
        return self.info['device'] == deviceid

    def basecall_model_was(self, modelid):
        return self.info['model'] == modelid

    def full_comparison(self, readtype, minlen=0, maxlen=float('inf'), minq=0, maxq=float('inf'), channel=None, readnum=None, asic=None, runid=None, deviceid=None, modelid=None):
        t1 = self.is_readtype(readtype)
        t2 = self.is_length(minlen,maxlen)
        t3 = self.is_quality(minq,maxq)
        t4 = self.is_from_channel(channel) or channel == None
        t5 = self.is_read_number(readnum) or readnum == None
        t6 = self.is_from_asic(asic) or asic == None
        t7 = self.is_from_run(runid) or runid == None
        t8 = self.is_from_device(deviceid) or deviceid == None
        t9 = self.basecall_model_was(modelid) or modelid == None
        return [t1,t2,t3,t4,t5,t6,t7,t8,t9]

    def full_comparison_score(self, readtype, minlen=0, maxlen=float('inf'), minq=0, maxq=float('inf'), channel=None, readnum=None, asic=None, runid=None, deviceid=None, modelid=None):
        fclist = self.full_comparison(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid)
        score = sum(fclist)
        topscore = len(fclist)
        return score, topscore

    def passes_intersect_of(self, readtype, minlen=0, maxlen=float('inf'), minq=0, maxq=float('inf'), channel=None, readnum=None, asic=None, runid=None, deviceid=None, modelid=None):
        score, topscore = self.full_comparison_score(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid)
        return score == topscore

    def passes_union_of(self, readtype, minlen=0, maxlen=float('inf'), minq=0, maxq=float('inf'), channel=None, readnum=None, asic=None, runid=None, deviceid=None, modelid=None, minscore=1):
        score, topscore = self.full_comparison_score(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid)
        return score >= minscore

    def passes_filter(self, readtype, minlen=0, maxlen=float('inf'), minq=0, maxq=float('inf'), channel=None, readnum=None, asic=None, runid=None, deviceid=None, modelid=None, minscore=1, filter_rule="and"):
        if filter_rule == "and":
            return self.passes_intersect_of(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid)
        elif filter_rule == "or":
            return self.passes_union_of(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore)

    def is_2d(self):
        if self._is_2d == None:
            self.determine_read_type()
        return self._is_2d

    def is_template(self):
        if self._is_template == None:
            self.determine_read_type()
        return self._is_2d

    def is_complement(self):
        if self._is_complement == None:
            self.determine_read_type()
        return self._is_2d

    def stored_seq_len_matches_fxlen(self):
        if self.fxlen == None:
            self.fxlen = len(self.fx.seq)
        if self._len_correct == None:
            self._len_correct = self.info['len'] == self.fxlen
        else:
            return self._len_correct

    def fix_stored_seq_len(self):
        ## This should never have to be used. If stored Seq len doesnt match, there is a problem somewhere in pipeline.
        if not stored_seq_len_matches_fxlen():
            self.info['len'] = self.fxlen

    def _check_quals(self):
        try:
            self.fx.letter_annotations['phred_quality']
            self._quals_are_floats = False
        except KeyError: ## if seq was fasta record
            self._quals_are_floats = True
            self.fx.letter_annotations['phred_quality'] = [int(round(self.get_mean_qscore()))]*self.get_seq_len()
                    ## if started out as fasta, when quals need conversion will have to be rounded to nearest int
            #
        self.qualschecked = True


    def int_quals(self):
        if not self.qualschecked:
            self._check_quals()
        return self.fx.letter_annotations['phred_quality']

    def sanger_quals(self):
        if not self.qualschecked:
            self._check_quals()
        return SeqIO.QualityIO._get_sanger_quality_str(self.fx)


        
    def get_fastq(self):
        if not self.qualschecked:
            self._check_quals()
        return '\n'.join(['@'+self.fx.description, self.get_seq(), '+', self.sanger_quals()])

    def get_fasta(self):
        return '\n'.join(['>'+self.fx.description, self.get_seq()])


    def falconize_name(self, zmw_num, style="old"):
        if style == "old":
            moviename = "m000_000"
            otherinfo = self.fx.description
        elif style == "new":
            moviename = "asic:"+self.get_asic_id() + "|run:"+self.get_run_id() + "|device:"+self.get_device_id() + "|model:"+self.get_model_type()
            otherinfo = self.get_read_type() + "|Q:"+self.info['Q'] + "|Read:"+self.info['Read'] + "|channel:"+self.info['channel']
        return moviename + "/" + str(zmw_num) + "/0_"+self.info['len'] + " " + otherinfo
    
    def get_falcon_fasta(self, zmw_num=None, style="old"):
        if zmw_num == None:
            zmw_num = randint(0,1000000000000)
        return '\n'.join([">"+self.falconize_name(zmw_num, style), self.get_seq()])


    def get_sanger_quals(self):
        if not self.qualschecked:
            self._check_quals()
        return '\n'.join(['>'+self.fx.description, self.sanger_quals()])

    def get_int_quals(self):
        if not self.qualschecked:
            self._check_quals()
        return '\n'.join(['>'+self.fx.description, (" ").join([str(e) for e in self.int_quals()])])


    def get_fastx_entry(self, outtype, i):
        if outtype == "fasta":
            return self.get_fasta()
        elif outtype == "fastq":
            return self.get_fastq()
        elif outtype == "qual":
            return self.get_sanger_quals()
        elif outtype == "intqual":
            return self.get_int_quals()
        elif outtype == "falcon":
            return self.get_falcon_fasta(i, "old")
        elif outtype == "oldfalcon":
            return self.get_falcon_fasta(i, "old")
        elif outtype == "newfalcon":
            return self.get_falcon_fasta(i, "new")


    def __str__(self):
        return str(self.info)

    
class FastXFile(object):
    def __init__(self, fxname, intype=None):
        if intype == None or intype == "input":
            self.intype = self.determine_type(fxname)
        else:
            self.intype = intype
        self.fx = SeqIO.parse(fxname, self.intype)
        self.filename = fxname
        
    def parse_as_seqio(self):
        return self.fx

    def determine_type(self, fxname):
        for sfx in ['.fasta', '.fa', '.fna']:
            if fxname.endswith(sfx):
                return 'fasta'
        for sfx in ['.fastq', '.fq']:
            if fxname.endswith(sfx):
                return 'fastq'
        ## qual and intqual may need to be added.
        return None

    def __iter__(self):
        return self

    def next(self):
        try:
            seqio_fx = self.fx.next()
            fx = FastXSeqFromFast5(seqio_fx)
            return fx
        
        except Exception as e:
            raise StopIteration
    

## May need to add random seeds to avoid cross-talk between different operations in same dir
F5_TMP_DIR = ".fast5tools_tmp_dir"
F5_TMP_DIR = "fast5tools_tmp_dir"



class FastXFileList(object):
    def __init__(self, fxlist, intypes=['fasta'], tar_filenames_only=False, keep_tar_footprint_small=False):
        #ensure type is list
        if isinstance(fxlist, list):
                self.fxlist = fxlist
        elif isinstance(filelist, str):
                self.fxlist = [fxlist]
        self._define_allowable_files(intypes)
        self._tars_detected = False
        self.tar_filenames_only = tar_filenames_only
        self.keep_tar_footprint_small = keep_tar_footprint_small
        self.nfiles = None
        self.allfiles = None
        self.n_tars = 0
        self.tars = {} ## for small footprint method
##        self._expand_list() ##TODO(?) - make expand fofn optional, so not wasting time in most situations
        self._extract_fastx_files()

        
    def _define_allowable_files(self, intypes):
        self.intypes = []
        if 'fasta' in intypes:
            self.intypes.append('.fasta')
            self.intypes.append('.fa')
            self.intypes.append('.fna')
        if 'fastq' in intypes:
            self.intypes.append('.fastq')
            self.intypes.append('.fq')

        
    def __iter__(self):
        return self

    def next(self):
        try:
            newfile = self.allfiles.next()
            if self.keep_tar_footprint_small:
                if newfile.startswith("f5tar|"):
                    f5tar, key, tar_member = newfile.split("|")
                    tarkey = "f5tar|" + key + "|"
                    self.tars[tarkey].extract(tar_member, path=F5_TMP_DIR)
                    newfile = os.path.join(F5_TMP_DIR, tar_member)
                    fx = FastXFile(newfile)
##                    fx = newfile
                    os.remove(newfile)
            else:
                fx = FastXFile(newfile)
##                fx = newfile
            return fx
        
        except Exception as e:
            if self._tars_detected and os.path.exists(F5_TMP_DIR):
                shutil.rmtree(F5_TMP_DIR)
            raise StopIteration
        
    
    def _initialize_tar_tmp_dir(self):
        if os.path.isdir(F5_TMP_DIR):
            shutil.rmtree(F5_TMP_DIR)
        self._tars_detected = True
        os.mkdir(F5_TMP_DIR)

    def allowable_file(self, fname):
        #fname is string filename
        for sfx in self.intypes:
            if fname.endswith(sfx):
                return True
        return False
    
    def _expand_fofn(self, fofn):
        # FOFN must have single file per line - where a file can be of 3 types:
        # allows for .fast5 files to be listed
        # allows for directories to be listed
        # allows for fofn to include tars
        # FOFNs, for now, should not include other FOFNs - an fofn listng itself could create an infinite loop
        f = open(fofn,'r')
        lines = [line.strip() for line in f.readlines()]
        f.close()
        files = []
        for fname in lines:
            if self.allowable_file(fname):
                files.append(fname)
            elif os.path.isdir(fname):
                files += self._expand_dir(fname)
            elif tarfile.is_tarfile(fname):
                files += self._expand_tar(fname)
            else: ## line in FOFN is ignored
                pass
        return files

    def _expand_dir(self, d):
        files = []
        for sfx in self.intypes:
            pattern = d + '/' + '*' + sfx
            files.append(glob(pattern))
        return files

    def _expand_tar(self, tarball):
        if not self._tars_detected and not self.tar_filenames_only:
            self._initialize_tar_tmp_dir()
        f = tarfile.open(tarball)
        self.n_tars += 1
        if self.tar_filenames_only:
            files = [os.path.join(tarball,fname) for fname in f.getnames() if self.allowable_file(fname)]
            f.close()
        else: ## will be using files in tarball
            if self.keep_tar_footprint_small:
                tarkey = "f5tar|"+str(self.n_tars)+"|"
                self.tars[tarkey] = f
                files = [tarkey+fname for fname in f.getnames() if self.allowable_file(fname)]
                ## purposely do not close tarfile
            else:                  
                f.extractall(path=F5_TMP_DIR) ## in situations where tarball includes many big non-fast5 files, this may not be best way.
                ## os.path.basename(filename).endswith('.fast5') and not os.path.basename(filename).startswith('.')
                files = [os.path.join(F5_TMP_DIR, fname) for fname in f.getnames() if self.allowable_file(fname)]
                f.close()
        return files

    def _extract_fastx_files(self):
        # can be
        # [f.fast5] -- list of a single file
        # [f5dir/] -- list of a single dir
        # [f1.fast5, ..., fN.fast5] -- list of a multiple files
        # [f5dir1, ..., f5dirN/] -- list of multiple dirs
        # [f1.fast5, f5dir1/, f5dir2/, f2.fast5, ..., fN.fast5, f5dirM/] -- mixed list of N files and M dirs
        # FOFNs can be thrown in as well
        self.files = []
        for e in self.fxlist:
            if self.allowable_file(e):
                self.files.append(e)
            elif e.endswith(".fofn"):
                self.files += self._expand_fofn(e)
            elif os.path.isdir(e):
                self.files += self._expand_dir(e)
            elif tarfile.is_tarfile(e):
                self.files += self._expand_tar(e)
            else: ## all else currently ignored by fast5tools
                pass
        self.nfiles = len(self.files)
        self.allfiles = iter(self.files)

    def __len__(self):
        return self.nfiles

    def __eq__(self, other):
        return set(self.files) == set(other.files)

    def get_filenames(self):
        return self.files

    def get_basenames(self):
        return [os.path.basename(e) for e in self.files]

    def get_dirnames(self):
        return [os.path.dirname(e) for e in self.files]











