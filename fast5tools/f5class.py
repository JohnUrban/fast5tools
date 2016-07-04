## JOHN URBAN (2015, 2016)

#info
import h5py, os, sys, tarfile, shutil
import cStringIO as StringIO
from Bio import SeqIO
from glob import glob

#logging
import logging
logger = logging.getLogger('fast5tools')
logging.basicConfig()

LOC_TEMP = "/Analyses/Basecall_2D_000/BaseCalled_template/" #01
LOC_TEMP2 = "/Analyses/Basecall_1D_000/BaseCalled_template/"
LOC_COMP = "/Analyses/Basecall_2D_000/BaseCalled_complement/" #01
LOC_COMP2 = "/Analyses/Basecall_1D_000/BaseCalled_complement/"
LOC_2D = "/Analyses/Basecall_2D_000/BaseCalled_2D/" #01
INPUT_EVENTS_PATH = "Analyses/EventDetection_000/Reads/" #01


ATTR_TEMP = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/" #01
ATTR_TEMP2 = "/Analyses/Basecall_1D_000/Summary/basecall_1d_template/"
ATTR_COMP = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/" #01
ATTR_COMP2 = "/Analyses/Basecall_1D_000/Summary/basecall_1d_complement/"
ATTR_2D = "/Analyses/Basecall_2D_000/Summary/basecall_2d" #01

BASECALL_TEST = '/Analyses/Basecall_2D_000' #01

GENERAL = "/Analyses/Basecall_2D_000/Configuration/general"
GENERAL001 = "/Analyses/Basecall_1D_000/Configuration/general"

TRACKING_ID = "/UniqueGlobalKey/tracking_id"
SPLIT_HAIRPIN = "/Analyses/Basecall_2D_000/Summary/split_hairpin/" #01
SPLIT_HAIRPIN2 = "/Analyses/Hairpin_Split_000/Summary/split_hairpin/"


ALIGNMENT_2D = "/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment" #01
FASTQ_2D = "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq" #01
EVENTS_COMP = "/Analyses/Basecall_2D_000/BaseCalled_complement/Events" #01
FASTQ_COMP = "/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq" #01
MODEL_COMP = "/Analyses/Basecall_2D_000/BaseCalled_complement/Model" #01
EVENTS_TEMP = "/Analyses/Basecall_2D_000/BaseCalled_template/Events"  #01
FASTQ_TEMP = "/Analyses/Basecall_2D_000/BaseCalled_template/Fastq" #01
MODEL_TEMP = "/Analyses/Basecall_2D_000/BaseCalled_template/Model" #01
ALIGNMENT_HAIRPIN = "/Analyses/Basecall_2D_000/HairpinAlign/Alignment" #01
LOG = "/Analyses/Basecall_2D_000/Log"


BASENAME = "/Analyses/Basecall_2D_000/Configuration/general/basename" #01

CHANNEL = "/Analyses/Basecall_2D_000/Configuration/general/channel" #01
CHANNEL_NUMBER = "/UniqueGlobalKey/channel_id/channel_number" #01

FILE_N = "/Analyses/Basecall_2D_000/Configuration/general/file_number" #01 
READ_ID = "/Analyses/Basecall_2D_000/Configuration/general/read_id" #01

TAG = "/Analyses/Basecall_2D_000/Configuration/general/tag" #01



MODEL_TYPE = "/Analyses/Basecall_2D_000/Configuration/general/model_type" #01


COMP_CALLED = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/called_events" #01
COMP_Q = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/mean_qscore" #01
COMP_N_EVENTS = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/num_events"  #01
COMP_N_SKIPS = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/num_skips" #01
COMP_N_STAYS = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/num_stays" #01
COMP_SEQ_LEN = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement/sequence_length" #01

TEMP_CALLED = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/called_events" #01
TEMP_Q = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/mean_qscore" #01
TEMP_N_EVENTS = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/num_events" #01
TEMP_N_SKIPS = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/num_skips" #01
TEMP_N_STAYS = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/num_stays" #01
TEMP_SEQ_LEN = "/Analyses/Basecall_2D_000/Summary/basecall_1d_template/sequence_length" #01

TWOD_Q = "/Analyses/Basecall_2D_000/Summary/basecall_2d/mean_qscore" #01
TWOD_SEQ_LEN = "/Analyses/Basecall_2D_000/Summary/basecall_2d/sequence_length" #01

EVENT_DETECTION_VERSION = "/Analyses/EventDetection_000/version" #01
TRACK_ID_VERSION = "/UniqueGlobalKey/tracking_id/version" #01
PROTOCOL_VERSION = "/UniqueGlobalKey/tracking_id/protocols_version_name" #01
VERSION_NAME = "/UniqueGlobalKey/tracking_id/version_name" #01



DEVICE_ID = "/UniqueGlobalKey/tracking_id/device_id" #01
ASIC_ID = "/UniqueGlobalKey/tracking_id/asic_id" #01
ASIC_TEMP = "/UniqueGlobalKey/tracking_id/asic_temp" #01
HEATSINK_TEMP = "/UniqueGlobalKey/tracking_id/heatsink_temp" #01

RUN_ID = "/UniqueGlobalKey/tracking_id/run_id" #01

POST_TEMP = '/Analyses/Basecall_2D_000/Summary/post_process_template' #01
POST_COMP = '/Analyses/Basecall_2D_000/Summary/post_process_complement' #01
DURATION = "/Analyses/Basecall_2D_000/BaseCalled_complement/Events/duration" #01
START = "/Analyses/Basecall_2D_000/BaseCalled_complement/Events/start_time" #01
MIN_EVENTS = "/Analyses/Basecall_2D_000/Configuration/general/min_events"  #01 #/Analyses/Basecall_2D_000/Configuration/split_hairpin/min_events 
MAX_EVENTS = "/Analyses/Basecall_2D_000/Configuration/general/max_events"  #01 #/Analyses/Basecall_2D_000/Configuration/split_hairpin/max_events 
SAMPLE_RATE = "/Analyses/Basecall_2D_000/Configuration/general/sampling_rate" #01


class Fast5(object):
    def __init__(self, filename):
        self.filename = filename
        self.filebasename = (".").join(filename.split("/")[-1].split(".")[:-1])
        self.abspath = os.path.abspath(filename)
        self.ATTR_2D = None
        self.ATTR_TEMP = None
        self.ATTR_COMP = None
        self.SPLIT_HAIRPIN = None
        self.is_open = self.open() #self.f5 is the f5 file object
        if self.is_open:# and self.is_nonempty():
            self._basecalling_attempted = self.basecalling_attempted()
            self.LOC_TEMP = None
            self.LOC_COMP = None
            self._has_read = {}
            self._has_read["template"] = self.has_template()
            self._has_read["complement"] = self.has_complement()
            self._has_read["2d"] = self.has_2d()
            self._has_read["input"] = True # arbitrarily true for now
            self.molecule = None
            self.log = None
            self.fq_path = None
            self.given_name = {"template":None, "complement":None, "2d":None}
            self.seq = {"template":None, "complement":None, "2d":None}
            self.fq_sep = {"template":None, "complement":None, "2d":None}
            self.quals = {"template":None, "complement":None, "2d":None}
            self.quals_as_int = {"template":None, "complement":None, "2d":None}
            self.info_name = {"template":None, "complement":None, "2d":None}
            self.base_info_name = None
            self.GENERAL_PATH = None
            self._get_general_path()

    def __str__(self):
        return self.filename
    
    def open(self):
        """
        Open an ONT Fast5 file, assuming HDF5 format
        """
        try:
            self.f5 = h5py.File(self.filename, 'r') 
            return True
        except Exception, e:
##            sys.stderr.write("Cannot open file: %s \n" % self.filename)
##            logger.warning("Cannot open file: %s" % self.filename)
            return False
        
    def is_nonempty(self):
        try:
            self.f5["Analyses"]
            #hdf5 may open but have no information in it, which is effectively corrupt/useless
            # self.f5['Analyses'] sees if the top most level is filled out in fast5, if not
            return True
        except Exception as e:
            return False

    def close(self):
        """
        Close an open an ONT Fast5 file, assuming HDF5 format
        """
        if self.is_open:
            self.f5.close()

    def is_not_corrupt(self):
        return self.is_open
		    
    def f5close(self):
        """
        Close an open an ONT Fast5 file, assuming HDF5 format
        """
        if self.is_open:
            self.hdf5file.close()
    
    def basecalling_attempted(self):
        try:
            self.f5[BASECALL_TEST]
            return True
        except KeyError:
            return False


    def _find_location_of_template(self):
        try:
            self.LOC_TEMP = LOC_TEMP
            self.f5[self.LOC_TEMP]
            return True
        except:
            pass
        try:
            self.LOC_TEMP = LOC_TEMP2
            self.f5[self.LOC_TEMP]
            return True
        except:
            return False

    def has_template(self):
        if self.LOC_TEMP == None:
            return self._find_location_of_template()
        else:
            try:
                self.f5[self.LOC_TEMP]
                return True
            except:
                return False

    def _find_location_of_complement(self):
        try:
            self.LOC_COMP = LOC_COMP
            self.f5[self.LOC_COMP]
            return True
        except:
            pass
        try:
            self.LOC_COMP = LOC_COMP2
            self.f5[self.LOC_COMP]
            return True
        except:
            return False

    def _get_location_path(self, readtype, dataset=""):
        '''readtpye in 2d, template, complement'''
        '''dataset in Fastq, Model, Events, Alignment'''
        if readtype == "2d":
            return LOC_2D + dataset
        elif readtype == "template":
            return self.LOC_TEMP + dataset
        elif readtype == "complement":
            return self.LOC_COMP + dataset
        elif readtype == "input":
##            return INPUT_EVENTS_PATH + "Read_" + self.get_read_id() +"/"+ dataset
            return INPUT_EVENTS_PATH + self.get_read_number() +"/"+ dataset

    def _get_general_path(self):
        try:
            self.GENERAL_PATH = GENERAL
            self.f5[self.GENERAL_PATH]
            return True
        except:
            pass
        try:
            self.GENERAL_PATH = GENERAL001
            self.f5[self.GENERAL_PATH]
            return True
        except:
            return False

    def has_complement(self):
        if self.LOC_COMP == None:
            return self._find_location_of_complement()
        else:
            try:
                self.f5[self.LOC_COMP]
                return True
            except:
                return False

    def has_2d(self):
        try:
            self.LOC_2D = LOC_2D
            self.f5[self.LOC_2D]
            return True
        except KeyError:
            return False
        
    def has_read(self, readtype):
        return self._has_read[readtype]

    def has_reads(self):
        ## the fast5 file is empty of sequence information if it does not even have template read information
        return self.has_read("template") 

    def _define_molecule(self):
        ## if 2d, use 2d as molecule
        ## if both t and c present, use the longer one
        ## if only t present, use t
        ### NOTE: as is - when fast5 is "empty" (lacking even template), this gives answer as template
        ### i.e. assumes non-empty
        if self.has_read("2d"):
            self.molecule = "2d"
        elif self.has_read("complement"):
            if int(self.get_seq_len("complement")) > int(self.get_seq_len("template")):
                self.molecule = "complement"
            else:
                self.molecule = "template"
        else:
            self.molecule = "template"

    def _build_log(self):
        self.log = {}
        self.log["Calibration Strand Log \n"] = ""
        self.log["Hairpin Split Log \n"] = ""
        self.log["Basecall 1D Log \n"] = ""
        self.log["Basecall 2D Log \n"] = ""
        try:
            self.log["Calibration Strand Log \n"] += self.f5["/Analyses/Calibration_Strand_000/Log"][()] + "\n\n"
        except:
            self.log["Calibration Strand Log \n"] += "No Calibration Log found \n\n"
        try:
            self.log["Hairpin Split Log \n"] = self.f5["/Analyses/Hairpin_Split_000/Log"][()] + "\n\n"
        except:
            self.log["Hairpin Split Log \n"] = "No Hairpin Split Log found \n\n"
        try:
            self.log["Basecall 1D Log \n"] = self.f5["/Analyses/Basecall_1D_000/Log"][()] + "\n\n"
        except:
            self.log["Basecall 1D Log \n"] = "No Basecall 1D Log found \n\n"
        try:
            self.log["Basecall 2D Log \n"] = self.f5["/Analyses/Basecall_2D_000/Log"][()] + "\n\n"
        except:
            self.log["Basecall 2D Log \n"] = "No Basecall 2D Log found \n\n"
            

    def has_time_error(self):
        ## repeated blocks of events
        read = [e for e in self.f5[INPUT_EVENTS_PATH]][0]
        path = INPUT_EVENTS_PATH + read + "/Events"
        if self._basecalling_attempted:
            i=2
        else:
            i=0
        t1=self.f5[path][0][i]
        for event in self.f5[path][1:]:
            t2 = event[i]
            if t2 < t1:
                return True
            t1 = t2
        return False


    def use_molecule(self):
        if self.molecule == None:
            self._define_molecule()
        return self.molecule

        
    def get_log_string(self):
        if self.log == None:
            self._build_log()
        logstring = ""
        for log in self.log.keys():
            logstring += log + self.log[log]
        return logstring
    
    def get_log(self):
        if self.log == None:
            self._build_log()
        return self.log
    

    def _get_attr(self, path, attr):
        return self.f5[path].attrs[attr]

    def _get_attr_path(self, readtype):
        '''readtpye in 2d, template, complement'''
        if self._basecalling_attempted:
            if self.ATTR_2D == None:
                self._find_attr_path()
            if readtype == "2d":
                return self.ATTR_2D
            elif readtype == "template":
                return self.ATTR_TEMP
            elif readtype == "complement":
                return self.ATTR_COMP
        else: #not basecalled
            return 

    def _find_attr_path(self):
        self.ATTR_2D = ATTR_2D
        try:
            self.f5[ATTR_TEMP]
            self.ATTR_TEMP = ATTR_TEMP
            self.ATTR_COMP = ATTR_COMP
        except:
            self.ATTR_TEMP = ATTR_TEMP2
            self.ATTR_COMP = ATTR_COMP2
        
    def _find_split_hairpin_path(self):
        try:
            self.f5[SPLIT_HAIRPIN]
            self.SPLIT_HAIRPIN = SPLIT_HAIRPIN
        except:
            self.SPLIT_HAIRPIN = SPLIT_HAIRPIN2

    def _get_split_hairpin_path(self):
        if self.SPLIT_HAIRPIN == None:
            self._find_split_hairpin_path()
        return self.SPLIT_HAIRPIN

    def get_mean_qscore(self, readtype):
        '''readtpye in 2d, template, complement'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = "mean_qscore")


    def get_seq_len(self, readtype):
        '''readtpye in 2d, template, complement'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = "sequence_length")


    def get_num_events(self, readtype):
        '''Assumes exists'''
        '''eventstpye in input, template, complement'''
        ''' if readtype is 2d, it just gives num of all events'''
        if readtype == "input":
            return self.f5["/Analyses/EventDetection_000/Reads/" + self.get_read_number() + "/Events"].shape[0]
        elif readtype == "2d":
            return self._get_attr(path = self._get_split_hairpin_path(), attr = "num_events")
        else:
            return self._get_attr(path = self._get_attr_path(readtype), attr = "num_events")


    def get_num_called_events(self, readtype):
        '''Assumes exists'''
        '''readtpye in template, complement (no input)'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = 'called_events')
        

    def get_num_skips(self, readtype):
        '''Assumes exists'''
        '''eventstpye in template, complement (no input)'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = 'num_skips')


    def get_num_stays(self, readtype):
        '''Assumes exists'''
        '''readtpye in template, complement (no input)'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = 'num_stays')

    def get_time_stamp(self):
        try:
            return self.f5['/Analyses/Basecall_2D_000'].attrs['time_stamp']
        except:
            try:
                return self.f5['/Analyses/Basecall_1D_000'].attrs['time_stamp']
            except:
                return False

    def get_channel_number(self):
        return self.f5['/UniqueGlobalKey/channel_id'].attrs['channel_number']
    #alt: return self.f5["/Analyses/Basecall_2D_000/Configuration/general"].attrs["channel"]

    def get_heatsink_temp(self):
        return self.f5[TRACKING_ID].attrs['heatsink_temp'] 


    def get_device_id(self):
        return self.f5[TRACKING_ID].attrs["device_id"]

    def get_asic_id(self):
        return self.f5[TRACKING_ID].attrs["asic_id"]

    def get_asic_temp(self):
        return self.f5[TRACKING_ID].attrs["asic_temp"]

    def get_run_id(self):
        return self.f5[TRACKING_ID].attrs["run_id"]

    def get_tracking_id_version_name(self):
        # usually the same as get_event_detection_version() and get_tracking_id_version()
        return self.f5[TRACKING_ID].attrs["version_name"]

    def get_event_detection_version(self):
        # usually the same as get_tracking_id_version_name() and get_tracking_id_version()
        return self.f5['/Analyses/EventDetection_000'].attrs["version"]

    def get_tracking_id_version(self):
        # usually the same as get_tracking_id_version_name() and get_event_detection_version()
        return self.f5[TRACKING_ID].attrs["version"]

    def get_protocol_version(self):
        # usually similar to others with out trailing barcode
        return self.f5[TRACKING_ID].attrs["protocols_version_name"]

    def get_basecall_version(self):
        try:
            return self.f5['/Analyses/Basecall_2D_000'].attrs['version']
        except:
            return "Basecaller Version not found."

    def get_file_version(self):
        return self.f5['/'].attrs['file_version']

    def get_file_number(self):
        return self.f5[self.GENERAL_PATH].attrs["file_number"]

    def get_read_id(self):
        try:
            return self.f5[self.GENERAL_PATH].attrs["read_id"]
        except:
            return self.f5["/Analyses/EventDetection_000/Reads/"+ self.get_read_number()].attrs["read_id"]

    def get_read_number(self): ## June 22, 2016 -- not tested for all yet -- developing for non-basecalled files
        return self.f5["/Analyses/EventDetection_000/Reads/"].keys()[0]

    def get_tag(self):
        return self.f5[self.GENERAL_PATH].attrs["tag"]

    def get_model_type(self):
        return self.f5[self.GENERAL_PATH].attrs['model_type']


    def get_basename(self):
        return self.f5[self.GENERAL_PATH].attrs["basename"]




    def _parse_fastq_info(self, readtype):
        if self.seq[readtype] == None:
            try:
                (self.given_name[readtype], self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]) = self.f5[self._get_location_path(readtype,"Fastq")][()].strip().split('\n')
                self.given_name[readtype] = self.given_name[readtype].lstrip('@')
            except:
                (self.given_name[readtype], self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]) = None, None, None, None


    def _get_pore_info_name(self, readtype):
        ## TODO -- this will store most read-pertinent info
        if self.base_info_name == None:
            self._get_base_info_name()
        if self.info_name[readtype] == None:
            info = []
            info.append(readtype)
            self.info_name[readtype] = ("_").join(info)
        return self.info_name[readtype] + "_" + self.base_info_name
    
    def _get_base_info_name(self):
        info = []
        info.append("channel"+self.get_channel_number())
        info.append("file"+self.get_file_number())
        info.append("run"+self.get_run_id())
        info.append("asic"+self.get_asic_id())
##            info.append(("-").join(self.get_time_stamp().split()))
        #could add start time, duration to help distinguish...
        self.base_info_name = ("_").join(info)
        
    def get_base_info_name(self):
        if self.base_info_name == None:
            self._get_base_info_name()
        return self.base_info_name
    
    def get_fastq(self, readtype):
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            return '\n'.join(['@'+self._get_pore_info_name(readtype), self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]])

    def get_fasta(self, readtype):
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            return '\n'.join(['>'+self._get_pore_info_name(readtype), self.seq[readtype]])

    def get_quals(self, readtype):
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            return '\n'.join(['>'+self._get_pore_info_name(readtype), self.quals[readtype]])

    def get_quals_as_int(self, readtype):
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            if self.quals_as_int[readtype] == None:
                self.quals_as_int[readtype] = StringIO.StringIO()
                SeqIO.convert(StringIO.StringIO('\n'.join(['@'+self._get_pore_info_name(readtype), self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]])), "fastq", self.quals_as_int[readtype], "qual")
            return self.quals_as_int[readtype].getvalue().rstrip()

    def get_start_time(self, which):
        ## which in experiment, input, template, complement
        if which == "experiment":
            return self.f5["/UniqueGlobalKey/tracking_id"].attrs["exp_start_time"]
        elif which == "input":
            try:
                return self.f5["/Analyses/EventDetection_000/Reads/Read_" + self.get_read_id()].attrs["start_time"]
            except:
                return self.f5["/Analyses/EventDetection_000/Reads/" + self.get_read_number()].attrs["start_time"]
        elif self.has_read(which) and (which == "template" or which == "complement"):
            return self.f5[self._get_location_path(which,"Events")].attrs["start_time"]


    def get_events(self, readtype):
        return self.f5[self._get_location_path(readtype,"Events")][()]

    def get_events_string(self, readtype):
        return ("\n").join([("\t").join((str(f) for f in e))  for e in self.get_events(readtype)])

    def get_events_header(self, readtype): ## 20160611-inprogress, make tests
        return self.f5[self._get_location_path(readtype,"Events")].dtype.names

    def get_events_header_string(self, readtype): ## 20160611-inprogress, make tests
        return ("\t").join([str(e) for e in self.get_events_header(readtype)])

    def get_model(self, readtype):
        return self.f5[self._get_location_path(readtype,"Model")][()]

    def get_model_string(self, readtype):
        return ("\n").join([("\t").join((str(f) for f in e))  for e in self.get_model(readtype)])

    def get_model_attrs(self, readtype): ## not in tests
        return self.f5[self._get_location_path(readtype,"Model")].attrs

    def get_model_attrs_string(self, readtype): ## not in tests
        return ("\n").join([ ("\t").join([readtype + "_" +str(attr), str(self.get_model_attrs(readtype)[attr])]) for attr in self.get_model_attrs(readtype) ])

    def get_2d_alignment(self):
        if self.has_read("2d"):
            return self.f5[self._get_location_path("2d","Alignment")][()]

    def get_2d_alignment_string(self):
        if self.has_read("2d"):
            return ("\n").join([("\t").join((str(f) for f in e))  for e in self.get_2d_alignment()])


    def get_strand_score(self, readtype):
        '''Assumes exists'''
        '''readtpye in template, complement (no input)'''
        if self.has_read(readtype):
            return self._get_attr(path = self._get_attr_path(readtype), attr = 'strand_score')
        else:
            return "-"

    def get_num_stutters(self, readtype):
        '''Assumes exists'''
        '''readtpye in template, complement (no input)'''
        try:
            if eventstype == "template" and self.has_read("template"):
                return self.f5[POST_TEMP].attrs['num_stutters_found']
            elif eventstype == "complement" and self.has_read("complement"):
                return self.f5[POST_COMP].attrs['num_stutters_found']
        except: ##numstutters only in newer (or some) fast5s
            return "-"
        
    def get_min_tc_ratio(self):
        try:
            return self.f5["/Analyses/Basecall_2D_000/Configuration/hairpin_align"].attrs['min_ratio']
        except:
            try:
               return self.f5["/Analyses/Basecall_1D_000/Configuration/hairpin_align"].attrs['min_ratio']
            except:
                return False

    def get_max_tc_ratio(self):
        try:
            return self.f5["/Analyses/Basecall_2D_000/Configuration/hairpin_align"].attrs['max_ratio']
        except:
            try:
               return self.f5["/Analyses/Basecall_1D_000/Configuration/hairpin_align"].attrs['max_ratio']
            except:
                return False

    def get_tc_ratio_limits(self):
        return self.get_min_tc_ratio(), self.get_max_tc_ratio()

    def get_basecaller_min_events(self):
        return self.f5[self.GENERAL_PATH].attrs['min_events']
        
    def get_basecaller_max_events(self):
        return self.f5[self.GENERAL_PATH].attrs['max_events']

    def get_basecaller_events_limits(self):
        return self.get_basecaller_min_events(), self.get_basecaller_max_events()


F5_TMP_DIR = ".fast5tools_tmp_dir"
F5_TMP_DIR = "fast5tools_tmp_dir"
class Fast5List(object):
    def __init__(self, fast5list, tar_filenames_only=False):
        #ensure type is list
        if isinstance(fast5list, list):
                self.fast5list = fast5list
        elif isinstance(fast5list, str):
                self.fast5list = [fast5list]
        self._tars_detected = False
        self.tar_filenames_only = tar_filenames_only
        self.nfiles = None
        self.allfiles = None
##        self._expand_list() ##TODO(?) - make expand fofn optional, so not wasting time in most situations
        self._extract_fast5_files()
        

    def __iter__(self):
        return self

    def next(self):
        try:
            return Fast5(self.allfiles.next())
        except Exception as e:
            if self._tars_detected and os.path.exists(F5_TMP_DIR):
                shutil.rmtree(F5_TMP_DIR)
            raise StopIteration

    ## I think ultimately one would need to be able to store info about which tarfile object a tar member belongs to (in case more than 1 tarfile)
    ## then when going over the list one would need to check if its an existing file (e.g. os.exists(file)) and if not, is it in one of the tars...
    ## if so, then extract it, load it into fast5. delete file. return fast5 object.
    ## ultimately would be better to just be able to extract into memory and use directly.... but cant find a way to do that.
    ## can also perhaps keep the tar unpacked.... not pre-expand... have the next function check:
    ##  Is this fast5? is it an fofn? is it 

    def _initialize_tar_tmp_dir(self):
        if os.path.isdir(F5_TMP_DIR):
            shutil.rmtree(F5_TMP_DIR)
        self._tars_detected = True
        os.mkdir(F5_TMP_DIR)

##    def _expand_list(self):
##        #Goal: find fofns and tarballs in filelist and expand
##        # Approach:
##        # make second list (empty to start)
##        # go through first list
##        # if item not fofn, and not a tar file, add to 2nd list
##        # ifit is either fofn or tar -- expand it -- then add all to second list
##        #overwrite first list with second
##        expanded_fast5list = []
##        for e in self.fast5list:
##            if e.endswith(".fofn"):
##                f = open(e,'r')
##                files = [line.strip() for line in f.readlines()]
##                f.close()
##                expanded_fast5list += files
##            elif tarfile.is_tarfile(e):
##                if not self._tars_detected and not self.tar_filenames_only:
##                    self._initialize_tar_tmp_dir()
##                f = tarfile.open(e)
##                if self.tar_filenames_only:
##                    files = [os.path.join(e,fname) for fname in f.getnames() if fname.endswith(".fast5")]
##                else: ## will be using files in tarball
##                    f.extractall(path=F5_TMP_DIR)
##                    ## os.path.basename(filename).endswith('.fast5') and not os.path.basename(filename).startswith('.')
##                    files = [os.path.join(F5_TMP_DIR, fname) for fname in f.getnames() if fname.endswith(".fast5")]
##                f.close()
##                expanded_fast5list += files
##            else:
##                expanded_fast5list.append(e)
##        self.fast5list = expanded_fast5list

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
            if fname.endswith(".fast5"):
                files.append(fname)
            elif os.path.isdir(fname):
                files += self._expand_dir(fname)
            elif tarfile.is_tarfile(fname):
                files += self._expand_tar(fname)
            else: ## line in FOFN is ignored
                pass
        return files

    def _expand_dir(self, d):
        pattern = d + '/' + '*.fast5'
        files = glob(pattern)
        return files

    def _expand_tar(self, tarball):
        if not self._tars_detected and not self.tar_filenames_only:
            self._initialize_tar_tmp_dir()
        f = tarfile.open(tarball)
        if self.tar_filenames_only:
            files = [os.path.join(tarball,fname) for fname in f.getnames() if fname.endswith(".fast5")]
        else: ## will be using files in tarball
            f.extractall(path=F5_TMP_DIR) ## in situations where tarball includes many big non-fast5 files, this may not be best way.
            ## os.path.basename(filename).endswith('.fast5') and not os.path.basename(filename).startswith('.')
            files = [os.path.join(F5_TMP_DIR, fname) for fname in f.getnames() if fname.endswith(".fast5")]
        f.close()
        return files

    def _extract_fast5_files(self):
        ##TODO -- also handle FOFNs in list -- i.e. can take FOFN with or without mixture of others
        ##TODO -- handle compressed files, tarchives, etc
        # can be
        # [f.fast5] -- list of a single file
        # [f5dir/] -- list of a single dir
        # [f1.fast5, ..., fN.fast5] -- list of a multiple files
        # [f5dir1, ..., f5dirN/] -- list of multiple dirs
        # [f1.fast5, f5dir1/, f5dir2/, f2.fast5, ..., fN.fast5, f5dirM/] -- mixed list of N files and M dirs
        self.files = []
        for e in self.fast5list:
            if e.endswith(".fast5"):
                self.files.append(e)
            elif e.endswith(".fofn"):
                self.files += self._expand_fofn(e)
            elif os.path.isdir(e):
##                pattern = e + '/' + '*.fast5'
##                files = glob(pattern)
##                self.files += files
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














