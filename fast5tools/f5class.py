## JOHN URBAN (2015, 2016)

#info
import h5py, os, sys, tarfile, shutil
import cStringIO as StringIO
from Bio import SeqIO
from glob import glob
from random import randint
import numpy as np

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

BASECALL_TEST1 = '/Analyses/Basecall_1D_000' #01
BASECALL_TEST2 = '/Analyses/Basecall_2D_000'

GENERAL = "/Analyses/Basecall_2D_000/Configuration/general"
GENERAL001 = "/Analyses/Basecall_1D_000/Configuration/general"
GENERAL_R94_01 = "/Analyses/Basecall_1D_000/Configuration/basecall_1d"

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


TOMBO="/Analyses/RawGenomeCorrected_000/"
TOMBO_BC=TOMBO + "BaseCalled_template/"
TOMBO_ALN=TOMBO_BC + 'Alignment/'
TOMBO_EVENTS=TOMBO_BC+'Events/'

class Fast5(object):
    def __init__(self, filename):
        self.filename = filename
        self.filebasename = (".").join(filename.split("/")[-1].split(".")[:-1]) ## takes entire basename except ".fast5"
        self.abspath = os.path.abspath(filename)
        self.ATTR_2D = None
        self.ATTR_TEMP = None
        self.ATTR_COMP = None
        self.SPLIT_HAIRPIN = None
        self.is_open = self.open() #self.f5 is the f5 file object
        if self.is_open:# and self.is_nonempty():
            self.file_version = None
            self._basecalling_attempted = self.basecalling_attempted()
            self.LOC_TEMP = None
            self.LOC_COMP = None
            self._has_read = {}
            self._has_read["template"] = self.has_template()
            self._has_read["complement"] = self.has_complement()
            self._has_read["2d"] = self.has_2d()
            self._has_read["input"] = True # arbitrarily true for now
            self.molecule = None
            self.molequal = None
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

##    def get_file_version(self):
##        return self.f5['/'].attrs['file_version']

    def get_file_version(self):
        if self.file_version is None:
            try:
                track_version = self.get_tracking_version()
                self.file_version = float(self.f5.attrs['file_version'])
                if float( track_version.split(".")[0] ) < 1 and self.file_version == float(1):
                    ## this is  very early file, call it version 0
                    self.file_version = 0.0
                
            except:
                #If file_version cannot be found, call it 0.0
                #This should not be the case as I see it in all files from 2014 to present
                self.file_version = 0.0
        return self.file_version

    def get_tracking_version(self):
        # This path exists in all file versions 2014 - Nov2017 so far
        return self.f5['/UniqueGlobalKey/tracking_id/'].attrs['version']
            
    def basecalling_detected(self):
        return self._basecalling_attempted

    def is_not_corrupt(self):
        return self.is_open
		    
    def f5close(self):
        """
        Close an open an ONT Fast5 file, assuming HDF5 format
        """
        if self.is_open:
            self.hdf5file.close()
    
    def basecalling_attempted(self):
        detected = 0
        try:
            self.f5[BASECALL_TEST1]
            detected += 1
        except KeyError:
            pass
        try:
            self.f5[BASECALL_TEST2]
            detected += 1
        except KeyError:
            pass
        if detected >= 1:
            return True
        else:
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
            pass
        try:
            self.GENERAL_PATH = GENERAL_R94_01
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

    def _define_molequal(self): ## added July 6, 2016 -- copied/modified from define_molecule
        ## if 2d, use 2d as molecule
        ## **if both t and c present, use the one of higher quality** <- how it differs from molecule
        ## if only t present, use t
        ### NOTE: as is - when fast5 is "empty" (lacking even template), this gives answer as template
        ### i.e. assumes non-empty
        if self.has_read("2d"):
            self.molequal = "2d"
        elif self.has_read("complement"):
            if int(self.get_mean_qscore("complement")) > int(self.get_mean_qscore("template")):
                self.molequal = "complement"
            else:
                self.molequal = "template"
        else:
            self.molequal = "template"


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

    def use_molequal(self):
        if self.molequal == None:
            self._define_molequal()
        return self.molequal
        
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
            try:
                return self.f5["/Analyses/EventDetection_000/Reads/" + self.get_read_number() + "/Events"].shape[0]
            except: ## possibly later file version when this was deprecated
                return '-'
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
        '''Note that in later file versions, num_skips is not provided'''
        '''Need to compute num skips from num moves > 1'''
        '''...'''
        try:
            return self.num_skips
        except:
            if self.get_file_version() < 0.6:
                self.num_skips = self._get_attr(path = self._get_attr_path(readtype), attr = 'num_skips')
            else: #later versions requires compute of number of moves > 1
                self.num_skips = sum(self.get_event_moves(readtype)>1)
        return self.num_skips

    def get_num_stays(self, readtype):
        '''Assumes exists'''
        '''readtpye in template, complement (no input)'''
        '''Note that in later file versions, num_stays is not provided'''
        '''Need to compute num stays from num moves == 0'''
        try:
            return self.num_stays
        except:
            if self.get_file_version() < 0.6:
                self.num_stays = self._get_attr(path = self._get_attr_path(readtype), attr = 'num_stays')
            else:
                self.num_stays = sum(self.get_event_moves(readtype) == 0)
        return self.num_stays

    def get_num_steps(self, readtype):
        '''Assumes exists'''
        '''readtpye in template, complement (no input)'''
        '''As far as I know - both early and late do not have this.'''
        '''Need to compute num steps from num moves == 1'''
        try:
            return self.num_steps
        except:
            self.num_steps = sum(self.get_event_moves(readtype) == 1)
        return self.num_steps


    def get_proportion_skips(self, readtype):
        ''' Summing up num skips, steps, and stays gives number of called events. Therefore, num_called_events is the denominator.'''
        return self.get_num_skips(readtype)/float(self.get_num_called_events(readtype))

    def get_proportion_steps(self, readtype):
        ''' Summing up num skips, steps, and stays gives number of called events. Therefore, num_called_events is the denominator.'''
        return self.get_num_steps(readtype)/float(self.get_num_called_events(readtype))

    def get_proportion_stays(self, readtype):
        ''' Summing up num skips, steps, and stays gives number of called events. Therefore, num_called_events is the denominator.'''
        return self.get_num_stays(readtype)/float(self.get_num_called_events(readtype))

    def get_skip_prob(self, readtype):
        '''Assumes exists'''
        '''eventstpye in template, complement (no input)'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = 'skip_prob')        

    def get_stay_prob(self, readtype):
        '''Assumes exists'''
        '''eventstpye in template, complement (no input)'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = 'stay_prob')        

    def get_step_prob(self, readtype):
        '''Assumes exists'''
        '''eventstpye in template, complement (no input)'''
        return self._get_attr(path = self._get_attr_path(readtype), attr = 'step_prob')        


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


    def get_file_number(self):
        return self.f5[self.GENERAL_PATH].attrs["file_number"]

    def get_read_id(self):
        try:
            return self.f5[self.GENERAL_PATH].attrs["read_id"]
        except:
            return self.f5["/Analyses/EventDetection_000/Reads/"+ self.get_read_number()].attrs["read_id"]

    def get_read_number(self): ## June 22, 2016 -- not tested for all yet -- developing for non-basecalled files
        try:
            return self.f5["/Analyses/EventDetection_000/Reads/"].keys()[0]
        except:
            pass
        
        try:
            ## Nov 14, 2017 -- deal with R9.4 and R9.5
            return str(self.f5['/Raw/Reads'].keys()[0])
                #.split("_")[-1]
        except:
            pass

    def get_tag(self):
        return self.f5[self.GENERAL_PATH].attrs["tag"]

    def get_model_type(self):
        if self.basecalling_detected():
            if self.GENERAL_PATH.endswith('general'):
                return self.f5[self.GENERAL_PATH].attrs['model_type']
            elif self.GENERAL_PATH.endswith('basecall_1d'): ## basically a versioning problem that I can rewrite now that I track that
                try:
                    return str(self.f5[self.GENERAL_PATH].attrs['template_model']) + "_" + str(self.f5[self.GENERAL_PATH].attrs['complement_model'])
                except:
                    pass
                try:
                    return str(self.f5[self.GENERAL_PATH].attrs['template_model'])
                except:
                    pass
                try:
                    return str(self.f5[self.GENERAL_PATH].attrs['model'])
                except:
                    raise Exception('Have not found model type in '+ self.GENERAL_PATH + ' for:' + self.filebasename)
                    quit() ## Don't allow it to continue
        else:
            return "NA"

    def get_model_path(self):
        if self.basecalling_detected():
            model_path = str(self.f5[self.GENERAL_PATH].attrs['model_path'])
            template =  str(self.f5[self.GENERAL_PATH].attrs['template_model'])
            complement = ''
            if self.has_complement:
                complement =  str(self.f5[self.GENERAL_PATH].attrs['complement_model'])
            return ("\t").join( [model_path, template, complement] )
        else:
            return "NA"    


    def get_basename(self):
        return self.f5[self.GENERAL_PATH].attrs["basename"]




    def _parse_fastq_info(self, readtype):
        if self.seq[readtype] == None:
            try:
                (self.given_name[readtype], self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]) = self.f5[self._get_location_path(readtype,"Fastq")][()].strip().split('\n')
                self.given_name[readtype] = self.given_name[readtype].lstrip('@')
            except:
                (self.given_name[readtype], self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]) = None, None, None, None


    def _get_info_name(self, readtype):
        info = []
        info.append(readtype)
        info.append("len:"+str(self.get_seq_len(readtype)))
        info.append("Q:"+str(self.get_mean_qscore(readtype)))
        self.info_name[readtype] = ("|").join(info)

    
    def get_info_name(self, readtype):
        if self.info_name[readtype] == None:
            self._get_info_name(readtype)
        return self.info_name[readtype]

    def _get_base_info_name(self):
        info = []
        info.append("channel:"+self.get_channel_number())
        info.append((":").join(self.get_read_number().split("_")))
##        info.append("file_"+self.get_file_number())
        info.append("asic:"+self.get_asic_id())
        info.append("run:"+self.get_run_id()) ## Though ASIC ID is different for different, if someone re-uses a flowcell, then it will not be. But I believe the RunID will be.... unless it is given the same tag dependent on ASIC ID and/or other constants
        info.append("device:"+self.get_device_id())
        info.append("model:"+self.get_model_type())
##            info.append(("-").join(self.get_time_stamp().split()))
        #could add start time, duration to help distinguish...
        self.base_info_name = ("|").join(info)
        
    def get_base_info_name(self):
        if self.base_info_name == None:
            self._get_base_info_name()
        return self.base_info_name       

    def _get_pore_info_name(self, readtype):
        ## TODO -- this will store most read-pertinent info
##        if self.base_info_name == None:
##            self._get_base_info_name()
##        if self.info_name[readtype] == None:
##            info = []
##            info.append(readtype)
##            info.append("len:"+str(self.get_seq_len(readtype)))
##            info.append("Q:"+str(self.get_mean_qscore(readtype)))
##            self.info_name[readtype] = ("|").join(info)
##        return self.info_name[readtype] + "|" + self.base_info_name
        return self.get_info_name(readtype) + "|" + self.get_base_info_name()

    def get_pore_info_name(self, readtype):
        ''' This is redundant with _get_pore_info_name() - wanted to transition to using this one w/o breaking things dependent on the other'''
        return self.get_info_name(readtype) + "|" + self.get_base_info_name()

    ## Nov 17 - I am adding other types of naming schemes, but I dont believe any of these need to be stored as I have been doing - they can just be created on the fly
    def get_read_stats_name(self, readtype):
        ''' Abridged version of pore_info_name: readtype:len:Q:channel:read. Missing asic:run:device:model.'''
        info = []
        info.append("channel:"+self.get_channel_number())
        info.append((":").join(self.get_read_number().split("_")))
        return self.get_info_name(readtype) + "|" + ("|").join(info)

    def get_event_stats_name(self, readtype):
        info = []
        info.append( readtype+"_events:"+str(self.get_num_events(readtype)))
        info.append( readtype+"_calledevents:"+str(self.get_num_called_events(readtype)))
        info.append( readtype+"_skips:"+str(self.get_num_skips(readtype)))
        info.append( readtype+"_stays:"+str(self.get_num_stays(readtype)))
        info.append( readtype+"_steps:"+str(self.get_num_steps(readtype)))
        info.append( readtype+"_prop_skips:"+str(self.get_proportion_skips(readtype)))
        info.append( readtype+"_prop_stays:"+str(self.get_proportion_stays(readtype)))
        info.append( readtype+"_prop_steps:"+str(self.get_proportion_steps(readtype)))
        info.append( readtype+"_skip_prob:"+str(self.get_skip_prob(readtype)))
        info.append( readtype+"_stay_prob:"+str(self.get_stay_prob(readtype)))
        info.append( readtype+"_step_prob:"+str(self.get_step_prob(readtype)))
        info.append( readtype+"_start_time:"+str(self.get_start_time(readtype)))
        info.append( readtype+"time_length:"+str( sum(self.get_event_lengths(readtype)) ) )
        info.append( 'raw_duration:'+ str(self.get_raw_duration() ) )
        return self.get_info_name(readtype) + "|" + ("|").join(info)


    def get_read_and_event_stats_name(self, readtype):
        ''' .'''
        info = []
        info.append( readtype+"_events:"+str(self.get_num_events(readtype)))
        info.append( readtype+"_calledevents:"+str(self.get_num_called_events(readtype)))
        info.append( readtype+"_skips:"+str(self.get_num_skips(readtype)))
        info.append( readtype+"_stays:"+str(self.get_num_stays(readtype)))
        info.append( readtype+"_steps:"+str(self.get_num_steps(readtype)))
        info.append( readtype+"_prop_skips:"+str(self.get_proportion_skips(readtype)))
        info.append( readtype+"_prop_stays:"+str(self.get_proportion_stays(readtype)))
        info.append( readtype+"_prop_steps:"+str(self.get_proportion_steps(readtype)))
        info.append( readtype+"_skip_prob:"+str(self.get_skip_prob(readtype)))
        info.append( readtype+"_stay_prob:"+str(self.get_stay_prob(readtype)))
        info.append( readtype+"_step_prob:"+str(self.get_step_prob(readtype)))
        info.append( readtype+"_start_time:"+str(self.get_start_time(readtype)))
        info.append( readtype+"time_length:"+str( sum(self.get_event_lengths(readtype)) ) )
        info.append( 'raw_duration:'+ str(self.get_raw_duration() ) )
        return self.get_read_stats_name(readtype) + "|" + ("|").join(info)
    
    def get_pore_info_name_with_abspath(self, readtype):
        return self.get_pore_info_name(readtype)+"|filename:"+self.abspath

    def get_pore_info_name_with_filebasename(self, readtype):
        return self.get_pore_info_name(readtype)+"|filename:"+self.filebasename

    def get_read_stats_name_with_abspath(self, readtype):
        return self.get_read_stats_name(readtype)+"|filename:"+self.abspath

    def get_read_stats_name_with_filebasename(self, readtype):
        return self.get_read_stats_name(readtype)+"|filename:"+self.filebasename

    def _get_fastx_name_and_comments(self, readtype, name=False,comments=False):
        if not name:
            name = self.get_pore_info_name(readtype)
        if not comments:
            comments = ''
        else:
            comments = "\t" + comments
        return name, comments
    
    def get_fastq(self, readtype, name=False, comments=False):
        #Nov 17 - transitioning to having this use any name given, pore_info by deault
        # Thus all extra fastq fxns below can be replicated by providing the appropriate name here
        # They will be kept as conveniences.
        # This new approach will allow much much mor eflexibility in the future.
        name, comments = self._get_fastx_name_and_comments(readtype, name, comments)
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            return '\n'.join(['@'+name+comments, self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]])

    def get_fastq_with_abspath(self, readtype, comments=False):
        name = self.get_pore_info_name_with_abspath(readtype)
        return self.get_fastq(readtype, name, comments=comments)

    def get_fastq_only_abspath(self, readtype, comments=False):
        return self.get_fastq(readtype, name=self.abspath, comments=comments)
    
    def get_fastq_with_filename(self, readtype, comments=False):
        name = self.get_pore_info_name_with_filebasename(readtype)
        return self.get_fastq(readtype, name, comments=comments)
    
    def get_fastq_only_filename(self, readtype, comments=False):
        return self.get_fastq(readtype, name=self.filebasename, comments=comments)
    
    def get_fasta(self, readtype, name=False, comments=False):
        #Nov 17 - transitioning to having this use any name given, pore_info by deault
        name, comments = self._get_fastx_name_and_comments(readtype, name, comments)
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            return '\n'.join(['>'+name+comments, self.seq[readtype]])

    def get_fasta_with_abspath(self, readtype, comments=False):
        name = self.get_pore_info_name_with_abspath(readtype)
        return self.get_fasta(readtype, name, comments=comments)

    def get_fasta_only_abspath(self, readtype, comments=False):
        return self.get_fasta(readtype, name=self.abspath, comments=comments)

    def get_fasta_with_filename(self, readtype, comments=False):
        name = self.get_pore_info_name_with_filebasename(readtype)
        return self.get_fasta(readtype, name, comments=comments)


    def get_fasta_only_filename(self, readtype, comments=False):
        return self.get_fasta(readtype, name=self.filebasename, comments=comments)


    def get_quals(self, readtype, name=False, comments=False):
        #Nov 17 - transitioning to having this use any name given, pore_info by deault
        name, comments = self._get_fastx_name_and_comments(readtype, name, comments)
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            return '\n'.join(['>'+name+comments, self.quals[readtype]])

    def get_quals_with_abspath(self, readtype, comments=False):
        name = self.get_pore_info_name_with_abspath(readtype)
        return self.get_quals(readtype, name, comments=comments)

    def get_quals_only_abspath(self, readtype, comments=False):
        return self.get_quals(readtype, name=self.abspath, comments=comments)


    def get_quals_with_filename(self, readtype, comments=False):
        name = self.get_pore_info_name_with_filebasename(readtype)
        return self.get_quals(readtype, name, comments=comments)

    def get_quals_only_filename(self, readtype, comments=False):
        return self.get_quals(readtype, name=self.filebasename, comments=comments)

    ## NOTE: did not add the abs path (or filename) options for thing below -- Add if needed.
    def get_quals_as_int(self, readtype):
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            if self.quals_as_int[readtype] == None:
                self.quals_as_int[readtype] = StringIO.StringIO()
                SeqIO.convert(StringIO.StringIO('\n'.join(['@'+self._get_pore_info_name(readtype), self.seq[readtype], self.fq_sep[readtype], self.quals[readtype]])), "fastq", self.quals_as_int[readtype], "qual")
            return self.quals_as_int[readtype].getvalue().rstrip()

    def falconize_name(self, readtype, zmw_num, style="old"):
        if style == "old":
            moviename = "m000_000"
            otherinfo = self._get_pore_info_name(readtype)
        elif style == "new":
            info = []
            info.append("asic:"+self.get_asic_id())
            info.append("run:"+self.get_run_id())
            info.append("device:"+self.get_device_id())
            info.append("model:"+self.get_model_type())
            moviename = ("|").join(info)
            info = []
            info.append(readtype)
            info.append("Q:"+str(self.get_mean_qscore(readtype)))
            info.append((":").join(self.get_read_number().split("_")))
            info.append("channel:"+self.get_channel_number())
            otherinfo = ("|").join(info)
        return moviename + "/" + str(zmw_num) + "/0_"+str(self.get_seq_len(readtype)) + " " + otherinfo
    

    def get_falcon_fasta(self, readtype, zmw_num=None, style="old"):
        if self.has_read(readtype):
            self._parse_fastq_info(readtype)
            if zmw_num == None:
                zmw_num = randint(0,1000000000000)
            return '\n'.join([">"+self.falconize_name(readtype, zmw_num, style), self.seq[readtype]])


    def get_start_time(self, which):
        ## which in experiment, input, template, complement
        if which == "experiment":
            return self.get_experiment_start_time()
        elif which == "input":
            return self.get_input_events_start_time()
        elif self.has_read(which) and (which == "template" or which == "complement"):
            return self.get_read_events_start_time(readtype=which,index=0)

    def get_experiment_start_time(self):
        return self.f5["/UniqueGlobalKey/tracking_id"].attrs["exp_start_time"]

    def get_input_events_start_time(self):
        ## Nov 2017: 'input' breaks this in latest files
        if self.get_file_version() < 0.6:
            try:
                return self.f5["/Analyses/EventDetection_000/Reads/Read_" + self.get_read_id()].attrs["start_time"]
            except:
                return self.f5["/Analyses/EventDetection_000/Reads/" + self.get_read_number()].attrs["start_time"]
        else:
            ## TODO: better solution here...
            print "Latest files have raw signal. Segmentation events only found after base-calling, and only for read types."
            quit()
##            return self.f5[self._get_location_path("input","Events")]["start"][0]

    def get_read_events_start_time(self, readtype, index=0):
        ''' By default this gives the START TIME -- i.e. time of first event'''
        ''' To get end time, use -1.'''
        ''' Any index between 0 and the vector length is accepted -- but not checked.'''
        if index == 0 and self.get_file_version() < 0.6:
            return self.f5[self._get_location_path(which,"Events")].attrs["start_time"]
        else:
            return self.get_event_start_times(readtype)[index]


    def get_events(self, readtype):
        ## Nov 2017: 'input' breaks this in latest files
        return self.f5[self._get_location_path(readtype,"Events")][()]

    def yield_events(self, readtype): ## 20160705 -- make tests
        yield self.f5[self._get_location_path(readtype,"Events")][()]

    def get_events_string(self, readtype):
        return ("\n").join([("\t").join((str(f) for f in e))  for e in self.get_events(readtype)])

    def get_events_header(self, readtype): ## 20160611-inprogress, make tests
        ## Nov 2017: 'input' breaks this in latest files
        return self.f5[self._get_location_path(readtype,"Events")].dtype.names

    def get_events_header_string(self, readtype): ## 20160611-inprogress, make tests
        return ("\t").join([str(e) for e in self.get_events_header(readtype)])

    def get_event_means(self, readtype):
        return self.f5[self._get_location_path(readtype,"Events")]["mean"]

    def get_event_stdevs(self, readtype):
        try: ##basecalled files report stdv
            return self.f5[self._get_location_path(readtype,"Events")]["stdv"]
        except:
            pass
        try: ## files that are not basecalled report variance
            return self.f5[self._get_location_path(readtype,"Events")]["variance"]**0.5
        except:
            return None
        
    def get_event_variance(self, readtype):
        try:  ## files that are not basecalled report 
            return self.f5[self._get_location_path(readtype,"Events")]["variance"]
        except:
            pass
        try:  ##basecalled files report stdv
            return self.f5[self._get_location_path(readtype,"Events")]["stdv"]**2
        except:
            return None

    def get_event_lengths(self, readtype):
        return self.f5[self._get_location_path(readtype,"Events")]["length"]

    def get_event_start_times(self, readtype):
        return self.f5[self._get_location_path(readtype,"Events")]["start"]

    def get_event_raw_starts(self, readtype):
        ''' Same as get_event_start_times. However, in newer files where raw data is provided,
            "start" is literally the data point it starts on...
            In HDFView I can see that datapoints start at 0 and end at length-1.
            i.e. python indexing'''
        return self.f5[self._get_location_path(readtype,"Events")]["start"]

    def get_event_flags(self): ## only present in non-basecalled files and seems to have shown up recently (e.g. first half of 2016)
        try:
            return self.f5[self._get_location_path(readtype="input",dataset="Events")]["flags"]
        except:
            return None
        
    def get_event_moves(self, readtype):
        return self.f5[self._get_location_path(readtype,"Events")]["move"]

    def get_events_dict(self, readtype):
        ## make it (store in external variables)
        events_dict = {}
        for key in self.get_events_header(readtype):
            events_dict[key] = self.f5[self._get_location_path(readtype,"Events")][key]
        return events_dict
        

    def get_model(self, readtype):
        if self.get_file_version() < 0.6:
            return self.f5[self._get_location_path(readtype,"Model")][()]
        else:
            msg = "Models_now_provided_with_Albacore.\tModels_used_locations:\t" + self.get_model_path() 
            return [[e] for e in msg.split()]

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


    def get_raw_path(self):
        readnum = str(self.get_read_number())
        RAW = '/Raw/Reads/' + readnum
        return RAW

    def get_raw_signal_path(self):
        RAW = self.get_raw_path() + '/Signal/'
        return RAW

    def get_raw_signal(self, median_normalized=False, scale=100.0):
        if median_normalized:
            raw = self.f5[self.get_raw_signal_path()][()]
            median = np.median(raw)
            return scale*raw/median
        else:
            return self.f5[self.get_raw_signal_path()][()]
    
    def get_raw_signal_string(self, delimiter='\n', median_normalized=False):
        return (delimiter).join([str(e) for e in self.get_raw_signal(median_normalized)])


    def get_raw_attribute(self, attr):
        return self.f5[self.get_raw_path()].attrs[attr]
    
    def get_raw_duration(self):
        return self.get_raw_attribute('duration')

    def get_raw_median_before(self):
        return self.get_raw_attribute('median_before')

    def get_raw_read_id(self):
        return self.get_raw_attribute('read_id')

    def get_raw_read_number(self):
        ## somewhat redundant with get_read_number -- but this one is not useful on older files
        ## also output is slightly different "Read_Number" vs "Number"
        return self.get_raw_attribute('read_number')

    def get_raw_start_mux(self):
        return self.get_raw_attribute('start_mux')

    def get_raw_start_time(self):
        return self.get_raw_attribute('start_time')

    def get_segmented_raw_signal(self, readtype='template', median_normalized=False):
        ''' Segments signal according to events.'''
        raw = self.get_raw_signal(median_normalized)
        events = self.get_events_dict(readtype) ## 12/12/17 For now with 1D reads, this is all that makes sense
        rawseg = {}
        rawseg['5primeclip'] = raw[:events['start'][0]]
        sumlen = events['start'][0]
        for i in range(len(events['start'])):
            start = events['start'][i]
            end = events['start'][i] + events['length'][i]
            sumlen += end-start
            rawseg[i] = raw[start:end]
        try:
            rawseg['3primeclip'] = raw[end:]
        except:
            rawseg['3primeclip'] = np.array([])
        return rawseg

    def get_segmented_raw_signal_string(self, includeclips=True, datadelim=' ', eventdelim='\n', readtype='template', median_normalized=False):
        rawseg = self.get_segmented_raw_signal(readtype, median_normalized)
        rawstring = ''
        #add 5'?
        if includeclips and rawseg['5primeclip'].any():
            rawstring += (datadelim).join( [str(e) for e in rawseg['5primeclip']] ) + eventdelim
        #get event-parsed raw
        for i in range(len(rawseg.keys())-2-1): #-2 b/c 5primeclip and 3primeclip; -1 to avoid eventdelim after last one
            rawstring += (datadelim).join( [str(e) for e in rawseg[i]] ) + eventdelim
        #add last event w/o eventdelim
        i+=1
        rawstring += (datadelim).join( [str(e) for e in rawseg[i]] )
        #add 3' ?
        if includeclips and rawseg['3primeclip'].any():
            rawstring +=  eventdelim
            rawstring += (datadelim).join( [str(e) for e in rawseg['3primeclip']] )
        return rawstring

    def get_event_from_raw_signal(self, start, raw, name=None):
        ''' start is calculated elsewhere
            raw is slice of raw signal -- numpy array
            name is optional'''
        mean = raw.mean()
        stdv = raw.std()
        length = len(raw)
        if name is not None:
            return (name, mean, stdv, start, length)
        else:
            return (mean, stdv, start, length)

    def get_segmented_raw_signal_stats(self, includeclips=True, readtype='template', median_normalized=False):
##        print ('\n').join([str(e) for e in 100.0*self.get_raw_signal()/np.median(self.get_raw_signal())])
        rawseg = self.get_segmented_raw_signal(readtype, median_normalized)
        rawstats = []
        sumlen = 0
        #add 5'?
        if includeclips and rawseg['5primeclip'].any():
            event = self.get_event_from_raw_signal(start=0, raw=rawseg['5primeclip'], name='5primeclip')
            rawstats.append( event )
            sumlen += event[4]
        #get event-parsed raw
        for i in range(len(rawseg.keys())-2): #-2 b/c 5primeclip and 3primeclip
            event = self.get_event_from_raw_signal(start=sumlen, raw=rawseg[i], name=i)
            rawstats.append( event )
            sumlen += event[4]
        #add 3' ?
        if includeclips and rawseg['3primeclip'].any():
            event = self.get_event_from_raw_signal(start=sumlen, raw=rawseg['3primeclip'], name='3primeclip')
            rawstats.append( event )
            sumlen += event[4]
##        assert sumlen == self.get_raw_duration()
        assert sumlen == self.get_raw_duration()
        return rawstats
                
    def get_segmented_raw_signal_stats_string(self, includeclips=True, readtype='template', median_normalized=False):
        eventstr = ''
        for event in self.get_segmented_raw_signal_stats(includeclips, readtype, median_normalized):
           eventstr += ('\t').join([str(e) for e in event]) + '\n'
        return eventstr.rstrip()

    def tombo_exists(self):
        try:
            self.f5[TOMBO]
            return True
        except:
            return False

    def tombo_aln_exists(self):
        try:
            self.f5[TOMBO_ALN]
            return True
        except:
            return False

    def tombo_events_exist(self):
        try:
            self.f5[TOMBO_EVENTS]
            return True
        except:
            return False

    def tombo_params_exist(self):
        try:
            self.f5[TOMBO_BC]
            return True
        except:
            return False

    def tombo_successful(self):
        return self.f5[TOMBO_BC].attrs['status'] == 'success'

    def get_tombo_version(self):
        print TOMBO
        return self.f5[TOMBO].attrs['tombo_version']
    
    def get_tombo_alignment_attribute(self, attr):
        return self.f5[TOMBO_ALN].attrs[attr]

    def get_tombo_parameters(self):
        return [(k,v) for k,v in self.f5[TOMBO_BC].attrs.iteritems()]

    def get_tombo_parameter_string(self, delim='\t'):
        return (delim).join([str(e[1]) for e in self.get_tombo_parameters()])

    def get_tombo_parameter_header_string(self, delim='\t'):
        return (delim).join([str(e[0]) for e in self.get_tombo_parameters()])

    def get_tombo_map_position(self):
        '''Must be a fast5 brought through tombo resquiggle/mapping process
        Returns 10-tuple of chr, start, end, strand.....'''
        chrom = self.get_tombo_alignment_attribute('mapped_chrom')
        start = self.get_tombo_alignment_attribute('mapped_start')
        end = self.get_tombo_alignment_attribute('mapped_end')
        strand = self.get_tombo_alignment_attribute('mapped_strand')
        n_match = self.get_tombo_alignment_attribute('num_matches')
        n_mismatch = self.get_tombo_alignment_attribute('num_mismatches')
        n_del = self.get_tombo_alignment_attribute('num_deletions')
        n_ins = self.get_tombo_alignment_attribute('num_insertions')
        clip_start = self.get_tombo_alignment_attribute('clipped_bases_start')
        clip_end = self.get_tombo_alignment_attribute('clipped_bases_end')
        return chrom, start, end, strand, n_match, n_mismatch, n_del, n_ins, clip_start, clip_end

    def get_tombo_map_position_string(self,delim='\t'):
        '''Must be a fast5 brought through tombo resquiggle/mapping process'''
        return (delim).join([str(e) for e in self.get_tombo_map_position()])
        

    def get_tombo_events(self):
        '''Must be a fast5 brought through tombo resquiggle/mapping process'''
        return self.f5[TOMBO_EVENTS][()]


    def get_tombo_events_string(self):
        return ("\n").join([("\t").join((str(f) for f in e))  for e in self.get_tombo_events()])

        
    def get_tombo_events_header(self): 
        return self.f5[TOMBO_EVENTS].dtype.names

    def get_tombo_events_header_string(self): 
        return ("\t").join([str(e) for e in self.get_tombo_events_header()])

    def get_tombo_genomic_events(self):
        '''Must be a fast5 brought through tombo resquiggle/mapping process'''
        chrom = self.get_tombo_alignment_attribute('mapped_chrom')
        start = self.get_tombo_alignment_attribute('mapped_start')
        end = self.get_tombo_alignment_attribute('mapped_end')
        strand = self.get_tombo_alignment_attribute('mapped_strand')
        events = self.get_tombo_events()
        nevents = len(events)
        assert nevents == start-end
        genomic_events = []
        for i in range(nevents):
            pos = start + i ##starts out as i=0, so start+0
            genomic_events.append( tuple(list(event) + [chrom, pos, strand]) )
        return genomic_events

    def get_tombo_genomic_events_string(self):
        return ("\n").join([("\t").join((str(f) for f in e))  for e in self.get_tombo_genomic_events()])
        
    def get_tombo_genomic_events_header_string(self): 
        return ("\t").join([str(e) for e in list(self.get_tombo_events_header()) + ['chr', 'pos', 'strand']])


F5_TMP_DIR = ".fast5tools_tmp_dir"
F5_TMP_DIR = "fast5tools_tmp_dir"
class Fast5List(object):
    def __init__(self, fast5list, tar_filenames_only=False, keep_tar_footprint_small=True):
        #ensure type is list
        if isinstance(fast5list, list):
                self.fast5list = fast5list
        elif isinstance(fast5list, str):
                self.fast5list = [fast5list]
        self._tars_detected = False
        self.tar_filenames_only = tar_filenames_only
        self.keep_tar_footprint_small = keep_tar_footprint_small
        self.nfiles = None
        self.allfiles = None
        self.n_tars = 0
        self.F5_TMP_DIR = None
        self.tars = {} ## for small footprint method
##        self._expand_list() ##TODO(?) - make expand fofn optional, so not wasting time in most situations
        self._extract_fast5_files()

        

    def __iter__(self):
        return self

    def next(self):
        try:
            newfile = self.allfiles.next()
            if self.keep_tar_footprint_small and newfile.startswith("f5tar|"):
                    f5tar, key, tar_member = newfile.split("|")
                    tarkey = "f5tar|" + key + "|"
                    self.tars[tarkey].extract(tar_member, path=self.F5_TMP_DIR)
                    newfile = os.path.join(self.F5_TMP_DIR, tar_member)
                    f5 = Fast5(newfile)
                    os.remove(newfile)
            else:
                f5 = Fast5(newfile)
            return f5
        
        except Exception as e:
            if self._tars_detected and os.path.exists(self.F5_TMP_DIR):
                shutil.rmtree(self.F5_TMP_DIR)
            raise StopIteration
        
           

    ## I think ultimately one would need to be able to store info about which tarfile object a tar member belongs to (in case more than 1 tarfile)
    ## then when going over the list one would need to check if its an existing file (e.g. os.exists(file)) and if not, is it in one of the tars...
    ## if so, then extract it, load it into fast5. delete file. return fast5 object.
    ## ultimately would be better to just be able to extract into memory and use directly.... but cant find a way to do that.
    ## can also perhaps keep the tar unpacked.... not pre-expand... have the next function check:
    ##  Is this fast5? is it an fofn? is it 

    def _initialize_tar_tmp_dir(self):
        if self.F5_TMP_DIR == None:
            seed = str(randint(10000000000,99999999999))
            self.F5_TMP_DIR = F5_TMP_DIR + "_" + seed
        if os.path.isdir(self.F5_TMP_DIR):
            shutil.rmtree(self.F5_TMP_DIR)
        self._tars_detected = True
        os.mkdir(self.F5_TMP_DIR)

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
##                    f.extractall(path=self.F5_TMP_DIR)
##                    ## os.path.basename(filename).endswith('.fast5') and not os.path.basename(filename).startswith('.')
##                    files = [os.path.join(self.F5_TMP_DIR, fname) for fname in f.getnames() if fname.endswith(".fast5")]
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
        self.n_tars += 1
        if self.tar_filenames_only:
            files = [os.path.join(tarball,fname) for fname in f.getnames() if fname.endswith(".fast5")]
            f.close()
        else: ## will be using files in tarball
            if self.keep_tar_footprint_small:
                tarkey = "f5tar|"+str(self.n_tars)+"|"
                self.tars[tarkey] = f
                files = [tarkey+fname for fname in f.getnames() if fname.endswith(".fast5")]
                ## purposely do not close tarfile
            else:                  
                f.extractall(path=self.F5_TMP_DIR) ## in situations where tarball includes many big non-fast5 files, this may not be best way.
                ## os.path.basename(filename).endswith('.fast5') and not os.path.basename(filename).startswith('.')
                files = [os.path.join(self.F5_TMP_DIR, fname) for fname in f.getnames() if fname.endswith(".fast5")]
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














