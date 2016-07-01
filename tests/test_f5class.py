import glob
import os
import unittest

from fast5tools.f5class import Fast5


data_path = "rundata"
fast5filename = "example.fast5"
data_dirs = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11",
             "11b", "12", "13", "14", "15", "16"]


class TestFast5Class(unittest.TestCase):

    def __init__(self, method, filename=""):
        super(TestFast5Class, self).__init__(method)
        self.filename = filename
        self.method = method

    def setUp(self):
        print("\n")
        print("test %s on %s" % (self.method, self.filename))
        self.f5 = Fast5(self.filename)

    def test_has_time_error(self):
        self.f5.has_time_error()

    def test_get_mean_qscore_template(self):
        if self.f5.has_read("template"):
            self.f5.get_mean_qscore("template")


    def test_get_mean_qscore_use_molecule(self):
        if self.f5.has_read("template"):
            self.f5.get_mean_qscore(self.f5.use_molecule())

    def test_get_mean_qscore_complement(self):
        if self.f5._has_read["complement"]:
            self.f5.get_mean_qscore("complement")

    def test_get_mean_qscore_2d(self):
        if self.f5._has_read["2d"]:
            self.f5.get_mean_qscore("2d")

    def test_get_seq_len_template(self):
        if self.f5.has_read("template"):
            self.f5.get_seq_len("template")

    def test_get_seq_len_use_molecule(self):
        if self.f5.has_read("template"):
            self.f5.get_seq_len(self.f5.use_molecule())

    def test_get_seq_len_complement(self):
        if self.f5._has_read["complement"]:
            self.f5.get_seq_len("complement")

    def test_get_seq_len_2d(self):
        if self.f5._has_read["2d"]:
            self.f5.get_seq_len("2d")

    def test_get_num_events_input(self):
        self.f5.get_num_events("input")

    def test_get_num_events_template(self):
        if self.f5.has_read("template"):
            self.f5.get_num_events("template")

    def test_get_num_events_complement(self):
        if self.f5._has_read["complement"]:
            self.f5.get_num_events("complement")
        
    def test_get_num_called_events_template(self):
        if self.f5.has_read("template"):
            self.f5.get_num_called_events("template")

    def test_get_num_called_events_complement(self):
        if self.f5._has_read["complement"]:
            self.f5.get_num_called_events("complement")

    def test_get_num_skips_template(self):
        if self.f5.has_read("template"):
            self.f5.get_num_skips("template")

    def test_get_num_skips_complement(self):
        if self.f5._has_read["complement"]:
            self.f5.get_num_skips("complement")

    def test_get_num_stays_template(self):
        if self.f5.has_read("template"):
            self.f5.get_num_stays("template")

    def test_get_num_stays_complement(self):
        if self.f5._has_read["complement"]:
            self.f5.get_num_stays("complement")

    def test_get_channel_number(self):
        self.f5.get_channel_number()

    def test_get_heatsink_temp(self):
        self.f5.get_heatsink_temp()

    def test_get_model_type(self):
        if self.f5._basecalling_attempted:
            self.f5.get_model_type()

    def test_get_basename(self):
        if self.f5._basecalling_attempted:
            self.f5.get_basename()

    def test_get_filenumber(self):
        if self.f5._basecalling_attempted:
            self.f5.get_file_number()

    def test_get_read_id(self):
        self.f5.get_read_id()

    def test_get_tag(self):
        if self.f5._basecalling_attempted:
            self.f5.get_tag()

    def test_get_fastq_template(self):
        if self.f5.has_read("template"):
            self.f5.get_fastq("template")

    def test_get_fastq_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_fastq("complement")

    def test_get_fastq_2d(self):
        if self.f5.has_read("2d"):
            self.f5.get_fastq("2d")

    def test_get_fasta_complement(self):
        if self.f5.has_read("template"):
            self.f5.get_fasta("template")

    def test_get_fasta_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_fasta("complement")

    def test_get_fasta_2d(self):
        if self.f5.has_read("2d"):
            self.f5.get_fasta("2d")

    def test_get_quals_template(self):
        if self.f5.has_read("template"):
            self.f5.get_quals("template")

    def test_get_quals_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_quals("complement")

    def test_get_quals_2d(self):
        if self.f5.has_read("2d"):
            self.f5.get_quals("2d")

    def test_get_quals_as_int_complement(self):
        if self.f5.has_read("template"):
            self.f5.get_quals_as_int("template")

    def test_get_quals_as_int_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_quals_as_int("complement")

    def test_get_quals_as_int_2d(self):
        if self.f5.has_read("2d"):
            self.f5.get_quals_as_int("2d")

    def test_get_minknow_version(self):
        self.f5.get_tracking_id_version_name()

    def test_get_event_detection_version(self):
        self.f5.get_event_detection_version()

    def test_get_tracking_id_version(self):
        self.f5.get_tracking_id_version()

    def test_get_protocol_version(self):
        self.f5.get_protocol_version()

    def test_get_basecall_version(self):
        self.f5.get_basecall_version()

    def test_get_file_version(self):
        self.f5.get_file_version()

    def test_get_device_id(self):
        self.f5.get_device_id()

    def test_get_asic_id(self):
        self.f5.get_asic_id()

    def test_get_asic_temp(self):
        self.f5.get_asic_temp()

    def test_get_run_id(self):
        self.f5.get_run_id()

    def test_get_time_stamp(self):
        self.f5.get_time_stamp()

    def test_get_start_time_exp(self):
        self.f5.get_start_time("experiment")

    def test_get_start_time_input(self):
        self.f5.get_start_time("input")

    def test_get_start_time_template(self):
        self.f5.get_start_time("template")

    def test_get_start_time_complement(self):
        self.f5.get_start_time("complement")

    def test_get_log(self):
        self.f5.get_log()

    def test_get_log_string(self):
        self.f5.get_log_string()

    def test_get_model_template(self):
        if self.f5.has_read("template"):
            self.f5.get_model("template")

    def test_get_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_model("complement")

    def test_get_model_string_template(self):
        if self.f5.has_read("template"):
            self.f5.get_model_string("template")

    def test_get_model_string_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_model_string("complement")

    def test_get_events_input(self):
        self.f5.get_events("input")

    def test_get_events_template(self):
        if self.f5.has_read("template"):
            self.f5.get_events("template")

    def test_get_events_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_events("complement")

    def test_get_events_string_input(self):
        self.f5.get_events_string("input")

    def test_get_events_string_template(self):
        if self.f5.has_read("template"):
            self.f5.get_events_string("template")

    def test_get_events_string_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_events_string("complement")

    def test_get_events_header_input(self):
        self.f5.get_events_header("input")

    def test_get_events_header_template(self):
        if self.f5.has_read("template"):
            self.f5.get_events_header("template")

    def test_get_events_header_complement(self):
        if self.f5.has_read("complement"):
            self.f5.get_events_header("complement")

    def test_get_events_header_string_input(self):
        self.f5.get_events_header_string("input")

    def test_get_events_header_string_template(self):
        if self.f5.has_read("template"):
            self.f5.get_events_header_string("template")

    def test_get_events_header_string(self):
        if self.f5.has_read("complement"):
            self.f5.get_events_header_string("complement")

    def test_get_2d_alignment(self):
        self.f5.get_2d_alignment()

    def test_get_2d_alignment_string(self):
        self.f5.get_2d_alignment_string()

    def test_get_strand_score_template(self):
        self.f5.get_strand_score("template")

    def test_get_strand_score_complement(self):
        self.f5.get_strand_score("complement")

    def test_get_num_stutters_template(self):
        self.f5.get_num_stutters("template")

    def test_get_num_stutters_complement(self):
        self.f5.get_num_stutters("complement")

    def test_get_min_tc_ratio(self):
        self.f5.get_min_tc_ratio()

    def test_get_max_tc_ratio(self):
        self.f5.get_max_tc_ratio()

    def test_get_tc_ratio_limits(self):
        self.f5.get_tc_ratio_limits() 

    def test_get_basecaller_min_events(self):
        if self.f5._basecalling_attempted:
            self.f5.get_basecaller_min_events()

    def test_get_basecaller_max_events(self):
            self.f5.get_basecaller_max_events()

    def test_get_basecaller_events_limits(self):
            self.f5.get_basecaller_events_limits()


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    fast5class_testnames = loader.getTestCaseNames(TestFast5Class)
    for d in data_dirs:
        f = os.path.join(data_path, d, fast5filename)
        for test in fast5class_testnames:
            suite.addTest(TestFast5Class(test, filename=f))
    return suite
