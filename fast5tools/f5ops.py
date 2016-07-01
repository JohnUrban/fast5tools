## JOHN URBAN (2015, 2016)

import os
######################### processing functions ######
#### e.g. used in:
#### fast5tofastx.py, 
#################################################

def get_single_read(f5, readtype, minlen, maxlen, minq, maxq, output):
    ''' f5 is Fast5 object.
        readtype in template, complement, 2d, molecule, all.
        minlen/maxlen - integers.
        minq/maxq - floats
        output is a function'''
    if f5.has_read(readtype):
        if f5.get_seq_len(readtype) >= minlen and f5.get_seq_len(readtype) <= maxlen:
            if f5.get_mean_qscore(readtype) >= minq and f5.get_mean_qscore(readtype) <= maxq:
                return output(f5, readtype)

def get_template_read(f5, minlen, maxlen, minq, maxq, output):
    return get_single_read(f5, "template", minlen, maxlen, minq, maxq, output)

def get_complement_read(f5, minlen, maxlen, minq, maxq, output):
    return get_single_read(f5, "complement", minlen, maxlen, minq, maxq, output)

def get_2d_read(f5, minlen, maxlen, minq, maxq, output):
    return get_single_read(f5, "2d", minlen, maxlen, minq, maxq, output)

def get_molecule_read(f5, minlen, maxlen, minq, maxq, output):
    return get_single_read(f5, f5.use_molecule(), minlen, maxlen, minq, maxq, output)

def get_all_reads(f5, minlen, maxlen, minq, maxq, output):
    allreads = ""
    for readtype in ["template", "complement", "2d"]:
        try:
            allreads += get_single_read(f5, readtype, minlen, maxlen, minq, maxq, output) + '\n'
        except:
            pass
    return allreads.rstrip()



######################### output functions ######
#### e.g. used in get_single_read()
#################################################

def fastq(f5, readtype):
    return f5.get_fastq(readtype)

def fasta(f5, readtype):
    return f5.get_fasta(readtype)

def qual(f5, readtype):
    return f5.get_quals(readtype)

def intqual(f5, readtype):
    return f5.get_quals_as_int(readtype)





#######
f5fxn = {}
f5fxn[1] = lambda f5: f5.get_base_info_name()
f5fxn[2] = lambda f5: f5.get_seq_len(f5.use_molecule()) if f5.has_reads() else "-"
f5fxn[3] = lambda f5: 1 if f5.has_read("complement") else 0
f5fxn[4] = lambda f5: 1 if f5.has_read("2d") else 0
f5fxn[5] = lambda f5: f5.get_seq_len("2d") if f5.has_read("2d") else "-"
f5fxn[6] = lambda f5: f5.get_seq_len("template") if f5.has_read("template") else "-"
f5fxn[7] = lambda f5: f5.get_seq_len("complement") if f5.has_read("complement") else "-"
f5fxn[8] = lambda f5: f5.get_mean_qscore("2d") if f5.has_read("2d") else "-"
f5fxn[9] = lambda f5: f5.get_mean_qscore("template") if f5.has_read("template") else "-"
f5fxn[10] = lambda f5: f5.get_mean_qscore("complement") if f5.has_read("complement") else "-"
f5fxn[11] = lambda f5: f5.get_num_events("input")
f5fxn[12] = lambda f5: f5.get_num_events("template") if f5.has_read("template") else "-"
f5fxn[13] = lambda f5: f5.get_num_events("complement") if f5.has_read("complement") else "-"
f5fxn[14] = lambda f5: f5.get_num_called_events("template") if f5.has_read("template") else "-"
f5fxn[15] = lambda f5: f5.get_num_called_events("complement") if f5.has_read("complement") else "-"
f5fxn[16] = lambda f5: f5.get_num_skips("template") if f5.has_read("template") else "-"
f5fxn[17] = lambda f5: f5.get_num_skips("complement") if f5.has_read("complement") else "-"
f5fxn[18] = lambda f5: f5.filename
f5fxn[19] = lambda f5: f5.abspath
f5fxn[20] = lambda f5: f5.get_log()
f5fxn[21] = lambda f5: f5.get_log_string()
##f5fxn[19] = lambda f5
##f5fxn[20] = lambda f5
##f5fxn[21] = lambda f5
##f5fxn[22] = lambda f5
##f5fxn[23] = lambda f5
##f5fxn[24] = lambda f5













