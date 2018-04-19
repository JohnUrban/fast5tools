## JOHN URBAN (2015, 2016)

import os
from random import shuffle, seed
from fast5tools.f5class import *

##
#### fast5tofastx.py,
def assert_readtype(readtype, legaloptions="tc2maM"):
    allowed = readtype[0] in legaloptions
    assert allowed
    if readtype[0] == "t" and allowed:
        readtype = "template"
    elif readtype[0] == "c" and allowed:
        readtype = "complement"
    elif readtype[0] == "2" and allowed:
        readtype = "2d"
    elif readtype[0] == "m" and allowed:
        readtype = "molecule"
    elif readtype[0] == "a" and allowed:
        readtype = "all"
    elif readtype[0] == "M" and allowed:
        readtype = "MoleQual"
    return readtype

#### fast5tofastx.py,
def details(f5, readtype):
    readstats = []
    readstats.append( f5._get_pore_info_name(readtype) )
    readstats.append( f5.get_seq_len(readtype) )
    readstats.append( f5.get_mean_qscore(readtype) )
    readstats.append( f5.get_num_events(readtype) )
    try:
        readstats.append( f5.get_num_called_events(readtype) )
    except:
        readstats.append("-")
    try:
        readstats.append( f5.get_num_skips(readtype) )
    except:
        readstats.append("-")
    try:
        readstats.append( f5.get_num_stays(readtype) )
    except:
        readstats.append("-")        
    return ("\t").join([str(e) for e in readstats])


#### fast5tofastx.py,
def get_fast5tofastx_output_fxn(outtype):
    if outtype == "fasta":
        output = fasta
    elif outtype == "fastq":
        output = fastq
    elif outtype == "qual":
        output = qual
    elif outtype == "intqual":
        output = intqual
    elif outtype == "details":
        output = details
    elif outtype == "falcon" or outtype == "oldfalcon":
        output = oldfalcon
    elif outtype == "newfalcon":
        output = newfalcon
    elif outtype == "fasta_with_abspath":
        output = fasta_with_abspath
    elif outtype == "fasta_only_abspath":
        output = fasta_only_abspath
    elif outtype == "fastq_with_abspath":
        output = fastq_with_abspath
    elif outtype == "fastq_only_abspath":
        output = fastq_only_abspath
    elif outtype == "qual_with_abspath":
        output = qual_with_abspath
    elif outtype == "qual_only_abspath":
        output = qual_only_abspath
    #
    elif outtype == "fasta_with_filename":
        output = fasta_with_filename
    elif outtype == "fasta_only_filename":
        output = fasta_only_filename
    elif outtype == "fastq_with_filename":
        output = fastq_with_filename
    elif outtype == "fastq_only_filename":
        output = fastq_only_filename
    elif outtype == "qual_with_filename":
        output = qual_with_filename
    elif outtype == "qual_only_filename":
        output = qual_only_filename
    #
    elif outtype == "fasta_readstatsname":
        output = fasta_readstatsname
    elif outtype == "fasta_readstatsname_with_abspath":
        output = fasta_readstatsname_with_abspath
    elif outtype == "fasta_readstatsname_with_filename":
        output = fasta_readstatsname_with_filename
        
    elif outtype == "fastq_readstatsname":
        output = fastq_readstatsname
    elif outtype == "fastq_readstatsname_with_abspath":
        output = fastq_readstatsname_with_abspath
    elif outtype == "fastq_readstatsname_with_filename":
        output = fastq_readstatsname_with_filename
    #
    elif outtype == "qual_readstatsname":
        output = qual_readstatsname
    elif outtype == "qual_readstatsname_with_abspath":
        output = qual_readstatsname_with_abspath
    elif outtype == "qual_readstatsname_with_filename":
        output = qual_readstatsname_with_filename
    #
    return output

def get_fast5tofastx_readtype_fxn(readtype):
    if readtype == "template":
        getread = get_template_read
    elif readtype == "complement":
        getread = get_complement_read
    elif readtype == "2d":
        getread = get_2d_read
    elif readtype == "molecule":
        getread = get_molecule_read
    elif readtype == "all":
        getread = get_all_reads
    elif readtype == "MoleQual":
        getread = get_molequal_read
    return getread

def get_fast5tofastx_fxns(outtype, readtype):
    ''' '''
    ### get outtype fxn ###
    output = get_fast5tofastx_output_fxn(outtype)
##    if args.outtype == "fasta":
##        output = fasta
##    elif args.outtype == "fastq":
##        output = fastq
##    elif args.outtype == "qual":
##        output = qual
##    elif args.outtype == "intqual":
##        output = intqual
##    elif args.outtype == "details":
##        output = details
##    elif args.outtype == "falcon" or args.outtype == "oldfalcon":
##        output = oldfalcon
##    elif args.outtype == "newfalcon":
##        output = newfalcon
##    elif args.outtype == "fasta_with_abspath":
##        output = fasta_with_abspath
##    elif args.outtype == "fasta_only_abspath":
##        output = fasta_only_abspath
##    elif args.outtype == "fastq_with_abspath":
##        output = fastq_with_abspath
##    elif args.outtype == "fastq_only_abspath":
##        output = fastq_only_abspath
##    elif args.outtype == "qual_with_abspath":
##        output = qual_with_abspath
##    elif args.outtype == "qual_only_abspath":
##        output = qual_only_abspath
##    #
##    elif args.outtype == "fasta_with_filename":
##        output = fasta_with_filename
##    elif args.outtype == "fasta_only_filename":
##        output = fasta_only_filename
##    elif args.outtype == "fastq_with_filename":
##        output = fastq_with_filename
##    elif args.outtype == "fastq_only_filename":
##        output = fastq_only_filename
##    elif args.outtype == "qual_with_filename":
##        output = qual_with_filename
##    elif args.outtype == "qual_only_filename":
##        output = qual_only_filename
##    #
##    elif args.outtype == "fasta_readstatsname":
##        output = fasta_readstatsname
##    elif args.outtype == "fasta_readstatsname_with_abspath":
##        output = fasta_readstatsname_with_abspath
##    elif args.outtype == "fasta_readstatsname_with_filename":
##        output = fasta_readstatsname_with_filename
##        
##    elif args.outtype == "fastq_readstatsname":
##        output = fastq_readstatsname
##    elif args.outtype == "fastq_readstatsname_with_abspath":
##        output = fastq_readstatsname_with_abspath
##    elif args.outtype == "fastq_readstatsname_with_filename":
##        output = fastq_readstatsname_with_filename
##        
##    elif args.outtype == "qual_readstatsname":
##        output = qual_readstatsname
##    elif args.outtype == "qual_readstatsname_with_abspath":
##        output = qual_readstatsname_with_abspath
##    elif args.outtype == "qual_readstatsname_with_filename":
##        output = qual_readstatsname_with_filename
        
    ### get readtype fxn ###
    getread = get_fast5tofastx_readtype_fxn(readtype)
##    if args.readtype == "template":
##        getread = get_template_read
##    elif args.readtype == "complement":
##        getread = get_complement_read
##    elif args.readtype == "2d":
##        getread = get_2d_read
##    elif args.readtype == "molecule":
##        getread = get_molecule_read
##    elif args.readtype == "all":
##        getread = get_all_reads
##    elif args.readtype == "MoleQual":
##        getread = get_molequal_read
    return output, getread



######################### processing functions ######
#### e.g. used in:
#### fast5tofastx.py, 
#################################################
def get_comments(request, f5, readtype, samflag=''):
    if not request:
        return False
    elif request in ('base_info', 'pore_info', 'read_stats', 'event_stats', 'read_event_stats', 'abspath', 'filename'):
        if request == 'base_info':
            return samflag + f5.get_base_info_name()
        elif request == 'pore_info':
            return samflag + f5.get_pore_info_name(readtype)
        elif request == 'read_stats':
            return samflag + f5.get_read_stats_name(readtype)
        elif request == 'event_stats':
            return samflag + f5.get_event_stats_name(readtype)
        elif request == 'read_event_stats':
            return samflag + f5.get_read_and_event_stats_name(readtype)
        elif request == 'abspath':
            return samflag + f5.abspath
        elif request == 'filename':
            return samflag + f5.filebasename
    else:
        return samflag + request
    
def get_single_read(f5, readtype, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    ''' f5 is Fast5 object.
        readtype in template, complement, 2d, molecule, all.
        minlen/maxlen - integers.
        minq/maxq - floats
        output is a function'''
    if f5.has_read(readtype):
        if f5.get_seq_len(readtype) >= minlen and f5.get_seq_len(readtype) <= maxlen:
            if f5.get_mean_qscore(readtype) >= minq and f5.get_mean_qscore(readtype) <= maxq:
                if kwargs['comments']:
                    kwargs['comments'] = get_comments(kwargs['comments'], f5, readtype, kwargs['samflag'])
                return output(f5, readtype, *args, **kwargs)

def get_template_read(f5, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    return get_single_read(f5, "template", minlen, maxlen, minq, maxq, output, *args, **kwargs)

def get_complement_read(f5, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    return get_single_read(f5, "complement", minlen, maxlen, minq, maxq, output, *args, **kwargs)

def get_2d_read(f5, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    return get_single_read(f5, "2d", minlen, maxlen, minq, maxq, output, *args, **kwargs)

def get_molecule_read(f5, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    return get_single_read(f5, f5.use_molecule(), minlen, maxlen, minq, maxq, output, *args, **kwargs)

def get_molequal_read(f5, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    return get_single_read(f5, f5.use_molequal(), minlen, maxlen, minq, maxq, output, *args, **kwargs)


def get_all_reads(f5, minlen, maxlen, minq, maxq, output, *args, **kwargs):
    allreads = ""
    for readtype in ["template", "complement", "2d"]:
        try:
            allreads += get_single_read(f5, readtype, minlen, maxlen, minq, maxq, output, *args, **kwargs) + '\n'
        except:
            pass
    return allreads.rstrip()



######################### output functions ######
#### e.g. used in get_single_read()
#################################################


def fastq(f5, readtype, *args, **kwargs):
    return f5.get_fastq(readtype, comments=kwargs['comments'])

def fastq_with_abspath(f5, readtype, *args, **kwargs):
    return f5.get_fastq_with_abspath(readtype, comments=kwargs['comments'])

def fastq_only_abspath(f5, readtype, *args, **kwargs):
    return f5.get_fastq_only_abspath(readtype, comments=kwargs['comments'])

def fastq_with_filename(f5, readtype, *args, **kwargs):
    return f5.get_fastq_with_filename(readtype, comments=kwargs['comments'])

def fastq_only_filename(f5, readtype, *args, **kwargs):
    return f5.get_fastq_only_filename(readtype, comments=kwargs['comments'])

def fastq_readstatsname(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name(readtype)
    return f5.get_fastq(readtype, name = name, comments=kwargs['comments'])

def fastq_readstatsname_with_filename(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name_with_filebasename(readtype)
    return f5.get_fastq(readtype, name = name, comments=kwargs['comments'])

def fastq_readstatsname_with_abspath(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name_with_abspath(readtype)
    return f5.get_fastq(readtype, name = name, comments=kwargs['comments'])




def fasta(f5, readtype, *args, **kwargs):
    return f5.get_fasta(readtype, comments=kwargs['comments'])

def fasta_with_abspath(f5, readtype, *args, **kwargs):
    return f5.get_fasta_with_abspath(readtype, comments=kwargs['comments'])

def fasta_only_abspath(f5, readtype, *args, **kwargs):
    return f5.get_fasta_only_abspath(readtype, comments=kwargs['comments'])

def fasta_with_filename(f5, readtype, *args, **kwargs):
    return f5.get_fasta_with_filename(readtype, comments=kwargs['comments'])

def fasta_only_filename(f5, readtype, *args, **kwargs):
    return f5.get_fasta_only_filename(readtype, comments=kwargs['comments'])

def fasta_readstatsname(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name(readtype)
    return f5.get_fasta(readtype, name = name, comments=kwargs['comments'])

def fasta_readstatsname_with_filename(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name_with_filebasename(readtype)
    return f5.get_fasta(readtype, name = name, comments=kwargs['comments'])

def fasta_readstatsname_with_abspath(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name_with_abspath(readtype)
    return f5.get_fasta(readtype, name = name, comments=kwargs['comments'])



def qual(f5, readtype, *args, **kwargs):
    return f5.get_quals(readtype, comments=kwargs['comments'])

def qual_with_abspath(f5, readtype, *args, **kwargs):
    return f5.get_quals_with_abspath(readtype, comments=kwargs['comments'])

def qual_only_abspath(f5, readtype, *args, **kwargs):
    return f5.get_quals_only_abspath(readtype, comments=kwargs['comments'])

def qual_with_filename(f5, readtype, *args, **kwargs):
    return f5.get_quals_with_filename(readtype, comments=kwargs['comments'])

def qual_only_filename(f5, readtype, *args, **kwargs):
    return f5.get_quals_only_filename(readtype, comments=kwargs['comments'])


def qual_readstatsname(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name(readtype)
    return f5.get_quals(readtype, name = name, comments=kwargs['comments'])

def qual_readstatsname_with_filename(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name_with_filebasename(readtype)
    return f5.get_quals(readtype, name = name, comments=kwargs['comments'])

def qual_readstatsname_with_abspath(f5, readtype, *args, **kwargs):
    name = f5.get_read_stats_name_with_abspath(readtype)
    return f5.get_quals(readtype, name = name, comments=kwargs['comments'])


def intqual(f5, readtype, *args, **kwargs):
    return f5.get_quals_as_int(readtype)

def oldfalcon(f5, readtype, *args, **kwargs):
    return f5.get_falcon_fasta(readtype, zmw_num=kwargs['falcon_i'], style="old")

def newfalcon(f5, readtype, *args, **kwargs):
    return f5.get_falcon_fasta(readtype, zmw_num=kwargs['falcon_i'], style="new")

#######
f5fxn = {}
f5fxn[1] = lambda f5: f5.get_base_info_name()
f5fxn[2] = lambda f5: f5.get_seq_len(f5.use_molecule()) if f5.has_reads() else "-"
f5fxn[3] = lambda f5: f5.get_mean_qscore(f5.use_molecule()) if f5.has_reads() else "-"
f5fxn[4] = lambda f5: 1 if f5.has_read("complement") else 0
f5fxn[5] = lambda f5: 1 if f5.has_read("2d") else 0
f5fxn[6] = lambda f5: f5.get_seq_len("2d") if f5.has_read("2d") else "-"
f5fxn[7] = lambda f5: f5.get_seq_len("template") if f5.has_read("template") else "-"
f5fxn[8] = lambda f5: f5.get_seq_len("complement") if f5.has_read("complement") else "-"
f5fxn[9] = lambda f5: f5.get_mean_qscore("2d") if f5.has_read("2d") else "-"
f5fxn[10] = lambda f5: f5.get_mean_qscore("template") if f5.has_read("template") else "-"
f5fxn[11] = lambda f5: f5.get_mean_qscore("complement") if f5.has_read("complement") else "-"
f5fxn[12] = lambda f5: f5.get_num_events("input")
f5fxn[13] = lambda f5: f5.get_num_events("template") if f5.has_read("template") else "-"
f5fxn[14] = lambda f5: f5.get_num_events("complement") if f5.has_read("complement") else "-"
f5fxn[15] = lambda f5: f5.get_num_called_events("template") if f5.has_read("template") else "-"
f5fxn[16] = lambda f5: f5.get_num_called_events("complement") if f5.has_read("complement") else "-"
f5fxn[17] = lambda f5: f5.get_num_skips("template") if f5.has_read("template") else "-"
f5fxn[18] = lambda f5: f5.get_num_skips("complement") if f5.has_read("complement") else "-"
f5fxn[19] = lambda f5: f5.filename
f5fxn[20] = lambda f5: f5.abspath
f5fxn[21] = lambda f5: f5.get_log()
f5fxn[22] = lambda f5: f5.get_log_string()
##f5fxn[19] = lambda f5
##f5fxn[20] = lambda f5
##f5fxn[21] = lambda f5
##f5fxn[22] = lambda f5
##f5fxn[23] = lambda f5
##f5fxn[24] = lambda f5



## F5 FILTERING/SELECTION OPERATIONS
def define_read_type(f5, readtype):
        if readtype in ('molecule'):
            return f5.use_molecule()
        else:
            return readtype

def meets_all_criteria(f5, readtype, minlen, maxlen, minq, maxq):
    if f5.is_not_corrupt() and f5.is_nonempty():
        readtype = define_read_type(f5, readtype)
        if f5.has_read(readtype):
            if f5.get_seq_len(readtype) >= minlen and f5.get_seq_len(readtype) <= maxlen:
                if f5.get_mean_qscore(readtype) >= minq and f5.get_mean_qscore(readtype) <= maxq:
                    return True
    else:
        return False


## Fast5List OPs
## The function below allows one to sample from the Fast5List
## get_fast5_list(args.nfiles, args.fast5, args.random, args.randomseed, args.notarlite, filemode='r', sort=True)
def get_fast5_list(nfiles, initial_list, random=False, randomseed=False, notarlite=False, filemode='r', sort=True):
    if nfiles <= 0:
        nfiles = len(initial_list)
    if random:
        if randomseed:
            seed(randomseed) ## This seed will only make things reproducible given same exact conditions - seed not re-used later.
        shuffle(initial_list) ## This only shuffles initial targets, not final files -- but can help simplify target expansion

    # Downsample as necessary
    initial_list_ds = initial_list[:nfiles] ## This only shrinks number of targets to simplify target expansion
    
    # Expand targets to get initial list
    f5list = Fast5List(initial_list_ds, keep_tar_footprint_small=(not notarlite), filemode=filemode)

    ## Take from expanded list: Shuffling and Downsampling actually happens here on indiv fast5s
    return f5list.get_sample(n=nfiles, random=random, sort=sort)






