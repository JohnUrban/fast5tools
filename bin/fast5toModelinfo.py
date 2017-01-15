#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
from fast5tools.f5class import *
from fast5tools.f5ops import *
import argparse
from glob import glob


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s, return info on base-calling model, its parameters, etc.

John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="template",
                   type= str, 
                   help='''Choose type of model to get.
Choices: 'template', 'complement', 'both'. (both is only for --get attr)
Default: template.
''')

parser.add_argument('-g', '--get', default="model",
                   type= str, 
                   help='''Get model, model attributes, or model type Specify:
model, attr, or type.
Default: model.
''')

parser.add_argument('-o', '--outdir', type=str, default="",
                    help = '''If single fast5 specified, it will be reported to stdout.
If multiple fast5s are specified, they will be saved to files in the working dir by default.
This flag allows you to specify a different output directory.
If model or attr specified to --get, filenames will be the the base info name of the fast5 file with .model.txt or .modelattr.txt appended.
If 'type' specified for --get, tab-delimited "base-info-name, modeltype" will be output to stdout no matter how many files.
''')

parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
Tarlite may become the default method after some testing if it performs at similar speeds.''')


args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
assert args.readtype[0] in "tcb"
if args.readtype == 'both':
    assert args.get == 'attr'
    
if args.readtype[0] == "t":
    args.readtype = "template"
elif args.readtype[0] == "c":
    args.readtype = "complement"
##elif args.readtype[0] == "i":
##    args.readtype = "input"
if args.outdir:
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    if args.outdir[-1] != "/":
        args.outdir += "/"



if args.get == "model":
    get_fxn = lambda f5, readtype: f5.get_model_string(readtype)
    
elif args.get == "attr":
    get_fxn = lambda f5, readtype: f5.get_model_attrs_string(readtype)
    
def get_model(f5list, args, get_fxn):
    if args.readtype == "both":
        readtypes = ("template", "complement")
    else:
        readtypes = (args.readtype,)
    if len(f5list) == 1:
        for f5 in f5list:
            for rt in readtypes:
                if f5.has_read(rt):
                    print get_fxn(f5,rt)
    elif len(f5list) > 1:
        for f5 in f5list:
            for rt in readtypes:
                if f5.has_read(rt):
                    out = open(args.outdir + f5.filebasename + "." + rt + "." + args.get + ".txt", 'w')
                    out.write(get_fxn(f5,rt)+"\n")
                    out.close()

                
def get_model_type(f5list):
    for f5 in f5list:
        print ("\t").join([f5.filename, f5.get_base_info_name(), f5.get_model_type() ])



#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    f5list = Fast5List(args.fast5, keep_tar_footprint_small=args.tarlite)

##    for f5 in f5list:
##        print f5.get_model_attrs_string("template")
##    quit()
    if args.get == "model" or args.get == "attr":
        get_model(f5list, args, get_fxn)
    elif args.get == "type":
        get_model_type(f5list)

        


