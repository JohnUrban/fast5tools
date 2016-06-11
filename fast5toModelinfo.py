#!/usr/bin/env python

import h5py, os, sys
import cStringIO as StringIO
from Bio import SeqIO
import f5class
from f5ops import *
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
Choices: 'template', 'complement'.
Default: input.
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



args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
assert args.readtype[0] in "tc"
if args.readtype[0] == "t":
    args.readtype = "template"
elif args.readtype[0] == "c":
    args.readtype = "complement"
elif args.readtype[0] == "i":
    args.readtype = "input"
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
    if len(f5list) == 1:
        for f5 in f5list:
            if f5.has_read(args.readtype):
                print get_fxn(f5,args.readtype)
    elif len(f5list) > 1:
        for f5 in f5list:
            if f5.has_read(args.readtype):
                out = open(args.outdir + f5.get_base_info_name() + "." + args.readtype + "." + args.get + ".txt", 'w')
                out.write(get_fxn(f5,args.readtype))
                out.close()

                
def get_model_type(f5list):
    for f5 in f5list:
        print ("\t").join([ f5.get_base_info_name(), f5.get_model_type() ])



#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    f5list = f5class.Fast5List(args.fast5)
##    for f5 in f5list:
##        print f5.get_model_attrs_string("template")
##    quit()
    if args.get == "model" or args.get == "attr":
        get_model(f5list, args, get_fxn)
    elif args.get == "type":
        get_model_type(f5list)

        


