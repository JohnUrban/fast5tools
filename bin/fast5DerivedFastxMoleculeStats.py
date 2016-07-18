#!/usr/bin/env python

import os, sys
from Bio import SeqIO
from fast5tools.fxclass import *
import argparse

## JOHN URBAN (2015,2016)
#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5-derived fastx file(s), return tab-delimited information for each molecule from each file.

1 = molecule_name
2 = molecule length
3 = has complement
4 = has 2d
5 = 2d seq len
6 = template seq len
7 = complement seq len
8 = 2d mean q score
9 = template mean q score
10 = filename

CAUTION:
Using fast5tofastx.py, it is possible to make fast5-derived fastx files that contain only specific read types.
For example, one may generate a file containing only 2D reads.
In that has has_complement is meaningless 
-->it did have a complement by default, but that will not be detected from the fastx file.
Moreover, since this was designed to give the same first 9 columns as the standard option for fast5stats.py,
it will never say whether or not templates were in the file.
All fast5 files that have any sequence information in them have template by default
-- so whether they were found in the file or not is not what is being interrogated here.


NOTE -- obtaining molecule stats vs read stats from fast5-derived fastx files:
If you are only interested in sequence lengths and/or mean quality scores on a per read basis,
fast5DerivedFastxMoleculeStats.py is not the appropriate tool/approach.
It should be used to translate a fast5-derived fastx file into information about each individual
fast5/molecule that those reads came from.
Recall each fast5 has a template, may have a complement, and if so, may also have a 2d.
Thus, this is best used on fastx files that include ALL reads from each fast5 made by:
fast5tofastx.py -r all

To get just seqlen and quality info for individual reads in a fast5-derived fastx file, all the info
you want is already in the read names. Using grep and awk is a better approach. For example, for mol_name, readtype, len, Q:

grep ">" file.fasta | awk 'OFS="\t" {sub(/>/,""); gsub(/\|/,"\t"); gsub(/:/,"\t"); print $11"_"$13"_"$7"_"$9, $1, $3, $5}'

I have added 'fast5DerivedFastaReadStats.sh' to the tools in the bin directory that does just this.

John Urban (2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fastx', metavar='fastx', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')


parser.add_argument('-d', '--delimiter', type=str, default='\t',
                    help='''Provide choice of delimiter. Default: tab.
For new line, type 'newline'.
Tab is default, but can also specify with "tab"''')


parser.add_argument('-f', '--filename', action='store_true', default=False,
                    help='''Append fastx filename each molecule was from to end of line (new column).''')


parser.add_argument('-i', '--intype', type=str, default='fasta',
                    help=''' Choices: fasta, fastq, input.
Default: fasta.
Note: input (one or both formats found in input).
Declare which input types are to be explored for filtering.
Since one can direct this script at directories that may contain both fasta and fastq,
this gives an extra level of awareness to explore only a given file type (or both).
One may also want to look at both fasta and fastq, but output only fasta (see -o).
In rare cases, one may want to read in both types, and return the same type (-i input, -o input).
For now all output is directed to stdout, so the latter case is not recommended.
In the future, if output from each given input file can automatically be directed to a similarly named
output file (with .filtered. added in, then it might make more sense.''')



args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################



if args.intype == 'input' or args.intype == "both":
    intypes = ['fasta', 'fastq']
elif args.intype == 'fasta':
    intypes = ['fasta']
elif args.intype == 'fastq':
    intypes = ['fastq']

readtype = "molecule"
minlen = 0
maxlen = float("inf")
minq = 0
maxq = float("inf")
channel = None
readnum = None
asic = None
runid = None
deviceid = None
modelid=None
minscore  = None
filter_rule = "and"

for fxfile in FastXFileList(args.fastx, intypes=intypes):
    fxmol = None
    if args.filename:
        filename = "\t" + fxfile.filename
    else:
        filename = ""
    for fx in fxfile:
        #if not initiated or current molecule not same as previous molecule, start new molecule
        ## elif initiated but new entry is from new molecule, print desired entry from previous molecule and start new molecule
        ## elif initiated and new entry is from current molecule, add it to the molecule info
        if fxmol == None:
            fxmol = FastXMolecule(fx)
        elif fx.get_molecule_name() != fxmol.get_molecule_name():
            rtype = fxmol.interpret(readtype)
            if fxmol.passes_filter(rtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule):
                print fxmol.get_molecule_stats_string() + filename
            fxmol = FastXMolecule(fx)
        elif fx.get_molecule_name() == fxmol.get_molecule_name():
            fxmol.add_fx(fx)
        else:
            quit("Error in code or file.....")
    ## process last molecule
    rtype = fxmol.interpret(readtype)
    if fxmol.passes_filter(rtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule):
        print fxmol.get_molecule_stats_string() + filename
            
