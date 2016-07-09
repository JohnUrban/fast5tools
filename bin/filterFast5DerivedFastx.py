#!/usr/bin/env python

import os, sys
from Bio import SeqIO
import argparse
from fast5tools.fxclass import *

#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fasta/fastq file(s) and/or directories containing them,

return sequences filtered by any of the following:
1. read type (2d, template, complement)
2. length
3. mean quality score
4. channel number
5. read number
6. ASIC ID
7. Run ID
8. Device ID
9. Base-calling model

Return in same format as input or choose:
fasta, fastq, qual, intqual, oldfalcon, newfalcon

Note: newfalcon output is fasta with falcon-compatible headers.
For newer Falcon:
>asic_run_device_basecallingmodel/i/0_readlen OtherInfo
Where
i is order it was encountered in
OtherInfo will include readtype,mean quality score, read number, channel number

For oldfalcon output:
>000_000/i/0_readlen OriginalFastaHeader
Where i is number read is encountered in.


TODO: fastaqual, fastaintqual

NOTE:
Fasta can be converted to fastq or quals, BUT the quals will not be correct per se.
First, they will be related to the mean q-score (Q).
Second, they will be rounded to the nearest int.
Thus, the mean q-score in the header/seqname will not be consistent with the mean of the quality scores.
It is related by int(round(Q)).

For now, input files are fasta, fastq, or dirs with them.
TODO: Allow tar/tarlite approach. Allow gzipped. Allow FOFN.


TODO: falconizeFast5DerivedFastx.py for more options and more description/info.

John Urban (2015, 2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fastx', metavar='fastx', nargs='+',
                   type= str, 
                   help='''Paths to as many fasta and/or fastq files and/or directories filled with them as you want.
Assumes all fasta files have '.fasta', '.fa', and/or '.fna' extensions (only accepts these).
Assumes all fastq files have '.fastq' or '.fq' extensions (only accepts these).
Assumes given that given one of the above extensions, the internal formatting is consistent with either fasta or fastq.
If inside dir of dirs with desired files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="all",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d', 'molecule', 'all', 'MoleQual'.
Default: all.
There is no need to write full word for options - can do: t, c, 2, m, a, M.
Molecule returns single fastx for when there is more than 1 record for a given Channel#/Read#/Run_ID/ASIC_ID:
if 2d present, return 2d.
elif complement present with no 2d, return longer of template or complement.
elif only template present, return template.
'MoleQual' is similar to molecule.
It differs only in choosing between template and complement when a 2D is not present.
Instead of choosing the longer one, it chooses the one with a higher quality mean quality score.

NOTE: it is assumed that reads from same molecule (shared Channel#/Read#/Run_ID/ASIC_ID)
are clustered together (i.e. occur consecutively) in given input.
If not, then molecule and MoleQual protocols will not work as expected.''')

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

parser.add_argument('-o', '--outtype', type=str, default='fasta',
                    help = '''Choices: input, fasta, fastq, qual, intqual, falcon.
Default: fasta.
Note: input (whatever format the file comes in as).
See -i for discussion on use cases.
falcon: returns fasta but with fasta headers compatible with FALCON assembler.
TODO:
fastaqual, fastaintqual''')

parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')

parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')

parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')

parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')

parser.add_argument('--channel', type=str, default=None, help='''Only report reads from given channel number. Default: reports from any/all channels present.''')
parser.add_argument('--readnum', type=str, default=None, help='''Only report reads with given read number. Default: reports from any/all read numbers.''')
parser.add_argument('--asic', type=str, default=None, help='''Only report reads with given asic ID. Default: reports from any/all ASIC IDs present.''')
parser.add_argument('--run', type=str, default=None, help='''Only report reads with given run ID. Default: reports from any/all Run IDs present.''')
parser.add_argument('--device', type=str, default=None, help='''Only report reads with given device ID. Default: reports from any/all Device IDs present.''')
parser.add_argument('--model', type=str, default=None, help='''Only report reads with given bas-calling model ID. Default: reports from any/all basecalling IDs present.''')

parser.add_argument('--rule', type=str, default='and', help='''Require each sequence to pass ALL the filters (use 'and') or pass at least N filters (use 'or')''')
parser.add_argument('--minscore', type=int, default=1, help='''If requiring sequences only pass at least N filters (--rule 'or'), then specify minimum number of filters to pass. Default: 1.''')



##parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
##The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
##The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
##The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
##Tarlite may become the default method after some testing if it performs at similar speeds.''')


args = parser.parse_args()


#################################################
## deal with some of the arguments
#################################################
assert args.outtype in ("fasta", "fastq", "qual", "intqual", "falcon", "oldfalcon", "newfalcon")
assert args.intype in ("input", "fasta", "fastq", "both")


assert args.readtype[0] in "tc2maM"
if args.readtype[0] == "t":
    args.readtype = "template"
elif args.readtype[0] == "c":
    args.readtype = "complement"
elif args.readtype[0] == "2":
    args.readtype = "2d"
elif args.readtype[0] == "m":
    args.readtype = "molecule"
elif args.readtype[0] == "a":
    args.readtype = "all"
elif args.readtype[0] == "M":
    args.readtype = "MoleQual"

    
if args.intype == 'input' or args.intype == "both":
    intypes = ['fasta', 'fastq']
elif args.intype == 'fasta':
    intypes = ['fasta']
elif args.intype == 'fastq':
    intypes = ['fastq']








def filter_by_entry(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule, intypes, outtype):
    ## if want all or specific read type, just filter by entry
    rtype = readtype
    falcon_i = 0
    for fxfile in FastXFileList(args.fastx, intypes=intypes):
        for fx in fxfile:
            falcon_i += 1
            if readtype == "all":
                rtype = fx.get_read_type()
            if fx.passes_filter(rtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule):
                print fx.get_fastx_entry(outtype, falcon_i)
                

def filter_by_molecule(readtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule, intypes, outtype):
    ## if filtering by molecule or molequal; else just filter by entry
    falcon_i = 0
    for fxfile in FastXFileList(args.fastx, intypes=intypes):
        fxmol = None
        for fx in fxfile:
            falcon_i += 1
            #if not initiated or current molecule not same as previous molecule, start new molecule
            ## elif initiated but new entry is from new molecule, print desired entry from previous molecule and start new molecule
            ## elif initiated and new entry is from current molecule, add it to the molecule info
            if fxmol == None:
                fxmol = FastXMolecule(fx)
            elif fx.get_molecule_name() != fxmol.get_molecule_name():
                rtype = fxmol.interpret(readtype)
                if fxmol.passes_filter(rtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule):
                    print fxmol.get_fastx_entry(rtype, outtype, falcon_i)
                fxmol = FastXMolecule(fx)
            elif fx.get_molecule_name() == fxmol.get_molecule_name():
                fxmol.add_fx(fx)
            else:
                quit("Error in code or file.....")
        ## process last molecule
        rtype = fxmol.interpret(readtype)
        if fxmol.passes_filter(rtype, minlen, maxlen, minq, maxq, channel, readnum, asic, runid, deviceid, modelid, minscore, filter_rule):
            print fxmol.get_fastx_entry(rtype, outtype, falcon_i)
            
if args.readtype in ('template', '2d', 'complement', 'all'):
    filter_by_entry(readtype=args.readtype, minlen=args.minlen, maxlen=args.maxlen, minq=args.minq, maxq=args.maxq, channel=args.channel, readnum=args.readnum, asic=args.asic, runid=args.run, deviceid=args.device, modelid=args.model, minscore=args.minscore, filter_rule=args.rule, intypes=intypes, outtype=args.outtype)
elif args.readtype in ('molecule', 'MoleQual'):
    filter_by_molecule(readtype=args.readtype, minlen=args.minlen, maxlen=args.maxlen, minq=args.minq, maxq=args.maxq, channel=args.channel, readnum=args.readnum, asic=args.asic, runid=args.run, deviceid=args.device, modelid=args.model, minscore=args.minscore, filter_rule=args.rule, intypes=intypes, outtype=args.outtype)



