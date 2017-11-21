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

Given path(s) to fast5 file(s) and/or directories of fast5s, return fasta, fastq, qual, or intqual for all fast5s found.

For files that are corrupt or empty, for now it silently skips them.
As an alternative, fast5stats will tell you all files skipped (in stderr or to specified file).

John Urban (2015, 2016, 2017)

TODO (11/17/2017): Allow customized name design with options for readtype, len, Q, channel, read num, asic, abspath, filename, etc

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="mol",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d', 'molecule', 'all', 'MoleQual'.
Default: molecule.
There is no need to write full word for options - can do: t, c, 2, m, a, M.
Molecule returns single fasta for each fast5 by following rules:
if 2d present, return 2d.
elif complement present with no 2d, return longer of template or complement.
elif only template present, return template.
'MoleQual' is similar to molecule.
It differs only in choosing between template and complement when a 2D is not present.
Instead of choosing the longer one, it chooses the one with a higher quality mean quality score.''')

parser.add_argument('-o', '--outtype', type=str, default="fasta",
                    help = '''Choices: fasta, fastq, qual, intqual, details, falcon, oldfalcon, newfalcon, fasta_readstatsname, fastq_readstatsname, qual_readstatsname.
Additional choices:
Add _with_abspath to fasta/fastq/qual options (and readstatname versions) to add absolute f5 file path to read name.
Add _with_filename to fasta/fastq/qual options (and readstatname versions) to add only the basename of each f5 file (excluding fast5 extension) to read name.
If only want abs path in name, add _only_abspath to fasta/fastq/qual options (and readstatname versions).
If only want basename of file in read name, add _only_filename to fasta/fastq/qualoptions (and readstatname versions).

Default: fasta.
If details, sequence not reported, but name, seqlen, and meanq are.
falcon/oldfalcon/newfalcon output fasta files that are compatible with FALCON assembler.
falcon and oldfalcon put out the same thing and might be safest choice as it should work with old and new FALCON versions.
newfalcon will only work with latest FALCON versions.
The real issue is fasta2DB, which is particular about fasta headers.
In older versions, it only allowed data from 1 SMRT cell per file.
Now it allows multiple SMRT cells per file,
only if all data from a given SMRT cell are grouped together.

NOTE: if 'all' is used, for now each will be given the same well number.
This could potentially have the side effect of using only the longest read in falcon,
if '-a' is not used in DBsplit.
To avoid this, just use '--outtype fasta', then use filterFast5DerivedFastx.py
to convert the nanopore fasta headers to falcon-compatible fasta headers.''')

parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')

parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')

parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')

parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')

parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
Tarlite may become the default method after some testing if it performs at similar speeds.''')

parser.add_argument('-c', '--comments', type=str, default=False, help='''Add comments to fastx names.
Comments are separated from main name (following > or @) with a tab.
Default: no comments/False.
Leave any desired string here (you may need to enclode in quotes is there are spaces).
Alternatively, specify one of the following options:
base_info
pore_info
read_stats
event_stats
read_event_stats

''')

parser.add_argument('-S', '--samflag', action='store_true', default=False, help='''Add sam flag to comments.

''')


args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################
legalouts = ("fasta", "fastq", "qual", "intqual", "details", "falcon", "oldfalcon", "newfalcon", "fasta_with_abspath", "fasta_only_abspath","fastq_with_abspath", "fastq_only_abspath", "qual_with_abspath", "qual_only_abspath", "fasta_with_filename", "fasta_only_filename","fastq_with_filename", "fastq_only_filename", "qual_with_filename", "qual_only_filename")
legalouts += ("fasta_readstatsname", "fasta_readstatsname_with_abspath", "fasta_readstatsname_with_filename")
legalouts += ("fastq_readstatsname", "fastq_readstatsname_with_abspath", "fastq_readstatsname_with_filename")
legalouts += ("qual_readstatsname", "qual_readstatsname_with_abspath", "qual_readstatsname_with_filename")
assert args.outtype in legalouts
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


#################################################
### uses output functions from f5ops.py
### fasta(), fastq(), qual(), intqual()
### 
#################################################



#################################################
## fast5tofastx specific "output" functions
#################################################

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



    
#################################################
####### argument processing functions ###########
#################################################

def get_fast5tofastx_fxns(args):
    ### get outtype fxn ###
    if args.outtype == "fasta":
        output = fasta
    elif args.outtype == "fastq":
        output = fastq
    elif args.outtype == "qual":
        output = qual
    elif args.outtype == "intqual":
        output = intqual
    elif args.outtype == "details":
        output = details
    elif args.outtype == "falcon" or args.outtype == "oldfalcon":
        output = oldfalcon
    elif args.outtype == "newfalcon":
        output = newfalcon
    elif args.outtype == "fasta_with_abspath":
        output = fasta_with_abspath
    elif args.outtype == "fasta_only_abspath":
        output = fasta_only_abspath
    elif args.outtype == "fastq_with_abspath":
        output = fastq_with_abspath
    elif args.outtype == "fastq_only_abspath":
        output = fastq_only_abspath
    elif args.outtype == "qual_with_abspath":
        output = qual_with_abspath
    elif args.outtype == "qual_only_abspath":
        output = qual_only_abspath
    #
    elif args.outtype == "fasta_with_filename":
        output = fasta_with_filename
    elif args.outtype == "fasta_only_filename":
        output = fasta_only_filename
    elif args.outtype == "fastq_with_filename":
        output = fastq_with_filename
    elif args.outtype == "fastq_only_filename":
        output = fastq_only_filename
    elif args.outtype == "qual_with_filename":
        output = qual_with_filename
    elif args.outtype == "qual_only_filename":
        output = qual_only_filename
    #
    elif args.outtype == "fasta_readstatsname":
        output = fasta_readstatsname
    elif args.outtype == "fasta_readstatsname_with_abspath":
        output = fasta_readstatsname_with_abspath
    elif args.outtype == "fasta_readstatsname_with_filename":
        output = fasta_readstatsname_with_filename
        
    elif args.outtype == "fastq_readstatsname":
        output = fastq_readstatsname
    elif args.outtype == "fastq_readstatsname_with_abspath":
        output = fastq_readstatsname_with_abspath
    elif args.outtype == "fastq_readstatsname_with_filename":
        output = fastq_readstatsname_with_filename
        
    elif args.outtype == "qual_readstatsname":
        output = qual_readstatsname
    elif args.outtype == "qual_readstatsname_with_abspath":
        output = qual_readstatsname_with_abspath
    elif args.outtype == "qual_readstatsname_with_filename":
        output = qual_readstatsname_with_filename
        
    ### get readtype fxn ###
    if args.readtype == "template":
        getread = get_template_read
    elif args.readtype == "complement":
        getread = get_complement_read
    elif args.readtype == "2d":
        getread = get_2d_read
    elif args.readtype == "molecule":
        getread = get_molecule_read
    elif args.readtype == "all":
        getread = get_all_reads
    elif args.readtype == "MoleQual":
        getread = get_molequal_read
    return output, getread




                
            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    output, getread = get_fast5tofastx_fxns(args)

    samflag=""
    if args.samflag:
        samflag = "F5:i:"

    falcon_i = 0
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=args.tarlite):
        if f5.is_not_corrupt() and f5.is_nonempty:
            ## counter in case using falcon options
            falcon_i += 1
            ## Process args.comments
            read = getread(f5, args.minlen, args.maxlen, args.minq, args.maxq, output, comments=args.comments, falcon_i=falcon_i, samflag=samflag)
            if read:
                print read


