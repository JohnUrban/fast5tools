#!/usr/bin/env python2.7
import h5py, os, sys
from fast5tools.f5class import *
from fast5tools.f5ops import *
from fast5tools.helperops import *
from fast5tools.fileListClass import *
import argparse
from glob import glob
import re
import sys
import string
from cStringIO import StringIO


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s,
find regular expressions in the base-called sequences.

Example:
G4 motif

Search sequences in set of fast5 files (or in a FASTA/FASTQ file) for a regular expression.
    Output BED file has default columns:
    1. Name of sequence \n
    2. Start of the match \n
    3. End of the match
    4. Strand (+/- relative to sequence given, NOT to be confised with template/complement reads.)
    5. Optional Matched sequence (--reportseq/-s)

    These can be changed with --outformat/-O which allows you to report name,start,end,strand,seq in any order.

    If --counts is used, default columns are:
    1. name
    2. pos strand count
    3. neg strand count
    4. total count
    
    This script will write out all positive strand entries of a given sequence followed by all negative strand entries.
    If name,start,end are used as first 3 columns, sortBed from BEDtools (or unix sort) can sort the BED file based on coordinates if needed.

John Urban (2015, 2016, 2017, 2018)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')


parser.add_argument('-n', '--nfiles', type=int, default=1000000000000,
                    help = '''This defaults to 1000000000000 in order to use all files (will you ever need to look at more than that?).
However, you can downsample with this option by adjusting this number to get info from the first N files.
Use --random to select N at random from the list.
Aim this script at a specific file for that file's contents.''')

parser.add_argument('-R', '--random', action='store_true', default=False,
                    help = '''Randomize what files are looked at.''')

parser.add_argument('-S', '--randomseed', type=int, default=False,
                    help = '''Randomize what files are looked at, but use given seed for reproducibility.''')

parser.add_argument('--filesused', type=str, default='qual_v_pos', help='''
''')

parser.add_argument('-o', '--outdir', type=str, default="./",
                    help = '''....''')

parser.add_argument('--prefix', type=str, default=None, help='''
Depending on options, a few files can be made.
1. BED coordinates (adjustable)
2. Plots...?
...TODO...''')

parser.add_argument('--filename', type=str, default=None, help='''
For output. Default None (stdout).
If a filename is given, filesused will be reported in a similarly named file ending with .filesused.fofn''')


parser.add_argument('-r', '--readtype', default='template',
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d', 'molecule', 'all', 'MoleQual'.
Default: template.
Molecule returns single fasta for each fast5 by following rules:
    if 2d present, return 2d.
    elif complement present with no 2d, return longer of template or complement.
    elif only template present, return template.
'MoleQual' is similar to molecule.
    It differs only in choosing between template and complement when a 2D is not present.
    Instead of choosing the longer one, it chooses the one with a higher quality mean quality score.''')



parser_not_fast5 = parser.add_mutually_exclusive_group()
parser_not_fast5.add_argument('--fasta', '-fa', action='store_true', default=False,
                           help='''Looking at a FASTA file or list of FASTA files, not FAST5s''')
parser_not_fast5.add_argument('--fastq', '-fq', action='store_true', default=False,
                           help='''Looking at a FASTQ file or list of FASTQ files, not FAST5s''')

parser_regex = parser.add_mutually_exclusive_group(required=True)

parser_regex.add_argument('--regex', 
               type= str, default=None, 
               help='''Required: Regex to be searched in the fasta input.
Matches to this regex will have + strand. This string passed to python
re.compile().''')

parser.add_argument('--regexrev',
               type= str, default=None, required=False,
               help='''The second regex to be searched in fasta input.
Matches to this regex will have - strand.
By default (None), --regexrev will be --regex complemented by replacing
'actguACTGU' with 'tgacaTGACA'.                                    
               ''')

parser.add_argument('--noreverse',
               action= 'store_true',
               help='''Do not search for any complement regular expression in the given sequences.
In each sequence, search only for regex given on the given strand.
Note: this does NOT mean to search only template reads for regex and it does NOT mean complement reads are ignored.
It means for all reads, only pay attention to the read sequence, not the inferred reverse complement of that sequence.
               ''')

parser.add_argument('--reportseq', '-s',
                action= 'store_true', default=False,
                help='''Report sequence of reg exp match in output.                                   
               ''')

parser.add_argument('--outformat', '-O',
                type=str, default='name,start,end,strand',
                help='''Provide comma-separated list of desired output infomation.
Options are name (sequence name), start (start of match), end (end of match),
strand (strand of match +/-), seq (sequence of match).
Default = 'name,start,end,strand'. --reportSeq/-s option changes default to: 'name,start,end,strand,seq'
Any other combination can be provided.
When using --counts, defaults to name,pos,neg
               ''')

parser.add_argument('--counts', '-c',
                action= 'store_true', default=False,
                help='''Report count for number of matches in each sequence instead of individually reporting all occurences in the sequence. 
               ''')

parser_regex.add_argument('--G4', 
                action= 'store_true', default=False,
                help='''Use the G4 motif as the regular expression.
A G4 is typically defined as: ([gG]{3,}\w{1,7}){3,}[gG]{3,}''')

parser.add_argument('--minG', type= int, default=3,
               help='''Only applies if --G4 flag used.
minG is the minimum number of Gs in a G tract.
The default minG value is 3.
This is typically the shortest allowable G-tract, but 2 is used in some cases to increase sensitivity.
Requiring longer G-tracts has more specificity, but lower sensitivity.
               ''')

parser.add_argument('--maxN', type= int, default=7,
               help='''Only applies if --G4 flag used.
maxN is the maximum number of number of Ns in loops between G tracts.
The default maxN value is 7.
Recently people have also often used maxN=15 -- i.e. ([gG]{3,}\w{1,15}){3,}[gG]{3,}
In general, allowing longer loops have more sensitivity, but lower specificity.
Some literature suggests that the probability of forming a G4 decreases with length.
               ''')
    
parser.add_argument('--numtracts',
                action= 'store_true', default=False,
                help='''Only applies if --G4 flag used.
For each G4 location, also report number of poly-G tracts inside G4 motif (and poly-C tracts in G4 complement motif). 
               ''')

parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')

parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')

parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')

parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')


parser.add_argument('--notarlite', action='store_true', default=False, help=''' The default methof (called tarlite) extracts 1 file from a given tarchive at a time, processes, and deletes it.
This options says to turn tarlite off resulting in extracting entire tarchive before proceeding (and finally deleting).
It is possible that --notarlite is faster, but at the expense of exceeding file number limits or disk storage quotas.
Nonetheless, the difference in speed is a lot smaller than the difference in space needed.
For example, not using tarlite will require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
A corollary is that tarlite just needs to be allowed to form 1 (or a few) files compared to what could be thousands to millions.
''')

parser.add_argument('--tarlite', action='store_true', default=False, help='''This legacy option is outdated.
However, it is kept here to avoid breaking pipelines that make use of it.
The tarlite approach is now default. Specifying this will not change that default behavior.
It will just prevent pipelines from breaking.
However, not specifying this will still also result in the tarlite approach.
Use --notarlite to turn it off.''')

args = parser.parse_args()





## should these return them compiled?

if __name__ == "__main__":
    
    # Process Args
    args.outdir = process_outdir(args.outdir)
    outfile = args.outdir + args.filename if (args.filename is not None) else None
    
    if args.G4:
        args.regex = get_g4_regex(args.minG, args.maxN)
        args.regexrev = get_g4_revregex(args.minG, args.maxN)
        
    ## define regexrev as complement of regex if none given
    if args.regexrev is None: ## then use complement of regex
        args.regexrev = get_complement_regex(args.regex)

    ## compile fwd and rev reg expressions
    re_f = re.compile(args.regex)
    re_r = re.compile(args.regexrev)
    fwd_re = re.compile(args.regex)
    rev_re = re.compile(args.regexrev)

    ## process automated reporting/output format requests
    if args.reportseq:
        args.outformat = 'name,start,end,strand,seq'                                                      
    if args.counts:
        args.outformat = 'name,pos,neg'
    if args.G4:
        if args.numtracts:
            args.numtracts = int(args.minG)
            if not args.counts:
                args.outformat += ",gtracts"

    ## allow fasta and fastq files be piped in as stdin
    if args.fast5 == '-':
        fastx_fh = sys.stdin
        if len(args.fast5) > 1:
            sys.exit('\nFast5tools g4/regex: Only one input file at a time can be processed:\n--fasta: %s\n' %(args.fast5))
    elif not args.fast5:
        fastx_fh = open(args.fast5)

    ## determine function to use
    if args.counts:
        regex_function = countSequence
    else:
        regex_function = parseSequence

    ## Tracking files used 
    if args.filename is not None:
        filesused_h = '.'.join(args.filename.split('.')[:-1]) + '.filesused.fofn'
    filesused = ''
    
    ## Execute:
    ## TODO: Can just use my FastXFileList class... for fa/fq, but would need to change regex_parseFasta/q stuff
    if args.fasta:
        for fa in FileList(args.fast5, extension=('.fa','.fasta', '.fna'), keep_tar_footprint_small=(not args.notarlite), downsample=args.nfiles, random=args.random, randomseed=args.randomseed):
            filesused += os.path.abspath(fa) + '\n'
            with open(fa) as fh:
                regex_parseFasta(fh, re_f, re_r, regex_function, outformat=args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)
    elif args.fastq:
        for fq in FileList(args.fast5, extension=('.fq','.fastq'), keep_tar_footprint_small=(not args.notarlite), downsample=args.nfiles, random=args.random, randomseed=args.randomseed):
            filesused += os.path.abspath(fq) + '\n'
            with open(fq) as fh:
                regex_parseFastq(fh, re_f, re_r, regex_function, outformat=args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)
    else:
        #regex_parseFast5(args.fast5, args.type, re_f, re_r, regex_function, outformat=args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)
        ## Iterate over fast5s
        for f5 in Fast5List(args.fast5, keep_tar_footprint_small=(not args.notarlite), filemode='r', downsample=args.nfiles, random=args.random, randomseed=args.randomseed):
            if meets_all_criteria(f5, args.readtype, args.minlen, args.maxlen, args.minq, args.maxq):
                filesused += f5.abspath + '\n'
                readtype = define_read_type(f5, args.readtype)
                seq = f5.get_seq(readtype)
                #seq_name = f5.abspath
                seq_name = f5.get_pore_info_name(readtype)
                regex_function(seq, seq_name, fwd_re, rev_re, args.outformat, count_gtract=args.numtracts, noreverse=args.noreverse)

    ## Files used
    process_filesused(trigger=args.filename, filesused=filesused, outdir=args.outdir)


