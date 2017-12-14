#!/usr/bin/env python2.7

import h5py, os, sys, argparse
import cStringIO as StringIO
from Bio import SeqIO
from fast5tools.f5class import *
from fast5tools.barcodeclass import *
from fast5tools.f5ops import *
from fast5tools.barcodeops import *
from fast5tools.helperops import *
from glob import glob


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to fast5 file(s) and/or directories of fast5s,
return info on barcodes.

Running time proportional to barcode lengths and specified search lengths along reads.
That needs to be balanced with the fact that longer barcodes will perform better.

This is complementary to fast5_sw_bardecoder.py

It will allow you to visualize alignments when you want to spot check what the bardecoding process is selecting.

John Urban (2015, 2016, 2017)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('fast5', metavar='fast5', nargs='+',
                   type= str, 
                   help='''Paths to as many fast5 files and/or directories filled with fast5 files as you want.
Assumes all fast5 files have '.fast5' extension.
If inside dir of dirs with .fast5 files, then can just do "*" to get all files from all dirs.''')

parser.add_argument('-r', '--readtype', default="template",
                   type= str, 
                   help='''Choose type of fasta to get.
Choices: 'template', 'complement', '2d'.
Default: template.''')


parser.add_argument('-o', '--outprefix', default="fast5_sw_bardecoded",
                   type= str, 
                   help='''Choose an outprefix for files generated. Default: fast5_sw_bardecoded''')

parser_barcode = parser.add_mutually_exclusive_group(required=True)
parser_barcode.add_argument('-b', '--barcodetable', type=str, default=False,
                   help='''Path to file with barcode names and sequences.
This can be any tab-separated table file.
Name and sequence are expected to be in column 3 and 5 by default.
Change columns using barcolumns option.
Lines starting with # (e.g. header) are ignored.''')
parser_barcode.add_argument('-B', '--barcodefasta', type=str,  default=False,
                   help='''Path to fasta file with barcode names and sequences.''')

parser.add_argument('--barcolumns', type=str, default='3,5', 
                   help='''If using a barcode table file, define which colums have the barcode name and sequence.
Provide comma-separated, 1-based pair of columns.
Default: 3,5''')


parser.add_argument('--barstart', type=int, default=0, help='''Start location within read sequence to begin searching for barcode.
Can also think about this option as how much of the read to clip off or ignore in the search.
Default: 0.
We have used 60 as well to ignore adapter sequences that remain.''')

parser.add_argument('--barlength', type=int, default=False, help='''Expected maximum length of barcodes.
By default, barlength is detected as longest barcode sequence.
Specifying this gives a hard/constant number to use.
It is best to ensure it is at least as long as the longest barcode.
This needs to be higher than barmin option, else barmin is returned.''')

parser.add_argument('--barflank', type=int, default=False, help='''Amount of extra sequence to add on to barlength.
By default, barlength is detected as longest barcode sequence.
By default, the amount of flank or extra sequence length to search is the given or detected barlength (i.e. search len = 2*barlen).
Specifying this gives a hard/constant number to use.''')

parser.add_argument('--barmin', type=int, default=100, help='''Absolute minimum length to use for barcode search. Default: 100''')

parser.add_argument('--match', type=int, default=4, help='''Match parameter. Should be >= 0. Default: 4.
After lots of messing around, using 4/-2/-1/-1 seemed to work well.
Original settings were 1,-1,-1,-1 somewhat modeled after BWA ont2d type.
However, in Jain et al paper indels cause ~10 pct error and subs cause ~5 pct - meaning when trying to align it here gaps should be penalized less than mismatches.
Also, increasing the match reward encouraged the alignments to cover larger parts of barcodes, allowing better discrimination.
The --full_query option of python swalign attempts to do this, but with some weird results - so I think this is better.''')

parser.add_argument('--mismatch', type=int, default=-2, help='''Mismatch penalty. Should be <= 0. Default: -2.''')
parser.add_argument('--gap_open', type=int, default=-1, help='''Gap open penalty. Should be <= 0.Default: -1.''')
parser.add_argument('--gap_ext', type=int, default=-1, help='''Gap extension penalty. Should be <= 0.Default: -1.''')
parser.add_argument('--gap_decay', type=int, default=0, help='''Gap extend decay parameter. Should be >= 0. Default: 0.''')

parser_alntype = parser.add_mutually_exclusive_group()
parser_alntype.add_argument('--global_aln', action='store_true', default=False, help='''Default is local smith-waterman alignment. Set this to use global alignment as implemented in swalign.
Experimental. Cannot use with --full_query.
This is not recommended, especially with large search spaces for barcodes -- as barcodes begin to look equally unlikely.''')
parser_alntype.add_argument('--full_query', action='store_true', default=False, help='''Default is local smith-waterman alignment. Set this to use full query alignment as implemented in swalign.
Experimental. Cannot use with --global_aln.
This is not recommended as it seems to give weird results.
Instead try setting the match/mismatch/gap parameters to encourage full barcode alignments - which I attempted to do aready.''')

parser.add_argument('--maxscore', action='store_true', default=False, help='''By default read name and max score probability returned. Add max score to output.''')
parser.add_argument('--meanscore', action='store_true', default=False, help='''By default read name and max score probability returned. Add mean score to output.''')
parser.add_argument('--allscores', action='store_true', default=False, help='''By default read name and max score probability returned. Add all scores to output.''')
parser.add_argument('--allprobs', action='store_true', default=False, help='''By default read name and max score probability returned. Add all probabilities to output.''')
parser.add_argument('--barcode_coords', action='store_true', default=False, help='''By default read name and max score probability returned. Add the start and end coordinates of barcode in alignment.''')
parser.add_argument('--read_coords', action='store_true', default=False, help='''By default read name and max score probability returned. Add the start and end coordinates of read in alignment.''')
parser.add_argument('--cigar', action='store_true', default=False, help='''By default read name and max score probability returned. Add cig string from swalign.''')

parser.add_argument('-s', '--sequence', action='store_true', default=False, help='''Add read sequence to output.''')
parser.add_argument('-q', '--quals', action='store_true', default=False, help='''Add qual string to output.''')


parser.add_argument('--all', action='store_true', default=False, help='''Is to only show best alignment. This shows all.''')


## When two different barcodes are expected in a row at 5' end
parser_barcode2 = parser.add_mutually_exclusive_group()
parser_barcode2.add_argument('--barcodetable2', type=str, default=False,
                   help='''When two different barcodes are expected in a row at 5' end.''')
parser_barcode2.add_argument('--barcodefasta2', type=str,  default=False,
                   help='''When two different barcodes are expected in a row at 5' end.''')
parser.add_argument('--barstart2', type=int, default=0, help='''Start location within read sequence to begin searching for the second barcode.
Default: 0.
Will use match, mismatch, gap, barcolumns, barlength, barflank, and barmin parameters from barcode1.''')


## Search both ends
parser.add_argument('--bothstrands', action='store_true', default=False, help='''Search both strands of read for barcodes.
This will actually reverse complement the end of reads (only the amount of search length needed), and look for barcodes in them.''')

parser.add_argument('--blocksize', type=int, default=100, help='''Length of lines to print out for alignments. Default=100.''')


parser.add_argument('-c', '--comments', type=str, default=False, help='''Add fast5 info to output.
Default: no comments/False.
Leave any desired string here (you may need to enclode in quotes is there are spaces).
Alternatively, specify one of the following options:
base_info
pore_info
read_stats
event_stats
read_event_stats

''')



parser.add_argument('--minlen', type=int, default=0, help='''Only report reads >= minlen. Default: 0 bp.''')

parser.add_argument('--maxlen', type=int, default=int(3e9), help='''Only report reads <= maxlen. Default: 3 billion bp.''')

parser.add_argument('--minq', type=float, default=0, help='''Only report reads with mean quality scores >= Q. Default: 0.''')

parser.add_argument('--maxq', type=float, default=int(10e3), help='''Only report reads with mean quality scores <= Q.
Default: 10000 (this is orders of magnitude higher than normal max which are always < 20)''')



parser.add_argument('--outfile', type=str, default=False, help='''Default is to print to stdout. This will redirect into a file.''')



parser.add_argument('--tarlite', action='store_true', default=False, help=''' This method extracts 1 file from a given tarchive at a time, processes, and deletes it.
The older still-default routine extracts the entirety of all given tarchives at once, then processes files.
The default method will therefore require >2*tarchive amount of disk space (i.e. the tar.gz and its extracted contents).
The tarlite method only requires the disk space already taken by the tarchive and enough for 1 additional file at a time.
Tarlite may become the default method after some testing if it performs at similar speeds.''')


## TODO:
## 1. Test option -- provide full length known sequences. When barcode is chosen, alignment score of entire read
##      against entire known sequence is then computed (as well as maybe coordinates of the alignment)...
##      Perhaps global options will be used...
## 2. Accept adapter sequences to trim off before barcode alignment OR to use with the barcode alignment.

args = parser.parse_args()


#################################################
## Require read type to be set correctly
#################################################
args.readtype = assert_readtype(args.readtype, legaloptions="tc2")

#################################################
## Obtain barcodes
#################################################



                
            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    ## print to stdout?
    if not args.outfile:
        OUT = sys.stdout
    else:
        OUT = open(args.outfile, 'w')
        
    ## For now, this is controlled to grab the fastq with only the filename
    output = get_fast5tofastx_output_fxn('fastq_only_filename')
    getread = get_fast5tofastx_readtype_fxn(args.readtype)

    #get barcodes and search parameters 
    barcodes = read_in_barcodes(barcodetable=args.barcodetable, barcodefasta=args.barcodefasta, barcolumns=args.barcolumns)
    search_len = get_barcode_search_length(barcodes, barlength=args.barlength, barflank=args.barflank, minbarlen=args.barmin)
    search_start = args.barstart
    search_end = search_start + search_len

    ## Two positive strand barcodes True?
    two_pos_barcodes = args.barcodetable2 or args.barcodefasta2
    if two_pos_barcodes:
        barcodes2 = read_in_barcodes(barcodetable=args.barcodetable2, barcodefasta=args.barcodefasta2, barcolumns=args.barcolumns)
        search_start2 = args.barstart2
        search_end2 = search_start2 + search_len
        
    #Looking on both strands? Then rev comp barcodes
    if args.bothstrands:
        revbarcodes = rev_comp_seq_dict(barcodes)
        rev_search_end = -search_start
            
        if two_pos_barcodes:
            revbarcodes2 = read_in_barcodes(barcodetable=args.barcodetable2, barcodefasta=args.barcodefasta2, barcolumns=args.barcolumns)
           
        
    #create alignment objects
    sw = get_alignment_object(args.match, args.mismatch, args.gap_open, args.gap_ext, args.gap_decay, verbose=False, globalalign=args.global_aln, full_query=args.full_query)


    # What to print from best alignment
    # always give name and probability
    get = [0,1]
    header = ['readname']
    headeradds = ['barcode', 'probability']
    if args.maxscore:
        get.append( 2 )
        headeradds.append('maxscore')
    if args.meanscore:
        get.append( 3 )
        headeradds.append('meanscore')
    if args.allscores:
        get.append( 4 )
        headeradds.append('allscores')
    if args.allprobs:
        get.append( 5 )
        headeradds.append('allprobs')
    if args.barcode_coords:
        get.append( 6 )
        get.append( 7 )
        headeradds.append('bc_start')
        headeradds.append('bc_end')
    if args.read_coords:
        get.append( 8 )
        get.append( 9 )
        headeradds.append('read_start')
        headeradds.append('read_end')
    if args.cigar:
        get.append( 10 )
        headeradds.append('cigar')
    for add in headeradds:
        header.append( 'set1_topstrand_' + add )
    if two_pos_barcodes:
        for add in headeradds:
            header.append( 'set2_topstrand_' + add )
    if args.bothstrands:
        for add in headeradds:
            header.append( 'set1_bottomstrand_' + add )
        if two_pos_barcodes:
            for add in headeradds:
                header.append( 'set2_bottomstrand_' + add )
    if args.comments:
        header.append( 'f5info' )
    if args.sequence:
        header.append('f5_sequence')
    if args.quals:
        header.append('f5_quals')
        
        
    ## Set these to empty
    seq=[]
    quals=[]
    aln2 = []
    revaln1 = []
    revaln2 = []


    # execute for loop
    for f5 in Fast5List(args.fast5, keep_tar_footprint_small=args.tarlite):
        if f5.is_not_corrupt() and f5.is_nonempty:
            read = getread(f5, args.minlen, args.maxlen, args.minq, args.maxq, output, comments=args.comments) #, samflag=samflag)
            if read:
                #Remove '@' from name and newline from end while breaking fastq up into list of 4 fields.
                read = read[1:].strip().split('\n')
                name = read[0].split()
                if len(name) == 1:
                    readname = [name[0]]
                    comments = []
                elif len(name) == 2:
                    readname = [name[0]]
                    comments = [name[1]]
                if args.sequence:
                    seq = [read[1]]
                if args.quals:
                    quals = [read[3]]

                # Look for first (maybe only) barcode on top strand
                bcaln = BarcodeChoice(barcodes, sw, read[1], search_start, search_end, ref_name=read[0], use_entire_barcode=True, compute_all=True)
                aln1 = [str(e) for e in bcaln.get_subset_maxbar_list(get)]
                out = ('\t').join( [">"+readname[0]] + aln1 + comments + seq + quals) + '\n\n'
                out += bcaln.get_all_formatted_pairwise_alignments(blocksize=args.blocksize) ## add in args to toggle minion error params
                OUT.write( out + '\n' )
##                quit()
                
                ##12/12/17: TODO: update the bottom 3 alignment sections as done for aln1


                
##                # If opted, look for second barcode on top strand
##                if two_pos_barcodes:
##                    ans2, all2 = choose_barcode(barcodes2, sw, read[1], search_start2, search_end2, ref_name=read[0], use_entire_barcode=False, return_alns=args.all)               
##                    aln2 = [str(ans2[i]) for i in get]
##
##                # Potentially look at the bottom strand by using reverse complement barcodes
##                if args.bothstrands:
##                    revcompseq = revcomp(read[1][-search_end:])
##                    revans, revall1 = choose_barcode(barcodes, sw, revcompseq, search_start, search_end, ref_name='revcomp_'+read[0], use_entire_barcode=False, return_alns=args.all)
##                    revaln1 = [str(revans[i]) for i in get]
##                    # If opted, look for second barcode on top strand
##                    if two_pos_barcodes:
##                        revans2, revall2 = choose_barcode(barcodes2, sw, revcompseq, search_start2, search_end2, ref_name='revcomp_'+read[0], use_entire_barcode=False, return_alns=args.all)               
##                        revaln2 = [str(revans2[i]) for i in get]
##                        
##                out = ('\t').join( ">" + readname + aln1 + aln2 + revaln1 + revaln2 + comments + seq + quals)
##                OUT.write( out + '\n' )
                
                
    ## closing out
    if args.outfile:
        OUT.close()




