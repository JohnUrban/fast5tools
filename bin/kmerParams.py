#!/usr/bin/env python2.7

import os, sys, itertools
from fast5tools.fileListClass import *
import argparse
from collections import defaultdict
import numpy as np


#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to event file(s) and/or directories of event files,

...summarize kmers seen.

For now, for k>1 it will take:
- average of means (either uniform or weighted by num data points)
- average of SDs (either uniform or weighted by num data points)
- sum of num data points

In future - maybe also allow just adding all component data points...
e.g. if dimer, then instead of the average of the 2 means, add both means to list for that dimer...

John Urban (2015, 2016, 2017, 2018)


    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('files', metavar='files', nargs='+',
                   type= str, 
                   help='''Paths to as many files and/or directories filled with files as you want.
Assumes all fast5 files have '.txt' extension.
If inside dir of dirs with .txt files, then can just do "*" to get all files from all dirs.
IF file extension is anything other than .txt, change it with --extension/-e''')

parser.add_argument('-e','--extension', type=str, default=".txt", help='''Provide the common file extension for target files.
Default is .txt.
Note that no other file type in target directory, fofn, etc should have this extension.
For K>6, it will just report information on kmers seen.''')



parser.add_argument('-k','--kmersize', type=int, default=1, help='''Provide kmer size. Default: 1. Common: 5 or 6.
It will report information on all possible ACGT kmers up to k=6, reporting None if no information found.
For K>6, it will just report information on kmers seen.''')

out_parser = parser.add_mutually_exclusive_group()
out_parser.add_argument('-w','--withdist', action='store_true', default=False,
                        help='''By default, the output is summary statistics for each kmer.
This says to also include a comma-separated list of the values that make up the distribution.
If --kmer_out_fprefix provided, then dist is split up into 3 files for mu, sd, n.
The stats are split up into 3 corresponding files as well.''')

out_parser.add_argument('-d','--onlydist', action='store_true', default=False,
                        help='''By default, the output is summary statistics for each kmer.
This says to instead produce a comma-separated list of the values that make up the distribution (without the summary stats).
If --kmer_out_fprefix provided, then dist is split up into 3 files for mu, sd, n.
Stats not reported.''')

parser.add_argument('-r','--rewrite', action='store_true', default=False, help='''By default, the output is
summary statistics for each kmer. This says to instead output the events re-written as kmers of length k.
Note that in rewritten events the start and length (i.e. n) columns do not follow the same logic as the input.
When using weighting,
interpret as the starting data point used to calculate mu and sd, and total number data points in weighting.
When using uniform,
interpret as the starting data point used to calculate mu and sd, and total number data points that were seen underlying bases used (though they were collapsed for uniform weighting).''')


##LOOK AT THIS....
parser.add_argument('-o', '--outdir', type=str, default='./',
                    help = '''If rewriting events or writing stats/dist to file w/ -O, store them in this directory. Default is working dir.''')


parser.add_argument('-O', '--kmer_out_fprefix', type=str, default='stdout',
                    help = '''By default, the kmer summary stats and/or distributions are collapsed to a single line and printed to stdout.
This will redirect into mu_stats/sd_stats/n_stats/mu_dist/sd_dist/n_dist files of given prefix.
These files can be collapsed together with 'paste' afterward if need be.''')


approach_parser = parser.add_mutually_exclusive_group(required=True)
approach_parser.add_argument('-W','--weighted', action='store_true', default=False, help=''' Mu and sd values for a kmer are calculated as a weighted mean by
weighting the constituent MUs and SDs assigned by tombo to inidividual bases according to the number of data points they each had.''')


approach_parser.add_argument('-U','--uniform', action='store_true', default=False, help='''Mu and sd values for a kmer are calculated as a uniformly-weighted mean by
averaging the values of each constituent MUs and SDs assigned by tombo to inidividual bases (w/o regard to how many datapoints each had).''')

approach_parser.add_argument('-P','--kmer_pos', type=int, default=False, help='''Mu and sd values for a kmer are taken from just one of
the values of the constituent MUs and SDs assigned by tombo to inidividual bases (instead of averaging). The underlying base that is used
is determined by the integer given w/ this option that denotes the 0-based index within the kmer that contains values to assign to the whole kmer.
This is in contrast to trying to average values of each base that makes up the kmer.''')

parser.add_argument('-G','--genomic_approach', action='store_true', default=False, help='''
Not yet implemented.

The 3 other options collect all values associated with a kmer from all reads.
Different coverage levels across the different places in the genome where that kmer occurs will give differential weights to the distribution of values for that kmer.
This is not necessarily bad, but might develop a model that is overfit to the particular dataset.
They are certainly fine solutions when there are few contexts (genomic sites) for a few or many of the kmers.
However, Tombo takes an approach where each genomic site that a kmer appears is treated like a "biological replicate" and all the values
assigned to it are treated like "technical replicates". Overall, it takes only the median value over each position.
The distribution of values for that kmer then is the size-N distribution of median values over the N genomic sites the kmer appears.
The genomic_position approach can, in theory, be applied to all 3 approaches above (weighted means, uniform means, and kmer_position).
In the future, this option will be a modifier for all 3 options.
For now, it is only set up to be a modifier for the --kmer_pos approach - i.e. mu and sd values will be taken from a given pos inside the kmer at each genomic site.

If using this approach, you need to provide genomic event text files as opposed to regular (tombo) event files.
''')


parser.add_argument('--notarlite', action='store_true', default=False, help=''' Default method extracts 1 file from a given tarchive at a time, processes, and deletes it.
This says to turn that off resulting in extracting entire tarchive before proceeding (and finally deleting).
It is possible that --notarlite is faster, but at the expense of exceeding file number limits or disk storage quotas.''')

parser.add_argument('-T', '--targzout', type=str, default=False,
                    help='''Only relevant for --rewrite.
                    This will make a file of the rewritten events, add the file to a growing tarchive,
                    followed by deleting the orignal file. This is a good option to prevent exceeding file number quotas, etc.
                    Provide a tarchive name - this script insists on using .tar.gz at the end, and will add it if absent.
                    This will put the tarchive into the specified outdir whether or not the outdir is included as part of the given name here.''')



args = parser.parse_args()

def tarpit(tarchive, fn, arcname=None):
    ''' Assumes tarhive is an already existing tarfile object for 'w|gz'
        fn is filename to add.
        Will add file to the arhive, then delete the original...
        Be careful.
        arcname is an alternate name to give the file inside the tarchive.'''
    if arcname is None:
        arcname = fn
    tarchive.add(name=fn, arcname=arcname)
    os.remove(fn)



def get_all_kmers(k):
    if k == 1:
        MERS = [''.join(e) for e in itertools.product("ACGT")]
    elif k == 2:
        MERS = [''.join(e) for e in itertools.product("ACGT","ACGT")]
    elif k == 3:
        MERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT")]
    elif k == 4:
        MERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT")]
    elif k == 5:
        MERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT")]
    elif k == 6:
        MERS = [''.join(e) for e in itertools.product("ACGT","ACGT","ACGT","ACGT","ACGT","ACGT")]
    else:
        MERS = []
        ## for all these, can do:
        ## MERS= [''.join(e) for e in itertools.product("ACGT", repeat=k)]
    return MERS

def initiateDict(MERS):
    kmers = defaultdict(dict)
    for mer in MERS:
        kmers[mer] = {'mu':[], 'sd':[], 'n':[]}
    return kmers

def initiateReadDict(MERS):
    kmers = defaultdict(int)
    for mer in MERS:
        kmers[mer] = 0
    return kmers

def prepareDict(kmers, mer):
    try:
        kmers[mer]['mu']
    except:
        kmers[mer] = {'mu':[], 'sd':[], 'n':[]}
    return kmers


def initiateGenomicDict(MERS):
    kmers = defaultdict(dict)
    for mer in MERS:
        kmers[mer] = {}
    return kmers

def prepareGenomicDict(kmers, mer, genomic_coords):
    #has kmer been added?
    try:
        kmers[mer]
    except:
        kmers[mer] = {}
    #have genomic_coords been added
    try:
        kmers[mer][genomic_coords]
    except:
        kmers[mer][genomic_coords] = {'mu':[], 'sd':[], 'n':[]}
    return kmers

def get_events(flines, nlines):
    mu = np.zeros(shape=nlines, dtype=float)
    sd = np.zeros(shape=nlines, dtype=float)
    start = np.zeros(shape=nlines, dtype=float)
    n = np.zeros(shape=nlines, dtype=float)
    b = np.zeros(shape=nlines, dtype=str)
    for i in range(nlines):
        line = flines[i].strip().split()
        mu[i] = float(line[0])
        sd[i] = float(line[1])
        start[i] = float(line[2])
        n[i] = float(line[3])
        b[i] = str(line[4])

    return mu, sd, start, n, b



def get_genomic_events(flines, nlines):
    mu = np.zeros(shape=nlines, dtype=float)
    sd = np.zeros(shape=nlines, dtype=float)
    start = np.zeros(shape=nlines, dtype=float)
    n = np.zeros(shape=nlines, dtype=float)
    b = np.zeros(shape=nlines, dtype=str)
    chrom = np.zeros(shape=nlines, dtype=str)
    pos = np.zeros(shape=nlines, dtype=float)
    strand = np.zeros(shape=nlines, dtype=str)
    for i in range(nlines):
        line = flines[i].strip().split()
        mu[i] = float(line[0])
        sd[i] = float(line[1])
        start[i] = float(line[2])
        n[i] = float(line[3])
        b[i] = str(line[4])
        chrom[i] = str(line[5]) ## can make this a single variable
        pos[i] = float(line[6])
        strand[i] = str(line[7]) ## can make this a single variable

    return mu, sd, start, n, b, chrom, pos, strand


def get_kmer_dist(mu, sd, start, n, b, k, kmers, nlines, uniform=False, rewrite_events=False, rewrite_name=False):
    '''
    mu = list of event means
    sd = list of event SDs
    start = list of event starts
    n = list of event lengths
    kmers = defaultdict that looks like {kmer:{mu:[], sd:[], n:[], ...}
    nlines = number events in txt file
    uniform = true/false on how to calculate mean and SD of kmer from individual tombo bases. Default False.
    rewrite_events = Create file that contains the collapsed kmer events.
    rewrite_name = name to give output file containing rewritten events.
    '''
    ## will return 0 or 1 for if kmer was detected
    kmerdetected = initiateReadDict(kmers.keys())
    ##
    if rewrite_events:
        rwh = open(rewrite_name, 'w')
    ## 1. Calculate first window
    ## 1a. Identify kmer
    mer = ('').join( [str(e) for e in b[0:k]] )
    kmerdetected[mer] = 1
    kmers = prepareDict(kmers, mer)
    ## 1b. Idenity total num data points in kmer
    nsum = n[0:k].sum()
    ## 1c. Identify kmer mean and sd - either uniform or weighted
    if uniform:
        musum = mu[0:k].sum() 
        sdsum = sd[0:k].sum()
        mumu = musum / float(k)
        sdmu = sdsum / float(k)
    else: #weighted
        musum = (mu[0:k] * n[0:k]).sum()
        sdsum = (sd[0:k] * n[0:k]).sum()
        mumu = musum / float(nsum)
        sdmu = sdsum / float(nsum)
    ## 1d. Add n/mu/sd values to growing lists
    kmers[mer]['n'].append( nsum )
    kmers[mer]['mu'].append( mumu )
    kmers[mer]['sd'].append( sdmu )
    ## 1e. Optionally add kmer event to re-written events file.
    if rewrite_events:
        re_event = ('\t').join([str(e) for e in [mumu, sdmu, start[0], nsum, mer]])
        rwh.write(re_event + '\n')
    
    ## 2. Calculate remaining windows
    for i in range(1, nlines-k+1, 1):
        ## 2a. Identify kmer
        mer = mer[1:] + b[i+k-1]
        kmerdetected[mer] = 1
        kmers = prepareDict(kmers, mer)
        ## 2b. Idenity total num data points in kmer
        nsum += n[i+k-1] - n[i-1]
        ## 2c. Identify kmer mean and sd - either uniform or weighted
        if uniform:
            musum += mu[i+k-1] - mu[i-1]
            sdsum += sd[i+k-1] - sd[i-1]
            mumu = musum / float(k)
            sdmu = sdsum / float(k)
        else: #weighted
            musum += (mu[i+k-1] * n[i+k-1]) - (mu[i-1] * n[i-1])
            sdsum += (sd[i+k-1] * n[i+k-1]) - (sd[i-1] * n[i-1])
            mumu = musum / float(nsum)
            sdmu = sdsum / float(nsum)
        ## 3d. Add n/mu/sd values to growing lists
        kmers[mer]['n'].append( nsum )
        kmers[mer]['mu'].append( mumu )
        kmers[mer]['sd'].append( sdmu )
        ## 2e. Optionally add kmer event to re-written events file.
        if rewrite_events:
            re_event = ('\t').join([str(e) for e in [mumu, sdmu, start[i], nsum, mer]])
            rwh.write(re_event + '\n')

        
    if rewrite_events:
        rwh.close()
        
    return kmers, kmerdetected



def get_kmer_pos_dist(mu, sd, start, n, b, k, kmers, nlines, kmer_pos=None, rewrite_events=False, rewrite_name=False):
    '''
    Different approach than get_kmer_dist().
    Does not take average of all underlying events.
    Only takes values from event at given position in kmer.
    This is more like tombo - although tombo only takes median value for all values assigned to each genomic pos.
    
    mu = list of event means
    sd = list of event SDs
    start = list of event starts
    n = list of event lengths
    kmers = defaultdict that looks like {kmer:{mu:[], sd:[], n:[], ...}
    nlines = number events in txt file
    kmer_pos is the 0-indexed pos w/in the kmer that is used to grab data (like in tombo) - 
    rewrite_events = Create file that contains the tombo events re-written w/ kmer instead of base.
                    --> Note that information is lost at this step on the read ends
    rewrite_name = name to give output file containing rewritten events.
    '''
    assert kmer_pos < k
    ## will return 0 or 1 for if kmer was detected
    kmerdetected = initiateReadDict(kmers.keys())
    ##
    if rewrite_events:
        rwh = open(rewrite_name, 'w')
    ##Calculate all windows
    for i in range(nlines-k+1):
        ## a. Identify kmer
        mer = ('').join( [str(e) for e in b[i:i+k]] )
        kmerdetected[mer] = 1
        kmers = prepareDict(kmers, mer)
        ## b. Append mean, sd, and num data points in kmer (given kmer pos) 
        kmers[mer]['n'].append( n[i+kmer_pos] )
        kmers[mer]['mu'].append( mu[i+kmer_pos] )
        kmers[mer]['sd'].append( sd[i+kmer_pos] )
        ## c. Optionally add kmer event to re-written events file.
        if rewrite_events:
            re_event = ('\t').join([str(e) for e in [mu[i+kmer_pos], sd[i+kmer_pos], start[i+kmer_pos], n[i+kmer_pos], mer]])
            rwh.write(re_event + '\n')
    if rewrite_events:
        rwh.close()
        
    return kmers, kmerdetected



def get_genomic_kmer_pos_dist(mu, sd, start, n, b, k, kmers, chrom, pos, strand, nlines, kmer_pos=None, rewrite_events=False, rewrite_name=False):
    '''
    
    '''
    assert kmer_pos < k
    ## will return 0 or 1 for if kmer was detected
    kmerdetected = initiateReadDict(kmers.keys())
    ##
    if rewrite_events:
        rwh = open(rewrite_name, 'w')
    ##Calculate all windows
    for i in range(nlines-k+1):
        ## a. Identify kmer
        mer = ('').join( [str(e) for e in b[i:i+k]] )
        kmerdetected[mer] = 1
        ## b. Identify genomic site
        genomic_coords = chrom[i] + '_' + str(pos[i+kmer_pos]) + '_' + strand ## genomic position of kmer is not nec first base of kmer, it is the "central base"
        kmers = prepareGenomicDict(kmers, mer, genomic_coords)
        
        ## c. Append mean, sd, and num data points in kmer (given kmer pos) 
        kmers[mer][genomic_coords]['n'].append( n[i+kmer_pos] )
        kmers[mer][genomic_coords]['mu'].append( mu[i+kmer_pos] )
        kmers[mer][genomic_coords]['sd'].append( sd[i+kmer_pos] )
        ##TODO____________________--------
        ## c. Optionally add kmer event to re-written events file.
        if rewrite_events:
            re_event = ('\t').join([str(e) for e in [mu[i+kmer_pos], sd[i+kmer_pos], start[i+kmer_pos], n[i+kmer_pos], mer, chrom[i+kmer_pos], pos[i+kmer_pos], strand[i+kmer_pos]]])
            rwh.write(re_event + '\n')
    if rewrite_events:
        rwh.close()

    
    return kmers, kmerdetected


def convert_genomic_kmers_to_kmers(gkmers):
    kmers = initiateDict()
    for kmer in gkmers.keys():
        kmers = prepareDict(kmers, kmer)
        for genomic_coord in gkmers[kmer].keys():
            for stat in ('mu', 'sd', 'n'):
                kmers[kmer][stat].append( np.median( gkmers[kmer][genomic_coord][stat] ) )
    return kmers
    

def get_stats(l):
    lnp = np.array(l, dtype=float)
    N = len(lnp)
    if N > 0:
        mu = np.mean(lnp)
        if N > 1:
            sd = np.std(lnp, ddof=1)
        else:
            sd = np.std(lnp, ddof=0)
        med = np.median(lnp)
        medsub = np.absolute(lnp - med)
        mad = np.median( medsub )
        mad2 = np.mean( medsub )
        scaled_sd = sd*0.67449
        scaled_mad = mad*1.4826
        scaled_mad2 = mad2*1.253314
        return mu, sd, med, mad, mad2, scaled_sd, scaled_mad, scaled_mad2, N
    else:
        return None, None, None, None, None, None, None, None, N



##def calc_weighted(mu, sd, start, n, b, k, kmers, nlines):
##    ## Calculate first window
##    musum = (mu[0:k] * n[0:k]).sum()
##    sdsum = (sd[0:k] * n[0:k]).sum()
##    nsum = n[0:k].sum()
##    mer = ('').join( [str(e) for e in b[0:k]] )
##    kmers = prepareDict(kmers, mer)
##    
##    kmers[mer]['mu'].append( musum / float(nsum) )
##    kmers[mer]['sd'].append( sdsum / float(nsum) )
##    kmers[mer]['n'].append( nsum )
##    
##    ## Calculate remaining windows
##    for i in range(1, nlines-k+1, 1):
##        musum += (mu[i+k-1] * n[i+k-1]) - (mu[i-1] * n[i-1])
##        sdsum += (sd[i+k-1] * n[i+k-1]) - (sd[i-1] * n[i-1])
##        nsum += n[i+k-1] - n[i-1]
##        mer = mer[1:] + b[i+k-1]
##
##        kmers[mer]['mu'].append( musum / float(nsum) )
##        kmers[mer]['sd'].append( sdsum / float(nsum) )
##        kmers[mer]['n'].append( nsum )
##
##    return kmers



##TODO:
## 1. Optionally write out distribution for each kmer
## 2. Optionally get mean/sd/med/mad of Mu, Sd, and N for each kmer
## 3. Optionally re-write events collapsed as kmers (instead of reporting dist or summary stats)

##def get_kmers_string(kmers):
##    string = ''
##    for mer in sorted(kmers.keys()):
##        if len(kmers[mer]['mu']) > 0:
##            string += ('\t').join([str(e) for e in [mer, kmers[mer]['mu']
            
#################################################
#### EXECUTE @@@@@@@@@@@@
#################################################

if __name__ == "__main__":
    tarlite=True
    if args.notarlite:
        tarlite=False

        
    if args.outdir:
        if not os.path.isdir(args.outdir):
            try:
                os.mkdir(args.outdir)
            except OSError:
                os.system("mkdir -p " + args.outdir)
        if args.outdir[-1] != "/":
            args.outdir += "/"
            
    if args.targzout and args.rewrite:
        if not args.targzout.endswith('.tar.gz'):
            args.targzout += '.tar.gz'
        tarpath = args.targzout
        if not tarpath.startswith(args.outdir):
            tarpath = args.outdir + tarpath
        tarchive = tarfile.open(tarpath, 'w|gz')

        
    if args.uniform:
        weighting = "uniform"
    elif args.weighted:
        weighting = "weighted"
    elif args.kmer_pos:
        weighting = "kmerpos-k" + str(args.kmersize) + '-pos' + str(args.kmer_pos)
    if args.genomic_approach:
        approach = 'genomicKmerMedians'
    else:
        approach = 'allKmersFound'
    
    k = args.kmersize
    kpos = args.kmer_pos
    MERS = get_all_kmers(k)
    if args.genomic_approach:
        kmers = initiateGenomicDict(MERS)
    else:
        kmers = initiateDict(MERS)
    nreadsdict = initiateReadDict(MERS)
    



    ## FOR EACH EVENTS TXT FILE, READ IN AND PROCESS KMER INFO
    filelist = FileList(args.files, extension=args.extension, keep_tar_footprint_small=tarlite)
    for fh in filelist:
        rewrite_name = args.outdir + ('.').join(os.path.basename(fh).split('.')[:-1]) + ".rewrite_as_" + weighting + '_' + str(k) + "mers.txt"
        ## READ IN AND ADD TO SUMMARY
        with open(fh, 'r') as f:
            flines = f.readlines()
        nlines = len(flines)
        if nlines > k:
            # get_kmer_pos step also does event re-writing to file inside function
            if args.genomic_approach:
                mu, sd, start, n, b, chrom, pos, strand = get_genomic_events(flines, nlines)
                if args.uniform or args.weighted:
                    pass
                elif args.kmer_pos:
                    kmers, kmerdetected = get_genomic_kmer_pos_dist(mu, sd, start, n, b, k, kmers, nlines, args.kmer_pos, args.rewrite, rewrite_name)
            else:
                mu, sd, start, n, b = get_events(flines, nlines)
                if args.uniform or args.weighted:
                    kmers, kmerdetected = get_kmer_dist(mu, sd, start, n, b, k, kmers, nlines, args.uniform, args.rewrite, rewrite_name)
                elif args.kmer_pos:
                    kmers, kmerdetected = get_kmer_pos_dist(mu, sd, start, n, b, k, kmers, nlines, args.kmer_pos, args.rewrite, rewrite_name)
            if args.targzout and args.rewrite:
                arcname = ('.').join(os.path.basename(fh).split('.')[:-1]) + ".rewrite_as_" + weighting + '_' + str(k) + "mers.txt"
                tarpit(tarchive, rewrite_name, arcname)
            for kmer, detected in kmerdetected.iteritems():
                if detected:
                    nreadsdict[kmer] += 1



    ## AFTER COLLECTING KMER MU/SD/N DISTS FROM ALL FILES, IF GENOMIC APPROACH, THEN TAKE MEDIAN VALUE FOR EACH GENOMIC SITE
    if args.genomic_approach:
        kmers = convert_genomic_kmers_to_kmers(kmers)
                    

    ## AFTER COLLECTING KMER MU/SD/N DISTS FROM ALL FILES, OPTIONALLY COMPUTE STATS
    if not args.onlydist:
        #calculate stats
        kmerstats = initiateDict(kmers.keys())
        #For mu, sd, and n values - calculate the mean, SD, median, mad, mad2, scaled_SD (by 0.67449), scaled_mad (by 1.4826), scaled_mad2 (by 1.253314), and num values used to calc stats for each
        #current out: mu, sd, med, mad, mad2, scaled_sd, scaled_mad, scaled_mad2, N
        #   where mad = median abs dev from median; and mad2 = mean abs dev from median
        # essentially - if using median and mad - one can use scaled_mad (or mad) if mad != 0 or scaled_mad2 (or mad2) otherwise.
        for kmer in kmers.keys():
            kmerstats[kmer]['mu'] = get_stats(kmers[kmer]['mu'])
            kmerstats[kmer]['sd'] = get_stats(kmers[kmer]['sd'])
            kmerstats[kmer]['n'] = get_stats(kmers[kmer]['n'])



            
    ## AFTER GETTING DIST VALS FROM FILES AND OPT COMPUTING STATS, OUTPUT DESIRED INFO
    if args.kmer_out_fprefix == 'stdout':
        ofh = sys.stdout
    if not args.onlydist and not (args.kmer_out_fprefix == 'stdout'):
        mu_ofh = open(args.outdir + args.kmer_out_fprefix + '_' + weighting + '_mu_stats_for_' + str(k) + 'mers' + '_from_' + approach + '.out', 'w')
        sd_ofh = open(args.outdir + args.kmer_out_fprefix + '_' + weighting + '_sd_stats_for_' + str(k) + 'mers' + '_from_' + approach + '.out', 'w')
        n_ofh = open(args.outdir + args.kmer_out_fprefix + '_n_stats_for_' + str(k) + 'mers.out', 'w')
    if (args.withdist or args.onlydist) and not (args.kmer_out_fprefix == 'stdout'):
        mu_dist_ofh = open(args.outdir + args.kmer_out_fprefix + '_' + weighting + '_mu_dist_for_' + str(k) + 'mers' + '_from_' + approach + '.out', 'w')
        sd_dist_ofh = open(args.outdir + args.kmer_out_fprefix + '_' + weighting + '_sd_dist_for_' + str(k) + 'mers' + '_from_' + approach + '.out', 'w')
        n_dist_ofh = open(args.outdir + args.kmer_out_fprefix + '_n_dist_for_' + str(k) + 'mers.out', 'w')
        
    for kmer in sorted(kmers.keys()):
        if args.withdist or args.onlydist:
        # get comma-separated list of distribution values
            mu_distlist = (',').join([str(e) for e in kmers[kmer]['mu']])
            sd_distlist = (',').join([str(e) for e in kmers[kmer]['sd']])
            n_distlist = (',').join([str(e) for e in kmers[kmer]['n']])
        if args.onlydist:
            if args.kmer_out_fprefix == 'stdout':
                stdout = ('\t').join( [kmer, mu_distlist, sd_distlist, n_distlist] )
                ofh.write(stdout + '\n')
            else:
                mu_dist_out = ('\t').join( [kmer, mu_distlist] )
                sd_dist_out = ('\t').join( [kmer, sd_distlist] )
                n_dist_out = ('\t').join( [kmer, n_distlist] )
                mu_dist_ofh(mu_dist_out + '\n')
                sd_dist_ofh(sd_dist_out + '\n')
                n_dist_ofh(n_dist_out + '\n')
        else:
            #current out: (mu, sd, med, mad, mad2, scaled_sd, scaled_mad, scaled_mad2, N) * 3, nreads, [mu_distlist, sd_distlist, n_distlist]
            # * 3 --> for mu, sd, and n values -- N will/should be the same for each
            mu = [str(e) for e in kmerstats[kmer]['mu']]
            sd = [str(e) for e in kmerstats[kmer]['sd']]
            n = [str(e) for e in kmerstats[kmer]['n']]
            nreads = [str(nreadsdict[kmer])]
            if args.kmer_out_fprefix == 'stdout':
                stdout_l = [kmer] + mu + sd + n + nreads
                if args.withdist:
                    stdout_l += [mu_distlist, sd_distlist, n_distlist]
                out = ('\t').join(stdout_l)
                ofh.write(out + '\n')
            else:
                muout = ('\t').join([kmer] + mu)
                sdout = ('\t').join([kmer] + sd)
                nout = ('\t').join([kmer] + n)
                mu_ofh.write(muout + '\n')
                sd_ofh.write(sdout + '\n')
                n_ofh.write(nout + '\n')
                if args.withdist:
                    mu_dist_out = ('\t').join( [kmer, mu_distlist] )
                    sd_dist_out = ('\t').join( [kmer, sd_distlist] )
                    n_dist_out = ('\t').join( [kmer, n_distlist] )
                    mu_dist_ofh.write(mu_dist_out + '\n')
                    sd_dist_ofh.write(sd_dist_out + '\n')
                    n_dist_ofh.write(n_dist_out + '\n')

    ## Close Files
    if args.kmer_out_fprefix != 'stdout':
        if not args.onlydist:
            mu_ofh.close()
            sd_ofh.close()
            n_ofh.close()
        if args.withdist or args.onlydist:
            mu_dist_ofh.close()
            sd_dist_ofh.close()
            n_dist_ofh.close()
    if args.targzout and args.rewrite:
        tarchive.close()
        
        


