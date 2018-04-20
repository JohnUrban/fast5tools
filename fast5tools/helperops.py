import os, sys, re
import string
from cStringIO import StringIO
from Bio import SeqIO
from collections import Counter, defaultdict
from math import log10, log

## Fast5Tools
from fast5tools.f5class import *
from fast5tools.f5ops import *
from fast5tools.fileListClass import *


## 2018-04-20
## Some of the contents derive from old code imported from my poreminion tools
## Needs to be tested, cleaned up, re-done

def gc(x):
    c=0
    for b in x:
        c += b.upper() in 'GC'
    return 100.0*c/len(x)

def complement(seq):
    c=""
    for b in seq:
            if b == "A": c += "T"
            elif b == 'a': c+= 't'
            elif b == "C": c += "G"
            elif b == "c": c += "g"
            elif b == "G": c += "C"
            elif b == "g": c += "c"
            elif b == "T": c += "A"
            elif b == "t": c += "a"
    return c    

def revcomp(seq):
    return complement(seq)[-1::-1]

def case_counter(seq):
    return sum(1 for b in seq if b.islower())


def rev_comp_seq_dict(seq_dict):
    '''given a dictionary with sequences as values,
        give back dictionary with same KEYS but reverse complemented VALUES.'''
    out = {}
    for key,value in seq_dict.iteritems():
        out[key] = revcomp(value)
    return out

def rc(seq):
    intab='ACGTacgtUuNn'
    outtab='TGCAtgcaAaNn'
    trantab = maketrans(intab, outtab)
    return seq.translate(trantab)[-1::-1]

##########################################################################
## PROCESS XYZ
##########################################################################
def process_outdir(outdir):
    if not outdir.endswith('/'):
        outdir += '/'
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
    return outdir

def process_filesused(trigger, filesused, outdir='./'):
    ## trigger is typically args.filename
    ## Regardless, it needs to be a prefix of sorts with a throw-away extension
    if trigger is not None:
        ## Parse name
        filesused_h_ = '.'.join(trigger.split('.')[:-1]) + '.filesused.fofn'
        if filesused_h_ == '.filesused.fofn':
            filesused_h_ = trigger + '.filesused.fofn'
        filesused_h = outdir + filesused_h_ if outdir.endswith('/') else outdir + '/' + filesused_h_
        ## Write out
        with open(filesused_h, 'w') as fofnout:
            fofnout.write(filesused)
    else:
        sys.stderr.write(filesused)



##########################################################################
## REGEX
##########################################################################
def get_g4_regex(minG=3, maxN=7):
    return '([gG]{' + str(minG) + ',}\w{1,' + str(maxN) + '}){3,}[gG]{' + str(minG) + ',}'

def get_g4_revregex(minG=3, maxN=7):
    return '([cC]{' + str(minG) + ',}\w{1,' + str(maxN) + '}){3,}[cC]{' + str(minG) + ',}'

def get_complement_regex(regex, transtab=string.maketrans('actguACTGU', 'tgacaTGACA')):
    return regex.translate(transtab)


def parseSequence(seq, seq_name, fwd_re, rev_re, outformat, count_gtract=False, noreverse=False):
    '''if count_gtract is not False, it should be an int > 1'''
    for m in re.finditer(fwd_re, seq):
        formatDict = {'name':seq_name, 'start':m.start(), 'end':m.end(), 'strand':'+', 'seq':m.group(0)}
        if count_gtract:
            formatDict['gtracts'] = get_regex_count(formatDict['seq'], re.compile("[gG]{"+str(count_gtract)+",}"))
        print '\t'.join(str(x) for x in [formatDict[e] for e in outformat.split(',')])
    if noreverse is False:
        for m in re.finditer(rev_re, seq):
            formatDict = {'name':seq_name, 'start':m.start(), 'end':m.end(), 'strand':'-', 'seq':m.group(0)}
            if count_gtract:
                formatDict['gtracts'] = get_regex_count(formatDict['seq'], re.compile("[cC]{"+str(count_gtract)+",}"))
            print '\t'.join(str(x) for x in [formatDict[e] for e in outformat.split(',')])

def countSequence(seq, seq_name, fwd_re, rev_re, outformat, count_gtract=False, noreverse=False):
    ##count_gtract is not useful for this function for now, but needed to add it to play nice with other functions
    poscount = 0
    negcount = 0
    for m in re.finditer(fwd_re, seq):
        poscount += 1
    if noreverse is False:
        for m in re.finditer(rev_re, seq):
            negcount += 1
    formatDict = {'name':seq_name, 'pos':poscount, 'neg':negcount}
    print '\t'.join(str(x) for x in [formatDict[e] for e in outformat.split(',')])
    

def regex_parseFasta(seq_fh, re_f, re_r, regex_function, outformat, count_gtract=False, noreverse=False):
    '''if count_gtract is not False, it should be an int > 1'''
    seq = []
    line = (seq_fh.readline()).strip()
    seq_name = re.sub('^>', '', line)
    line = (seq_fh.readline()).strip()
    while True:
        while line.startswith('>') is False:
            seq.append(line)
            line= (seq_fh.readline()).strip()
            if line == '':
                break
        seq = ''.join(seq)

        regex_function(seq=seq, seq_name=seq_name, fwd_re=re_f, rev_re=re_r, outformat=outformat, count_gtract=count_gtract, noreverse=noreverse)
        
        seq_name = re.sub('^>', '', line)
        seq= []
        line= (seq_fh.readline()).strip()
        if line == '':
            break


def regex_parseFastq(seq_fh, re_f, re_r, regex_function, outformat, count_gtract=False, noreverse=False):
    '''if count_gtract is not False, it should be an int > 1'''
    line = (seq_fh.readline()).strip()
    seq_name = re.sub('^@', '', line)
    seq = (seq_fh.readline()).strip()
    while True:
        regex_function(seq=seq, seq_name=seq_name, fwd_re=re_f, rev_re=re_r, outformat=outformat, count_gtract=count_gtract, noreverse=noreverse)
        seq_fh.readline() ## skip + line
        seq_fh.readline() ## skip base call quality line
        line = (seq_fh.readline()).strip()   
        seq_name = re.sub('^@', '', line)
        seq = (seq_fh.readline()).strip()
        if line == '':
            break


def regex_parseFast5(fast5dir, seqtype, re_f, re_r, regex_function, outformat, count_gtract=False, noreverse=False):
    for fast5 in Fast5FileSet(fast5dir):
        fas = fast5.get_fastas(seqtype)
        for fa in fas:
            if fa is None:			
                continue
            basename,filename = fa.name.strip().split()
            filename = filename.split(".")[0].split("/")[-1]
            if filename.endswith("strand"):
                filename = filename[:-6]
            basename = basename.split("_")[-1]
            name = filename + basename
            regex_function(seq=fa.seq, seq_name=name, fwd_re=re_f, rev_re=re_r, outformat=outformat, count_gtract=count_gtract, noreverse=noreverse)
        fast5.close()



 
def get_regex_count(seq, regex):
    ## regex is a re.compile(string) object
    return len(re.findall(regex, seq))

def get_regex_count_in_fast5_fasta(fast5fasta, regex):
    if fast5fasta is not None:
        return get_regex_count(fast5fasta.seq, regex)
    else:
        return "-"
    
def get_regex_counts_in_fast5(fast5, regex, regex2=None):
    ## regex and regex2 are re.compile(string) objects
    ## often regex2 is the complement of regex (not necessary though)
    twoD = fast5.fastas.get("twodirections")
    template = fast5.fastas.get("template")
    complement = fast5.fastas.get("complement")
    counts = []
    counts.append(get_regex_count_in_fast5_fasta(twoD, regex))
    counts.append(get_regex_count_in_fast5_fasta(template, regex))
    counts.append(get_regex_count_in_fast5_fasta(complement, regex))
    if regex2 is not None:
        counts.append(get_regex_count_in_fast5_fasta(twoD, regex2))
        counts.append(get_regex_count_in_fast5_fasta(template, regex2))
        counts.append(get_regex_count_in_fast5_fasta(complement, regex2))
    return counts
        

##########################################################################
## KMER
##########################################################################
def log2(x):
    return log(x)/float(log(2))

def logbase(x, base):
    return np.log(x)/float(np.log(base))

def basecounter(kmer, basedict):
    '''basedict is a defaultdict(int) -- provide empty or non-empty one for update'''
    for b in kmer.upper():
        basedict[b] += 1
    return basedict

def countbases(seq, uppper=True):
    return Counter(seq.upper()) if upper else Counter(seq)

def kmercount_in_string(string, kmerdict, k):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    string = string.upper()
    for i in range(len(string)-k+1):
        kmerdict[string[i:i+k]] += 1
    return kmerdict

def writekmer(kmerdict, outfile):
    ''' kmerdict is a dict/default dict object'''
    fhout = open(outfile, 'w') if outfile is not None else sys.stdout
    total = float(sum(kmerdict.values()))
    for kmer in sorted(kmerdict.keys()): ## this could be a big memory suck
        fhout.write(kmer + '\t' + str(kmerdict[kmer]) + '\t' + str(kmerdict[kmer]/total) + '\n')
    if outfile is not None:
        fhout.close()


def kmercount_in_fast5(f5, readtype, k=6, kmerdict=None, rev_comp=False):
    if kmerdict is None:
        kmerdict = defaultdict(int)
    kmerdict = kmercount_in_string(f5.get_seq(readtype), kmerdict, k)
    if rev_comp:
        kmerdict = kmercount_in_string(f5.get_rev_comp_seq(readtype), kmerdict, k)
    return kmerdict


def kmercount_in_fastx(fh, fastx='fasta', k=6, kmerdict=None, rev_comp=False):
    if kmerdict is None:
        kmerdict = defaultdict(int)
    for fa in SeqIO.parse(fh, fastx):
        if fa is not None:
            kmerdict = kmercount_in_string(str(fa.seq), kmerdict, k)
            if rev_comp:
                kmerdict = kmercount_in_string(str(fa.reverse_complement().seq), kmerdict, k)
    return kmerdict


def kmercount_in_table(kmerCountFile):
    ''' Input is a kmer count table.'''
    kmerdict = defaultdict(int)
    if not isinstance(kmerCountFile, file):
        if isinstance(kmerCountFile, dict):
            return kmerCountFile
        elif isinstance(kmerdict, str):
            kmerCountFile = open(kmerCountFile)        
    # Get kmer info
    for line in open(fh, 'r'):
        kmer, count, prop = line.rstrip().split('\t')
        kmerdict[kmer] = int(count)
    # Close file
    kmerCountFile.close()
    # Return counts
    return kmerdict

def ensureEqualKmerSets(kmerdict1, kmerdict2):
    ''' dicts are defaultdict(int)'''
    ## make sure they have same set of kmers, for any missing elements in set A or B, add it with count of 0
    for key in kmerdict1:
        kmerdict2[key]
    for key in kmerdict2:
        kmerdict1[key]
    return kmerdict1, kmerdict2


def kmerDictToPlotData(kmerdict):
    data = defaultdict(list)
    for k in sorted(kmerdict.keys()):
        data['kmers'].append(k)
        data['counts'].append(kmerdict[k])
    return data

def totalKmerCount(kmerdict):
    return sum(kmerdict.values())

def readInTwoKmerTables(kmercounts, refcounts):
    kmerdict = kmercounts if isinstance(kmercounts, dict) else kmercount_in_table(kmercounts)
    refdict = refcounts if isinstance(refcounts, dict) else kmercount_in_table(kmercounts)
    ## ensure equal kmer sets
    kmerdict, refdict = ensureEqualKmerSets(kmerdict, refdict)
    return kmerdict, refdict


def make_count_dict(parser, args):
    ## read in and equilibrate the 2 kmer count tables
    kmerdict1, kmerdict2 = readInTwoKmerTables(parser, args)
    ## get total counts
    total1 = totalKmerCount(kmerdict1)
    total2 = totalKmerCount(kmerdict2)
    ## add totals to respective kmerdicts
    kmerdict1["Total"] = total1
    kmerdict2["Total"] = total2
    ## add kmdericts to defaultdict(dict) that has keys condition1 and condition2
    counts = defaultdict(dict)
    counts['condition1'] = dict(kmerdict1)
    counts['condition2'] = dict(kmerdict2)
    return dict(counts)



def run_kmer_counting(initial_list, k, readtype, revcomp, nfiles, random, randomseed, notarlite, fasta, fastq, minlen, maxlen, minq, maxq):
    ## Tracking files used 
    filesused = ''
    
    #Initialize -- should this be initialized with all possible kmers (all the 0-counts could ruin diff representation analysis though)
    kmerdict = defaultdict(int)
    
    ## Execute:
    if fasta:
        for fa in FileList(initial_list, extension=('.fa','.fasta', '.fna'), keep_tar_footprint_small=(not notarlite), downsample=nfiles, random=random, randomseed=randomseed):
            filesused += os.path.abspath(fa) + '\n'
            with open(fa) as fh:
                kmerdict = kmercount_in_fastx(fh, fastx='fasta', k=k, kmerdict=kmerdict, rev_comp=revcomp)
    elif fastq:
        for fq in FileList(initial_list, extension=('.fq','.fastq'), keep_tar_footprint_small=(not notarlite), downsample=nfiles, random=random, randomseed=randomseed):
            filesused += os.path.abspath(fq) + '\n'
            with open(fq) as fh:
                kmerdict = kmercount_in_fastx(fh, fastx='fastq', k=k, kmerdict=kmerdict, rev_comp=revcomp)
    else:
        ## Iterate over fast5s
        for f5 in Fast5List(initial_list, keep_tar_footprint_small=(not notarlite), filemode='r', downsample=nfiles, random=random, randomseed=randomseed):
            if meets_all_criteria(f5, readtype, minlen, maxlen, minq, maxq):
                filesused += f5.abspath + '\n'
                readtype = define_read_type(f5, readtype)
                kmerdict = kmercount_in_fast5(f5, readtype, k=k, kmerdict=kmerdict, rev_comp=revcomp)
    return kmerdict, filesused


def median_normalize(counts):
    counts = np.array(counts)
    med = np.median( counts )
    diffs = counts - med
    mad = np.median( np.absolute( diffs ) )
    z = diffs / mad
    return z, med, mad
