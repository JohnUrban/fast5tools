#!/usr/bin/env python2.7
import argparse, sys, os,re
from collections import Counter, defaultdict
from scipy.stats import binom_test
from math import log
import numpy as np
#from fast5tools.samclass import *


parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Input = Two output files from samRefPosAlnComposition.py (with headers).
    Need to be sorted the same way.
    First data lines need to be from the same chromosome.
    Any chromosome not represented in one file should be removed from the other.

    Output = Two updated output RefPosAlnComp files.
        The same number of data points from both files are returned for each position.
        This means when one file lacks data for a position, the other file's data will not be returned.
        This more broadly means that files with more data over a position are always downsampled to the smaller number of data points.
        Downsampling is accomplished by sampling without replacement to mimic the effect of only seeing the first K reads.
            ...where K is the lower depth seen in the other file.
        
        Optional output:
        - Right now, only for mismatch probabilities
            - A table file with combined scores -- e.g. -log10(p-vals), -log2(fold changes)
            - Where there is data in both, a -log10(binomial p-val) is also calculated given a pro 


    Does not return stranded analysis.
        - Can hack it by giving RefPosAlnComp output of just the stranded info (for + or -) w/ colnames of nonstranded info.
        
    John Urban (2015, 2016, 2017, 2018)
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-f1', '--file1', 
                   type= str, required=True,
                   help='''Path to samRefPosAlnComposition file 1.''')
parser.add_argument('-f2', '--file2', 
                   type= str, required=True,
                   help='''Path to samRefPosAlnComposition file 1.''')
parser.add_argument('-c', '--compare', default=False, type=str, help='''If want comparison file, use this flag and provide file name to store in: --compare filename.ext.
For now, this only gives log2FoldChange and -log10Pval for marg2_pX between files (using "greater" in binomial test).
File1 is considered the test file and File2 is considered the control file.''')

parser.add_argument('-p', '--smallest_p', 
                   type=float, required=False, default=1e-300,
                   help='''Smallest p-value allowed. Default: 1e-300. Provide float.''')

parser.add_argument('-no', '--no-header', dest='no_header', default=False, action='store_true', help='''No header in output.''')
args = parser.parse_args()

##########################################################
'''FUNCTIONS'''
##########################################################
def newname(fh):
    b=os.path.basename(fh)
    d=os.path.dirname(fh)
    if d and not d.endswith('/'):
        d+='/'
    return d + 'updated_' + b

def next_line(fh):
    return parse_line(f.next().strip().split(), idict, typedict)

def parse_line(line, hdict, typedict):
    ''' input line is typically f.next().strip().split()
        hdict has colname->index
        idict has index->colname
    '''
    newdict = {}
    for idx,colname in idict.iteritems():
        #colname = idict[idx]
        newdict[colname] = typedict[colname](line[idx])
    return newdict

def same_chrom(fline, gline):
    return fline['chr'] == gline['chr']

def which_chrom_same_as_curr(fline, gline, currchr):
    ## Assumption is that only one changed -- if both did, we wouldn't be in this mess.
    assert fline['chr'] == currchr or gline['chr'] == currchr
    if fline['chr'] == currchr:
        return 0
    else:
        return 1

def posA_equals_posB(fline, gline):
    return fline['pos'] == gline['pos']

def posA_lower_than_posB(fline, gline):
    return fline['pos'] < gline['pos']

def posA_higher_than_posB(fline, gline):
    return fline['pos'] > gline['pos']

def get_stats(line1,line2,smallest_p=1e-300):
    log2fc = log((line1['marg2_pX']+1e-9)/(line2['marg2_pX']+1e-9))/log(2)
    ## Test if number mismatches in F are signif > num mismatches expected given prob_mismatch from G
    pval = max(binom_test(x=line1['X'], n=line1['X']+line1['='], p=line2['marg2_pX'], alternative='greater'), smallest_p)
    log10pval = (log(pval) / log(10))
    log10pval = 0 if log10pval == 0 else -1*log10pval
    return [log2fc, log10pval]
    
def make_comp_line(fline, gline, newfline, newgline, smallest_p=1e-300):
    ## chr, pos, log2fc, log10pval, new_log2fc, new_log10pval
    return ('\t').join([str(e) for e in [fline['chr'], fline['pos']] + get_stats(fline,gline,smallest_p) + get_stats(newfline,newgline,smallest_p)]) + '\n'

def make_line_out(line, colnames):
    return '\t'.join([str(out) for out in [line[e] for e in colnames]]) + '\n'

def make_header_out(colnames):
    return '\t'.join([str(i+1)+'='+colnames[i] for i in range(len(colnames))]) + '\n'

def make_compheader_out(colnames=['chr', 'pos', 'log2fc', 'log10pval', 'new_log2fc', 'new_log10pval']):
    return '\t'.join([str(i+1)+'='+colnames[i] for i in range(len(colnames))]) + '\n'

def downsample_line1_given_line2(line1, line2, colnames):
    charlist = [e for e in ''.join([base*int(line1[base]) for base in 'ACGTND'])]
    downsample = np.random.choice(charlist, size=int(line2['depth']), replace=False)
    counts = Counter(downsample)
    newline1 = new_line_given_counts(line1, counts, colnames)
    return newline1

def make_lambda_dict(allcolnames):
    lamdict = {}
    lamdict['prob'] = lambda x,y: str(float(x)/y)
    lamdict['id'] = lambda x: str(x)
    lambdadict = {}
##    for e in ['chr', 'pos', 'ref', 'depth', '='
##    lambdadict['chr'] = lambdict

def catch_zerodivisionerror(numerator, denominator):
    try:
        return numerator/denominator
    except ZeroDivisionError:
        return 'NA'


def get_insertion_info(line, newline):
    
    try:
        return newline['depth']*line['rI'], line['rI'], line['bI']
    except:
        return '.', '.','.'

def get_start_and_end_info(line, newline):
    scalar = newline['depth'] / float(line['depth'])
    try:
        return scalar * line['^'], scalar * line['$']
    except:
        return '.', '.'
    
def new_line_given_counts(line1, counts, colnames):
    newline = {}
    for e in ['chr', 'pos', 'ref']:
        newline[e] = line1[e]
    newline['depth'] = sum(counts.values())
    newline['='] = counts[line1['ref'].upper()]
    for e in 'ACGTDN':
        newline[e] = counts[e]
    newline['X'] = newline['depth'] - newline['='] - newline['D'] - newline['N']
    newline['I'], newline['rI'], newline['bI'] = get_insertion_info(line1, newline)
    marg2_total = float(sum([newline[e] for e in 'ACGT']))
    marg1_total = float(marg2_total + newline['D'])
    total = float(marg1_total + newline['N'])
    for e in '=XDNACGT':
        newline['p'+e] = catch_zerodivisionerror(newline[e], total)
    for e in '=XDACGT':
        newline['marg1_p'+e] = catch_zerodivisionerror(newline[e], marg1_total)
    for e in '=XACGT':
        newline['marg2_p'+e] = catch_zerodivisionerror(newline[e], marg2_total)
    newline['^'], newline['$'] = get_start_and_end_info(line1, newline)
    return newline
    
    
    
##########################################################
'''INPUTS'''
##########################################################
f = open(args.file1)
g = open(args.file2)
fheader = f.next()
gheader = g.next()
assert fheader == gheader

##########################################################
'''OUTPUTS'''
##########################################################
fout = open(newname(args.file1), 'w')
gout = open(newname(args.file2), 'w')
if args.compare:
    compareout = open(args.compare, 'w')

##########################################################
'''PARSE HEADER'''
##########################################################
#print fheader
header = [e.split('=',1) for e in fheader.split()]
#print header
#for idx,val in header:
#    print idx,val
hdict = {val:(int(idx)-1) for idx,val in header}
idict = {(int(idx)-1):val for idx,val in header}
typedict = {'chr':str, 'pos':int, 'ref':str, 'depth':float, '=':float, 'X':float, 'D':float, 'I':float, 'A':float, 'C':float, 'G':float, 'T':float, 'N':float, 'p=':float, 'pX':float, 'pD':float, 'pN':float, 'marg1_p=':float, 'marg1_pX':float, 'marg1_pD':float, 'marg2_p=':float, 'marg2_pX':float, 'pA':float, 'pC':float, 'pG':float, 'pT':float, 'marg1_pA':float, 'marg1_pC':float, 'marg1_pG':float, 'marg1_pT':float, 'marg2_pA':float, 'marg2_pC':float, 'marg2_pG':float, 'marg2_pT':float, 'rI':float, 'bI':float, '^':float, '$':float}
#Expecteding:
#1=chr   2=pos   3=ref   4=depth 5==     6=X     7=D     8=I     9=A     10=C    11=G    12=T    13=N    14=p=   15=pX   16=pD   17=pN   18=marg1_p=     19=marg1_pX     20=marg1_pD     21=marg2_p=     22=marg2_pX     23=pA   24=pC   25=pG   26=pT   27=marg1_pA     28=marg1_pC     29=marg1_pG     30=marg1_pT     31=marg2_pA     32=marg2_pC     33=marg2_pG     34=marg2_pT     35=rI   36=bI   37=^    38=$
#But could differ
allcolnames=['chr', 'pos', 'ref', 'depth', '=', 'X', 'D', 'I', 'A', 'C', 'G', 'T', 'N', 'p=', 'pX', 'pD', 'pN', 'marg1_p=', 'marg1_pX', 'marg1_pD', 'marg2_p=', 'marg2_pX', 'pA', 'pC', 'pG', 'pT', 'marg1_pA', 'marg1_pC', 'marg1_pG', 'marg1_pT', 'marg2_pA', 'marg2_pC', 'marg2_pG', 'marg2_pT', 'rI', 'bI', '^', '$']
colnames = [e for e in allcolnames if e in hdict.keys()]

##########################################################
'''PARSE INPUTS/WRITE OUTPUTS'''
##########################################################
# WRITE HEADERS
headerout = make_header_out(colnames)
fout.write(headerout)
gout.write(headerout)
if args.compare:
    compareout.write(make_compheader_out())

## INITIALIZE - Ensure they start out on the same chromosome -- else would need to make assumption on sort order to figure out which file to burn lines from.
fline = parse_line(f.next().strip().split(), idict, typedict)
gline = parse_line(g.next().strip().split(), idict, typedict)
i=0
while True:
    try:
        i+=1
        if i == 1:
            assert same_chrom(fline,gline)
            currchr = fline['chr']
        if same_chrom(fline,gline):
            if posA_equals_posB(fline, gline):                                                                                                                                       
                ## downsample if necessary
                if fline['depth'] > gline['depth']:
                    #downsample f - create modified fline
                    outfline = downsample_line1_given_line2(fline, gline, colnames)
                    outgline = gline
                elif fline['depth'] < gline['depth']:
                    #downsample g - create modified gline
                    outfline = fline
                    outgline = downsample_line1_given_line2(gline, fline, colnames)
                else:
                    outfline = fline
                    outgline = gline
                
                ## write out
                fout.write(make_line_out(outfline, colnames))
                gout.write(make_line_out(outgline, colnames))
                if args.compare:
                    #print fline['chr'], fline['pos'], fline['ref'], fline['depth'], gline['chr'], gline['pos'], gline['ref'], gline['depth'], "SH", log2fc, log10pval
                    compareout.write(make_comp_line(fline, gline, outfline, outgline, args.smallest_p))
 
                
                ## nextlines for f and g
                fline = parse_line(f.next().strip().split(), idict, typedict)
                gline = parse_line(g.next().strip().split(), idict, typedict)
                
            elif posA_lower_than_posB(fline,gline):
                while posA_lower_than_posB(fline,gline): ## while lower, nextline for f
                    fline = parse_line(f.next().strip().split(), idict, typedict)
            elif posA_higher_than_posB(fline,gline):
                ## while higher, nextline for g
                while posA_higher_than_posB(fline,gline):
                    gline = parse_line(g.next().strip().split(), idict, typedict)
            else:
                print "ERROR IN LOGIC."
                quit()
        else: #chrom changed in one of the files
            whichidx = which_chrom_same_as_curr(fline, gline, currchr)
            if whichidx == 0:
                while fline['chr'] != currchr:
                    fline = parse_line(f.next().strip().split(), idict, typedict)
            elif whichidx == 1:
                while gline['chr'] != currchr:
                    gline = parse_line(g.next().strip().split(), idict, typedict)
        
        
    except StopIteration:
        break


# CLOSE FILES
f.close()
g.close()
if args.compare:
    compareout.close()





