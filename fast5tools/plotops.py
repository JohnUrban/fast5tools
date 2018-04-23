import sys
from collections import defaultdict

# Plotting
import matplotlib
## may need following line for remote jobs (e.g. submitting batch scripts)
## matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


## Fast5Tools
from fast5tools.helperops import *


## 2018-04-20
## Some of the contents derive from old code imported from my poreminion tools
## Needs to be tested, cleaned up, re-done

def qualhist(quals, filename=None, minrange=0, maxrange=20, step=1, density=False, cumulative=False, text_only=True):
    ## TODO: histtype : bar, barstacked, step, stepfilled;; orientation : horizontal, vertical;;
    ## TODO: rwidth, color, 
    if quals:
        bins = range(minrange, maxrange, step)
        ylab = 'Density' if density else 'Frequency'
        ylab += ' (cumulative)' if cumulative else ''
        n, outbins, patches = plt.hist(x=quals, bins=bins, density=density, cumulative=cumulative)
        plt.xlabel("Base Quality")
        plt.ylabel(ylab)
        plt.xticks(rotation=65, fontsize=8)
    else:
        sys.stderr.write("No reads that meet criteria...\n")
    if text_only:
        hist_as_txt = '\n'.join([str(k)+'\t'+str(v) for k,v in zip(bins, n)]).strip()
        if filename is not None:
            with open(filename, 'w') as txtout:
                txtout.write(hist_as_txt)
        else:
            sys.stdout.write(hist_as_txt)
        
    else:
        if filename is not None:
            try:
                plt.savefig(filename)
                plt.close()
            except:
                sys.stderr.write("Unrecognized extension for %s!\nTry .pdf or .jpg or .png \n" % (plot_file))
        else:
                plt.show()


def update_qualpos(quals, qualpos, bin_width=1000, zscores=False, robust=False):
    ''' returns dictionary with keys=positions, values=lists of qual scores for that position
    qualpos = defaultdict(list)
    zpos = defaultdict(list)
    bin_width = integer
    '''
    if zscores or robust:
        if robust:
            qmean = np.median(quals)
            qSD = np.median( abs( np.array(quals) - qmean ) )
        else:
            qmean = np.mean(quals)
            qSD = np.std(quals)
    ctr = 0
    for q in quals:
        ctr += 1
        if zscores or robust:
            qualpos[1+int(ctr//bin_width)].append((q-qmean)/qSD) ## Extra 1 is to put in 1-based pos
        else:
            qualpos[1+int(ctr//bin_width)].append(q)  ## This is collecting information for bins default 1 kb
    return qualpos


def qualposplot(qualpos, bin_width, zscores=False, robust=False, filename=None):
    if zscores:
        ylab = "Quality Z-score"
        plotout = "qualZscore-Vs-pos"
    elif robust:
        ylab = "Robust Quality Z-score"
        plotout = "robustQualZscore-Vs-pos"
    else:
        ylab = "Quality score"
        plotout = "qual-Vs-pos"
    if qualpos.keys():
        data = [qualpos[e] for e in sorted(qualpos.keys())]
        plt.boxplot(data)
        xdetail = " (" + str(bin_width) + " bp bins)"
        plt.xlabel("Bin number in read" + xdetail)
        plt.ylabel(ylab)
        plt.xticks(rotation=65, fontsize=8)
    else:
        sys.stderr.write("No reads that meet criteria: cannot construct quality-Zscore  vs. position scatter plot...\n")
    if filename is not None:
        try:
            plt.savefig(filename)
            plt.close()
        except:
            sys.stderr.write("Unrecognized extension for %s!\nTry .pdf or .jpg or .png \n" % (plot_file))

    else:
            plt.show()





################################################################
##### 
################################################################

def choose_kmer_plot(kmerdict=False, refdict=False, gg=False):
    if kmerdict and not refdict:
        if not gg:
            singleTableKmerPlot(kmerdict)
        else:
            singleTablePlot_gg(parser, args)
    elif kmerdict and refdict:
        twoTableKmerScatterPlot(kmerdict, refdict)


def general_barplot(x, height, width=1.0, edgecolor='k', align='center', saveas=False):
    # For some, use: align='edge'
    plt.bar(x=x, height=height, width=width, edgecolor=edgecolor)
    if saveas:
        plt.savefig(saveas)
    else:
        plt.show()
    plt.close()


    
def singleTableKmerHist(kmercounts, density=False, cumulative=False, saveas=False):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    kmerdict = kmercounts if isinstance(kmercounts, dict) else kmercount_in_table(kmercounts)
    numKmers = len(kmerdict)
    data = kmerDictToPlotData(kmerdict)
    n, outbins, patches = plt.hist(x=data['counts'], density=density, cumulative=cumulative)
    if saveas:
        plt.savefig(saveas)
    else:
        plt.show()
    plt.close()
    
def singleTableKmerPlot(kmercounts, saveas=False):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    kmerdict = kmercounts if isinstance(kmercounts, dict) else kmercount_in_table(kmercounts)
    numKmers = len(kmerdict)
    data = kmerDictToPlotData(kmerdict)
    general_barplot(x=range(1,numKmers+1), height=data['counts'], width=1.0, saveas=saveas)

    


def twoTableKmerScatterPlot(kmercounts, refcounts, saveas=False):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    ## read in and equilibrate the 2 kmer count tables
    kmerdict, refdict = readInTwoKmerTables(kmercounts, refcounts)
    
    ## make approp data structures
    test = kmerDictToPlotData(kmerdict)
    reference = kmerDictToPlotData(refdict)

    plt.scatter(x=reference['counts'], y=test['counts'], s=10, marker='.')
    for i in range(len(test['kmers'])):
        plt.annotate(test['kmers'][i], (reference['counts'][i],test['counts'][i]), size='xx-small')

    plt.xlabel('Reference')
    plt.ylabel('Test')

    if saveas:
        plt.savefig(saveas)
    else:
        plt.show()
    plt.close()




def general_scatter(x, y, words=False, saveas=False, xlab="", ylab="", s=5, cex='xx-small', marker='.'):
    plt.scatter(x=x, y=y, s=s, marker=marker)
    if words:
        for i in range(len(words)):
            plt.annotate(words[i], (x[i], y[i]), size=cex)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.autoscale()
    if saveas:
        plt.savefig(saveas)
    else:
        plt.show()
    plt.close()



def twoTableKmer_MA_Plot(medNormObj, base=2, saveas=False, s=5, cex='xx-small'):
    x = medNormObj.get_logavg(base)
    y = medNormObj.get_logfc(base)
    k = medNormObj.get_genes()
    xlab = 'Average Log' + str(base) + ' Counts'
    ylab = 'Log' + str(base) + ' Fold Change'
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)
    

def volcanoPlot(logfc, p, k, saveas=False, xlab="log2(Fold Change)", ylab="-log10(p-value)", s=5, cex='xx-small'):
    '''
    logfc is expected to be log fold-changes
    p is expected to be p-values
    '''
    x = [e for e in logfc]
    y = [-1*log10(e) for e in p]
    #y = -1*logbase(p,base=10)
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)


def smearPlot(logfc, logcpm, k, saveas=False, xlab="log2(CPM)", ylab="log2(Fold Change)", s=5, cex='xx-small'):
    '''
    logfc is expected to be log fold-changes
    logcpm is expected to be log counts per million (I think)
    smearplot in poreminion was x=logfc, y=logcpm.
    But logcpm is the log average over both groups.
    Thus logcpm vs logfc is essentially an MA plot.
    So the old smear plot was just an MA plot forcing you to turn your head.
    
    '''
    x = [e for e in logcpm]
    y = [e for e in logfc]
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)


def alphabeticalPlot(y, k, saveas=False, xlab="kmer", ylab="log2(FC)", s=5, cex='xx-small'):
    '''Assumes given x is in same order as k.
        Example of y = logfc'''
    x = range(len(k))
    k, y = zip(*sorted(zip(k,y)))
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)

def gcPlot(y, k, saveas=False, xlab="Percent GC", ylab="log2(FC)", s=5, cex='xx-small'):
    '''Assumes given x is in same order as k.
        Example of y = logfc'''
    x = [gcbases(e) for e in k]
    x, k, y = zip(*sorted(zip(x,k,y)))
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)

def complexityPlot(y, k, saveas=False, xlab="Percent GC", ylab="log2(FC)", s=5, cex='xx-small'):
    '''Assumes given x is in same order as k.
        Example of y = logfc'''
    x = [gcbases(e) for e in k]
    x, k, y = zip(*sorted(zip(x,k,y)))
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)

def compressionPlot(y, k, saveas=False, xlab="Compression Length", ylab="log2(FC)", s=5, cex='xx-small'):
    '''Assumes given x is in same order as k.
        Example of y = logfc'''
    compress_lens = kmer_compression(k=len(k[0]))
    x = [compress_lens[e] for e in k]
    x, k, y = zip(*sorted(zip(x,k,y)))
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)


def twoTableKmerScatterPlotEdgeR(edgeRobj,saveas=False, xlab="TMM Norm Reference Count", ylab="TMM Norm Test Count", s=5, cex='xx-small'):
    nc = edgeRobj.get_normalized_counts()
    x = list(nc[:,0])
    y = list(nc[:,1])
    k = list(edgeRobj.get_dge_list_genes())
    #print x[:10], len(x)
    #print y[:10], len(y)
    #print k, len(k)
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex)





##### REPLACING SOON
    
def median_norm_analysis(kmercounts, refcounts, scale_to_ref=False, log_it=False, base=2):
    '''This function will/should be replaced by the MedNorm class'''
    ## read in and equilibrate the 2 kmer count tables
    kmerdict, refdict = readInTwoKmerTables(kmercounts, refcounts)
    
    ## make approp data structures
    test = kmerDictToPlotData(kmerdict)
    reference = kmerDictToPlotData(refdict)

    ## Scale Test to Reference
    test_z, test_med, test_mad = median_normalize(test['counts'])
    ref_z, ref_med, ref_mad = median_normalize(reference['counts'])


    ## Normalize to both ref_spread and ref_median -- this is analogous to comparing z-scores, so makes more sense to me
    test_z_to_ref = (test_z*ref_mad) + ref_med
    ref_z_to_ref = (ref_z*ref_mad) + ref_med
    ## Normalize to just ref_median, retaining test_spread
    ## test_z = (test_z*test_mad) + ref_med

    ## Can only take log of positive numbers.
    ## For now this assumes that scaling to the reference returns only positive numbers
    ## But this needs to be revisited as that assumption can be easily violated
    log_test_z_to_ref = logbase(test_z, base)
    log_ref_z_to_ref = logbase(ref_z, base)

    ## Return results - 3 tuples
    return (test_z, test_med, test_mad), (ref_z, ref_med, ref_mad), (test_z_to_ref, ref_z_to_ref)


def twoTableKmer_MA_Plot_(kmercounts, refcounts, saveas=False, scale_to_ref=False, logplot=False, base=2, s=5, cex='xx-small'):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    ## read in and equilibrate the 2 kmer count tables
    kmerdict, refdict = readInTwoKmerTables(kmercounts, refcounts)
    
    ## make approp data structures
    test = kmerDictToPlotData(kmerdict)
    reference = kmerDictToPlotData(refdict)

    ## Scale Test to Reference
    test_z, test_med, test_mad = median_normalize(test['counts'])
    ref_z, ref_med, ref_mad = median_normalize(reference['counts'])

    ## Optional: Scale Test to Reference
    if scale_to_ref:
        ## Normalize to both ref_spread and ref_median -- this is analogous to comparing z-scores, so makes more sense to me
        test_z = (test_z*ref_mad) + ref_med
        ref_z = (ref_z*ref_mad) + ref_med
        ## Normalize to just ref_median, retaining test_spread
        ## test_z = (test_z*test_mad) + ref_med

    ## Optional: Log
    if logplot:
        test_z = logbase(test_z, base)
        ref_z = logbase(ref_z, base)

    ## Get avg counts (x-axis)
    avg_z = (test_z + ref_z) / 2.0

    ## Get difference (y-axis)
    diffs = test_z - ref_z



    ## Handle axis labels
    ylab = 'Difference: Test - Reference'
    xlab = 'Average Counts'
    if logplot:
        ylab = 'Log'+str(base) + ' ' + ylab
        xlab = 'Average Log' +str(base) + ' Counts'
    if scale_to_ref:
        ylab += '\n(scaled to reference)'


    ## Scatter
    x = avg_z
    y = diffs
    k = test['kmers']
    general_scatter(x, y, k, saveas, xlab, ylab, s, cex) ## formerly s=10, now s=5
