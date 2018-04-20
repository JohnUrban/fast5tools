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
##### Non-EdgeR Kmer Plots
################################################################
## Need to test all of these.
def choose_kmer_plot(kmerdict=False, refdict=False, gg=False):
    if kmerdict and not refdict:
        if not gg:
            singleTableKmerPlot(kmerdict)
        else:
            singleTablePlot_gg(parser, args)
    elif kmerdict and refdict:
        twoTableKmerScatterPlot(kmerdict, refdict)



def singleTableKmerHist(kmercounts, density=False, cumulative=False):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    kmerdict = kmercounts if isinstance(kmercounts, dict) else kmercount_in_table(kmercounts)
    numKmers = len(kmerdict)
    data = kmerDictToPlotData(kmerdict)
    n, outbins, patches = plt.hist(x=data['counts'], density=density, cumulative=cumulative)
    plt.show()
    
def singleTableKmerPlot(kmercounts):
    ''' kmerdict is a defaultdict(int)
        It can take both empty and non-empty kmerdicts
        returns update of the input kmerdict given the input string and k'''
    kmerdict = kmercounts if isinstance(kmercounts, dict) else kmercount_in_table(kmercounts)
    numKmers = len(kmerdict)
    data = kmerDictToPlotData(kmerdict)
    plt.bar(x=range(1,numKmers+1), height=data['counts'], width=1.0)
    plt.show()
    


def twoTableKmerScatterPlot(kmercounts, refcounts, saveas=None):
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

    if saveas is not None and (saveas.endswith(".pdf") or saveas.endswith(".jpg")):
            plt.savefig(saveas)
    else:
        plt.show()

def twoTableKmer_MA_Plot(kmercounts, refcounts, saveas=None, scale_to_ref=False, logplot=False, base=2):
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

    ## Scatter
    plt.scatter(x=avg_z, y=diffs, s=10, marker='.')
    for i in range(len(test['kmers'])):
        plt.annotate(test['kmers'][i], (avg_z[i],diffs[i]), size='xx-small')

    ## Handle axis labels
    ylab = 'Difference: Test - Reference'
    xlab = 'Average Counts'
    if logplot:
        ylab = 'Log'+str(base) + ' ' + ylab
        xlab = 'Average Log' +str(base) + ' Counts'
    if scale_to_ref:
        ylab += '\n(scaled to reference)'

    plt.xlabel(xlab)
    plt.ylabel(ylab)

    ## Saving/Showing
    if saveas is not None and (saveas.endswith(".pdf") or saveas.endswith(".jpg")):
            plt.savefig(saveas)
    else:
        plt.show()


