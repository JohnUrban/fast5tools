import matplotlib
## may need following line for remote jobs (e.g. submitting batch scripts)
## matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

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
