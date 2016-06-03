import numpy as np
######################
#  Length stats
######################

def NX(l, x=[25,50,75], G=False):
        """
        Returns NX for all x for a list of numbers l.
        Default: N25, N50, N75
        Assumes all values in list x are between 0 and 100.
        Interpretation: When NX = NX_value, X% of data (in bp) is contained in reads at least NX_value bp long.
        """
        ## assumes both l and x are sorted
	if isinstance(l, list) and isinstance(x, list) and G:
            l = l[:]
            x = x[:]
            nxsum = 0
            nxvalues = {e:0 for e in x}
            for e in x:
                    xpct = G*e/100.0
                    while nxsum < xpct and l:
                            nxsum += l[-1]
                            lastsize = l.pop()
                    nxvalues[e] = lastsize
            return nxvalues

	else:
            return None

def e_size(l,G=False):
    #G is usually sum(l)
    return sum([e**2 for e in l])/float(G)

def basic_stats(l, msg):
    pass

def size_stats(l,x, msg=""):
    '''l and x are lists'''
    print msg
    # sort sizes
    l.sort()

    ## get X values for NX stats and sort
    x.sort()

    ## Get N reads
    N = len(l)
    print "N = %d" % (N)

    ## Get sum of data
    A = sum(l)
    print "Total length = %d" % (A)

    ## Get max length:
    MAX = max(l)
    print "Max length = %d" % (MAX)

    ## Get min length:
    MIN = min(l)
    print "Min length = %d" % (MIN)

    ## Get mean length
    MEAN = np.mean(l)
    print "Mean length = %d" % (MEAN)

    ## Get median contig size
    MEDIAN = np.median(l)
    print "Median length = %d" % (MEDIAN)

    ## Get NX values
    nxvalues = NX(l,x,G=A)
    for e in x:
        print "N%s length = %d" % (str(e), nxvalues[e])

    ## expected read length
    E = e_size(l,G=A)
    print "Expected size = %d" % (E)

    ##number reads >= X 
    print "N reads >= 10kb:"
    print sum(np.greater_equal(l,50e3))
    print "N reads >= 25kb:"
    print sum(np.greater_equal(l,50e3))
    print "N reads >= 50kb:"
    print sum(np.greater_equal(l,50e3))
    print "N reads >= 75kb:"
    print sum(np.greater_equal(l,75e3))
    print "N reads >= 100kb:"
    print sum(np.greater_equal(l,100e3))
    ##TODO: also want size data from reads >=X
    ## ALSO in diff fxn if both length and Q available - do longest with Q>x, etc
    print



def Q_stats(l,msg):
    print msg

    ## Get N 
    N = len(l)
    print "N = %d" % (N)

    ## Get max 
    MAX = max(l)
    print "Max Q = %f" % (MAX)

    ## Get min 
    MIN = min(l)
    print "Min Q = %f" % (MIN)

    ## Get mean 
    MEAN = np.mean(l)
    print "Mean Q = %f" % (MEAN)

    ## Get median 
    MEDIAN = np.median(l)
    print "Median Q = %f" % (MEDIAN)
    print

def read_type_stats(l,readtype="",msg=""):
    print msg
    ## Get N reads
    N = len(l)
    print "Total N = %d" % (N)
    n = sum(l)
    print "n with %s = %s" % (readtype, n)
    pct = 100.0*n/N
    print "Percent with %s = %s" % (readtype, pct)
    print
    


f5summary = {}
f5summary[1] = lambda l: "There were %d basenames found.\n" % (len(l))
f5summary[2] = lambda l,x: size_stats(l,x,msg="Molecule length stats:")

f5summary[3] = lambda l: read_type_stats(l,"Complement",msg="Molecules with complement reads:")
f5summary[4] = lambda l: read_type_stats(l,"2D",msg="Molecules with 2D reads:")

f5summary[5] = lambda l,x: size_stats(l,x,msg="2D length stats:")
f5summary[6] = lambda l,x: size_stats(l,x,msg="Template length stats:")
f5summary[7] = lambda l,x: size_stats(l,x,msg="Complement length stats:")
f5summary[8] = lambda l: Q_stats(l,msg="2D Q-score stats:")
f5summary[9] = lambda l: Q_stats(l,msg="Template Q-score stats:")
f5summary[10] = lambda l: Q_stats(l,msg="Complement Q-score stats:")

## pct of comp with 2D
## pct of 2D in pass
## 

    
