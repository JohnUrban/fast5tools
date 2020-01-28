import numpy as np
from collections import defaultdict
from fast5tools.helperops import *
try:
    Rless=False
    has_R=True
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    ##import rpy2.robjects.lib.ggplot2 as ggplot2
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    from fast5tools.edgeRops import *
except:
    Rless=True
    has_R=False
    sys.stderr.write('''
        Analyses that use R cannot be performed.
            Make sure to install:
                - Download R at: https://www.r-project.org/
                - the ggplot2 library in R: install.packagaes('ggplot2')
                    - ggplot2 fast5tools plots have been discontinued, so you can safely ignore this.
                - Install edgeR:
                    - See: https://bioconductor.org/packages/release/bioc/html/edgeR.html
                        source("https://bioconductor.org/biocLite.R")
                        biocLite("edgeR")
                - the rpy2 package from terminal/cmdline: pip install rpy2
                    - Install all necessary R stuff before rpy2, which may complain about R/library versions.
            Then try again.\n\n''')





## FUNCTIONS - some/many of these are no longer needed per se (the EdgeR class is the new approach)

if has_R:
    def load_edgeR():
        ## Load EdgeR Library
        robjects.r('''
            library(edgeR, quietly = TRUE)
        ''')

    def make_edgeR_function():
        ## Make EdgeR Function, "f"
        robjects.r('''
                f <- function(data, groups, genes, bcv){
                y <- DGEList(counts=data, group=groups, genes=as.matrix(genes, ncol=1))
                y <- calcNormFactors(y)
                y <- calcNormFactors(y)
                et <- exactTest(y, dispersion=bcv^2)
                return(et)
                }
                ''')

        ## Return EdgeR Function:
        return robjects.globalenv['f']

    def get_edgeR_function():
        load_edgeR()
        return make_edgeR_function()


    _get_DGE_list_ = get_edgeR_function()


    def get_DGE_list(data, groups, genes, bcv):
        return _get_DGE_list_(data, groups, genes, bcv)

## New approach
    def get_edgeR_functions():
        load_edgeR()
        return robjects.r['DGEList'], robjects.r['calcNormFactors'], robjects.r['exactTest']

    _DGEList, _calcNormFactors, _exactTest = get_edgeR_functions()

    def DGEList(counts, group, genes):
        '''counts is "data" in EdgeR class'''
        # This worked y = _DGEList(counts=np.array([[1,2],[1,3],[1,4],[1,5]]), group=np.array([1,2]), genes=np.array(['a','b','c','d']))
        return _DGEList(counts=counts, group=group, genes=genes)#

def calcNormFactors(dgelist):
    return _calcNormFactors(dgelist)

def exactTest(normdgelist, bcv=False, dispersion=False):
    ''' dispersion is just bcv**2'''
    if bcv is not False:
        return _exactTest(normdgelist, dispersion=bcv**2)
    elif dispersion is not False:
        return _exactTest(normdgelist, dispersion=dispersion)



def get_tag_table(dgelist, n_genes):
    ''' dgelist from get_DGE_list()
        n_genes is len(genes)'''
    return robjects.r.topTags(dgelist, n=n_genes)[0]

    

def get_counts(kmercounts, refcounts):
    '''Returns dict of dicts.
        ref = condition1.
        test = condition2.
    '''
    ## read in and equilibrate the 2 kmer count tables
    self.kmerdict, self.refdict = readInTwoKmerTables(kmercounts, refcounts)
    ## add total counts to respective kmerdicts
    kmerdict["Total"]  = totalKmerCount(kmerdict)
    refdict["Total"] = totalKmerCount(kmerdict2)
    ## add kmdericts to defaultdict(dict) that has keys condition1 and condition2
    counts = defaultdict(dict)
    counts['condition1'] = dict(refdict)
    counts['condition2'] = dict(kmerdict)
    return dict(counts)

def get_conditions(counts):
    ''' Simply returns condition names - which are keys in counts dict.'''
    return sorted(counts.keys())

def get_genes(counts, conditions=None):
    ''' Simply returns list of all unique kmers/genes found in the conditions to deal with 1 set.'''
    if conditions is None:
        conditions = get_conditions(counts)
    all_genes = []
    for condition in conditions:
        all_genes.extend(counts[condition].keys())
    return sorted(list(set(all_genes)))

def get_sizes(counts, conditions):
    ''' Simply returns total counts across all kmers/genes in each condition of counts dict.'''
    if conditions is None:
        conditions = get_conditions(counts)
    sizes = [work_counts[c]["Total"] for c in conditions]
    assert len(sizes) == 2
    return sizes

def get_conditions_and_genes(work_counts): 
    conditions = work_counts.keys()
    conditions.sort()
    all_genes = []
    for c in conditions:
        all_genes.extend(work_counts[c].keys())
    all_genes = list(set(all_genes))
    all_genes.sort()
    sizes = [work_counts[c]["Total"] for c in conditions]
    assert len(sizes) == 2
    all_genes.remove("Total")
    return conditions, all_genes, sizes
    
def edger_matrices(work_counts, conditions, all_genes, sizes):
    """Retrieve matrices for input into edgeR differential expression analysis.
    """
    groups = [1, 2]
    data = []
    final_genes = []
    for g in all_genes:
        cur_row = [int(work_counts[c][g]) for c in conditions]
        if sum(cur_row) > 0:
            data.append(cur_row)
            final_genes.append(g)
    return (numpy.array(data), numpy.array(groups), numpy.array(sizes),
            conditions, final_genes)


def run_kmer_diff(parser, args):
    counts = make_count_dict(parser, args)
    data, groups, sizes, conditions, genes = edger_matrices(counts)
    probs = run_edger(data, groups, sizes, genes, args)






#### CLASSES
if has_R:
    class EdgeR(object):
        def __init__(self, testdict, refdict, bcv=None):
            ''' testdict,refdict are both default_dictionaries_int with kmer/gene keys and counts as values.'''
            ## read in and equilibrate the 2 kmer count tables
            #self.testdict = testdict
            #self.refdict = refdict
            rpy2.robjects.numpy2ri.activate() ## This needs to be deactivated before the iter_rows() thing
            self.bcv = bcv
            self.conditions = ['condition1', 'condition2']
            self.groups = np.array([1,2])
            #self.test_size = sum(testdict.values())
            #self.ref_size = sum(refdict.values())
            self.sizes = np.array([sum(refdict.values()), sum(testdict.values())])
            self.counts = {'condition1':refdict, 'condition2':testdict}
            self.all_genes = sorted(list(set(testdict.keys() + refdict.keys())))
            self._define_data()
            self._ensureEqualKmerSets()
            self._debug_bcv=0.2
            #self.dgelist = get_DGE_list(self.data, self.groups, self.final_genes, self._debug_bcv)
            #self.tag_table = get_tag_table(self.dgelist, self.n_genes)
            self.dgelist = DGEList(self.data, self.groups, np.array(self.final_genes))
            self.dgelist = calcNormFactors(self.dgelist)
            self._determine_bcv(bcv)
            self.exactTest_results = exactTest(self.dgelist, dispersion=self.bcvsq)
            self.tag_table = get_tag_table(self.exactTest_results, self.n_genes)
            self.table_string = ''
            self.results = defaultdict(list)
            self._process_results()
            self.medianNorm = None

        def _define_data(self):
            self.data = []
            self.final_genes = []
            for gene in self.all_genes:
                curr_row =  [int(self.counts[condition][gene]) for condition in self.conditions]
                if sum(curr_row) > 0:
                    self.data.append(curr_row)
                    self.final_genes.append(gene)
            self.data = np.array(self.data)
            self.n_genes = len(self.final_genes)

        def _ensureEqualKmerSets(self):
            for key in self.final_genes:
                self.counts['condition1'][key]
                self.counts['condition2'][key]

        def _determine_bcv(self, bcv=None):
            if bcv is not None:
                self.bcvsq = bcv**2
            ## Else, get median sd from 
            else:
                #self.bcvsq = "auto"
                self.bcvsq = 0.2**2
                ## Could be set to a default value - e.g. 0.2
                ## Could be set to a diff value for each - e.g. a poisson-ish std given ref count, or std between conditions
                ## Could be set to other strings: "common", "trended", "tagwise"
                
        def _process_results(self):
            self.table_string = ''
            rpy2.robjects.numpy2ri.deactivate()
            for e in self.tag_table.iter_row():
                out = ("\t").join(str(e).split("\n")[1].split())
                ## Add to table_string
                self.table_string += out + '\n'

                ## Add elements to this class
                x = out.split("\t")
                self.results['k'].append(x[1])
                self.results['logfc'].append(float(x[2]))
                self.results['logcpm'].append(float(x[3]))
                self.results['p'].append(float(x[4]))
                self.results['fdr'].append(float(x[5]))
            rpy2.robjects.numpy2ri.activate()

        def get_countsdict(self):
            return self.counts

        def get_testdict(self):
            return self.counts['condition2']

        def get_refdict(self):
            return self.counts['condition1']

        def get_conditions(self):
            return self.conditions

        def get_groups(self):
            return self.goups

        def get_sizes(self):
            return self.get_sizes
        
        def get_data(self):
            return self.data

        def get_genes(self):
            assert self.final_genes == self.all_genes ## rm this line
            return self.final_genes

        def get_dge_list(self):
            return self.dgelist

        def get_norm_factors(self, condition=None):
            if condition is None:
                return self.dgelist[1][2]
            else:
                assert condition in [0,1]
                return self.dgelist[1][2][condition]

        def get_condition_sizes_from_dgelist(self):
            return self.dgelist[1][1]

        def get_dgelist_counts(self):
            return self.dgelist[0]

        def get_dge_list_genes(self):
            return self.dgelist[2][0]

        def get_normalized_counts(self):
            return self.get_dgelist_counts() * self.get_norm_factors()

        def get_average_normalized_counts_across_conditions(self):
            return self.get_normalized_counts().mean(1)
        
        def get_stdev_of_normalized_counts_across_conditions(self):
            return self.get_normalized_counts().std(1, ddof=1)
        
        def get_tag_table(self):
            return self.tag_table
        
        def get_table_string(self):
            return self.table_string

        def get_k_(self): #genes/names
            return list(self.exactTest_results[2][0])
        
        def get_k(self): #genes/names
            return self.results['k']
            #return list(self.exactTest_results[2][0])

        def get_logfc(self):
            return self.results['logfc']

        def get_logcpm(self):
            return self.results['logcpm']

        def get_pvalues(self):
            return self.results['p']

        def get_fdr(self):
            return self.results['fdr']
            
        def __str__(self):
            return self.get_table_string()
        





