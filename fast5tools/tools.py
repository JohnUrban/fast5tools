
## for bardecoder
def get_emissions_profile_probs(s, k, ket):
    ''' s is a string composed of ACGT
        k is kmer size (e.g. 5 or 6)
        ket = kmer_emissions_table = a dict with the kmer emissions read in'''
    npos = len(s)-k+1
    emit_probs = np.zeros([2,npos])
    for i in range(npos):
        emit_probs[0,i] = ket[s[i:i+k]][0] #mu
        emit_probs[1,i] = ket[s[i:i+k]][1]  #sigma
    return emit_probs

###end

def complement(DNAstring):
    DNAstring = DNAstring.upper()
    compString = ''
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    for base in DNAstring:
        compString = compString + complement[base]
    return compString

def reverseComplement(DNAstring):
    return complement(DNAstring[-1::-1])

def reverse_seq(DNAstring):
    return DNAstring[-1::-1]

def max_and_index(x):
    i=x.argmax()
    m=x[i]
    return i,m


def update_table_dict(d,l,keys, length=None):
    '''d is dict to update, l is list of values for k, the keys
        assumes l and d are of same length
        assumes k and l are paired by shared index'''
    if not length:
        length = len(l)
    for i in range(length):
        d[keys[i]].append(l[i])
    return d

def read_table(fh, keys, types=None):
    '''fh is a file path to a tsv file. keys are column names.
        lengths of types and keys = number colimns in table
        both keys and types should appear in same order as columns'''
    length = len(keys)
    if not types:
        types = [str]*length
        print types
    data = open(fh).readlines()
    table = defaultdict(list)
    for i in range(len(data)):
        line = data[i].strip().split("\t")
        line = [types[j](line[j]) for j in range(len(line))]
        table = update_table_dict(table,line,keys, length)
    return table




## These assume columns that may be from older fast5s/workflows
## -- can make them more flexible by reading the header from given files
def read_model_file(model_file, variantcolumn=False):
    if variantcolumn:
        keys = ["kmer","variant","level_mean","level_stdv","sd_mean","sd_stdv","weight"]
        types = [str] + [float]*6
    else:
        keys = ["kmer","level_mean","level_stdv","sd_mean","sd_stdv","weight"]
        types = [str] + [float]*5
    return read_table(model_file, keys, types)


def read_events_file(events_file, input_events=False):
    ''' file may contain input, template, or complement events '''
    if input_events:
        keys = ["mean", "stddev",  "start", "length"]
        types = [float]*4
    else:
        keys = ["mean", "stddev",  "start", "length", "model_state", "model_level", "move", "p_model_state", "mp_state", "p_mp_state", "p_A", "p_C", "p_G", "p_T"]
        types = [float]*4 + [str] + [float]*3 + [str] + [float]*5
    return read_table(events_file, keys, types)

