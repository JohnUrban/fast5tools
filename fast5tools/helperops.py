import os

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





def process_outdir(outdir):
    if not outdir.endswith('/'):
        outdir += '/'
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
