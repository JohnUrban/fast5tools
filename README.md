# fast5tools 0.4
Tools for working with fast5 files (Oxford Nanopore DNA/RNA sequencing data).


Add the fast5tools dir to your PYTHONPATH.
Add fast5tools/bin to your PATH

Try:

$ git clone https://github.com/JohnUrban/fast5tools.git
$ cd fast5tools
$ FAST5TOOLS=${PWD}
$ FAST5TOOLSPATHS=${FAST5TOOLS}:${FAST5TOOLS}/bin:${FAST5TOOLS}/fast5tools
$ export PATH=${FAST5TOOLSPATHS}:$PATH
$ export PYTHONPATH=${FAST5TOOLSPATHS}:$PYTHONPATH
$ python setup.py test

This will test various fast5 files in rundata/ that span from the beginning of MAP to recently. It is actually running tests/testf5\_class.py on them. 

Scripts are available in the bin directory and may be run as command line
programs once fast5tools is installed.


All scripts can read:
- a single fast5
- many fast5s listed in a row
- a file-of-file-names (with .fofn extension) of fast5s
- a tarchive with fast5s in dirs and subdirs
- acombination of those


Fast5 information:
- fast5attributes.py

Fast5 stats:
- fast5stats.py - generate a table with any fast5 info you want. 
- fast5TableSummary.py -- in progress -- The intention was to be able to give summary stats for any table with any fast5 info in any order (so long as you can tell it the order).
- fast5standardSummary.py

Plotting from stats tables:
- fast5plot.py -- in progress. 

Plotting (and kmer counting/plotting) from fast5s directly:
- fast5ReadLengthPlots.py
- fast5_qual_histogram.py
- fast5_qual_v_pos_boxplot.py
- fast5plotEvents.py
- fast5plotRaw.py
- kmerCountsFromFastx.py


Extracting Fasta/Fastq:
- fast5Tofastx -- convert to fasta, fastq, quals, integer quals -- or just details about the seqs (similar to fast5stats)
- fast5DerivedFastxMoleculeStats.py

Extracting events, raw data, model information:
- fast5ToEvents.py 
- fast5ToModelinfo.py 
- fast5toRaw.py


De-barcoding:
- fast5_sw_alignment_viz.py
- fast5_sw_bardecoder.py


Other:
- fast5staypos.py - get position of basecall "stays" in read (e.g. as performed in analyses in biorxiv paper below)
- fast5_regex_parser.py - identify positions of a regex in reads (e.g. as performed in analyses in biorxiv paper below)





Please cite fast5tools and/or poreminion as:
-------------------------------------------
Urban, J. M., Bliss, J., Lawrence, C. E. & Gerbi, S. A. 

Sequencing ultra-long DNA molecules with the Oxford Nanopore MinION. 

bioRxiv (Cold Spring Harbor Labs Journals, 2015). doi:10.1101/019281 

http://biorxiv.org/content/early/2015/05/20/019281

This paper contains the first descriptions and uses of poreminion, which evolved into fast5tools.



