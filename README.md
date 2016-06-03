# fast5tools
Tools for working with fast5 files (nanopore data).

Try:

$ bash test.sh

This will test various fast5 files in rundata/ that span from the beginning of MAP to recently. It is actually running testf5class.py on them. 
The tests are not comrehensive for all operations in the f5class. Initially, I was adding a test for each operation I added - but stopped. TODO.
Nonetheless, it has worked well so far. The idea is to run testf5class.py on a sample fast5 file from a new batch to see if they will cause problems and where.
When it tells us where, we can quickly fix it.

I have tried to make fast5tools a collection of pyscripts rather than a single command with subcommands like poreminion and poretools.
However, I am not adverse to switching to the "fast5tools sub-command" way. Ultimately, I'd like to leave room to parallelize the tools where possible.
It does not have all the functionality of poreminion or squiggler yet. I will merge stuff in -- with your help! :)
It can be pretty great though. 


- fast5stats.py - generate a table with any fast5 info you want. 
- fast5TableSummary.py -- in progress -- The intention was to be able to give summary stats for any table with any fast5 info in any order (so long as you can tell it the order).
- fast5plot.py us -- in progress. There are things in poreminion, squiggler, and poretools as well for plotting.
- fast5ToEvents.py -- self-explanatory
- fast5ToModelinfo.py -- self-explanatory
- fast5Tofastx.py -- convert to fasta, fastq, quals, integer quals -- or just details about the seqs (similar to fast5stats.py)
- testf5class.py -- run on a sample fast5 file from a new batch to see if they will cause problems and where. If we know where the new fast5 causes problems, we can quickly fix it.

All scripts can read:
- a single fast5
- many fast5s listed in a row
- a file-of-file-names (with .fofn extension) of fast5s
- a tarchive with fast5s in dirs and subdirs
- acombination of those
