#!/usr/bin/env python

##import os, sys
import argparse
import pandas
from fast5tools.f5dfclass import *


## JOHN URBAN (2015,2016)
#################################################
## Argument Parser
#################################################
parser = argparse.ArgumentParser(description = """

Given path(s) to standard table file output by either fast5stats.py or .

Assumes following columns (these and only these):
1 = molecule_name
2 = molecule length
3 = molecule mean quality score
4 = has complement
5 = has 2d
6 = 2d seq len
7 = template seq len
8 = complement seq len
9 = 2d mean q score
10 = template mean q score
11 = complement mean q score



Getting the top N longest reads and other sorting tips:
------------------------------------------------------
Although, this program offers some functionality for looking at the top N reads or molecules,
you may find that what you want is easily obtained by the following...

For 10 longest molecules:
sort -k2,2nr standardf5table.txt | head

For 10 longest 2d reads:
sort -k6,6nr standardf5table.txt | head

For 10 longest 2 reads, only report length, Q, and name:
sort -k6,6nr standardf5table.txt | head | awk 'OFS="\t" {print $6, $9, $1}'

For 10 highest quality template reads, only report qualities:
sort -k10,10nr standardf5table.txt | head | cut -f 10

Etc.

John Urban (2016)

    """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument('-f', '--file', 
                   type= str, required=True,
                   help='''Path to table file.''')


parser.add_argument('-d', '--delimiter', type=str, default='\t',
                    help='''Delimiter in table file. Default/assumes tab-delimited.''')


parser.add_argument('-H', '--header', action='store_true', default=False,
                    help='''Include header.''')

parser.add_argument('-V', '--vertical', action='store_true', default=False,
                    help='''Instead of horizontal 1-row/many-column output, do vertical many-row/1-column output.''')

parser.add_argument('-R', '--readtypeonly', action='store_true', default=False,
                    help='''Instead of full standard summary, return only the first 15 pieces of information about read types amongst molecules.''')

parser.add_argument('-L', '--longest', type=int, default=False,
                    help='''Instead of full standard summary, return only the top N longest of each read_type.
This overrides all other options. Provide integer for N.
See description above for other ways to get this type of infomation using unix sort.
''')

parser.add_argument('-LT', '--longestreadtype', type=str, default=all,
                    help='''Only show longest N for comma-separated read types. Choose: 'temp', 'comp', '2d', 'all'.  Default: 'all'.''')

parser.add_argument('-LQ', '--longestquality', type=float, default=0,
                    help='''Only show longest N for reads with mean quality > Q.  Default: 0.''')

parser.add_argument('-LC', '--longestcolumns', type=str, default=None,
                    help='''Show these columns only, in given order. Provide comma-separated list.
Use 'all' as short cut for all choices in order shown below.
Choices: 'name','mol_len','mol_q', 'hascomp','has2d','2d_len','temp_len','comp_len','2d_q','temp_q','comp_q'.
Default: length, quality, name (specific to given read type).
E.g., if temp chosen, defaults to: temp_len,temp_q,name''')

parser.add_argument('-LNR', '--longestnorank', action='store_true', default=False,
                    help='''Do not add rank column. Default: False.''')

parser.add_argument('-LNT', '--longestnotype', action='store_true', default=False,
                    help='''Do not add read type column. Default: False.''')

args = parser.parse_args()



#################################################
## deal with some of the arguments
#################################################


       


df = Fast5DataFrame(pandas.read_table(args.file, names=['name','mol_len','mol_q', 'hascomp','has2d','2d_len','temp_len','comp_len','2d_q','temp_q','comp_q'], na_values="-"))

if args.longest:
    addrank=True
    addrt=True
    if args.longestnorank:
        addrank=False
    if args.longestnotype:
        addrt=False
    allcolumns = ['name','mol_len','mol_q', 'hascomp','has2d','2d_len','temp_len','comp_len','2d_q','temp_q','comp_q']
    if "all" in args.longestreadtype:
        read_types = ["temp", "comp", "2d"]
    else:
        read_types = [e for e in args.longestreadtype.strip().split(",") if e in ["temp", "comp", "2d"]]
    if args.longestcolumns is not None:
        if 'all' in args.longestcolumns:
            columns = allcolumns
        else:
            columns = [e for e in args.longestcolumns.strip().split(",") if e in allcolumns]
    else:
        columns = None
    for read_type in read_types:
        df.print_N_longest(read_type, args.longest, args.longestquality, rank=addrank, columns=columns, add_read_type=addrt)
else:
    if args.readtypeonly:
        header, summary = df.read_type_summary()
    else:
        header, summary = df.standard_summary()
    assert len(header) == len(summary)
    if args.vertical:
        if args.header:
            for i in range(len(header)):
                print ("\t").join([header[i], str(summary[i])])
        else:
            for i in range(len(header)):
                print str(summary[i])
    else:
        if args.header:
            print ("\t").join(header)
        print ("\t").join([str(e) for e in summary])




