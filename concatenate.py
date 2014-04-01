#! /usr/local/bin/python
import sys
import os

path=sys.argv[1]
outfile = open(path+'all_contig_seqs.fa', 'w')
for fi in os.listdir(path):
    if fi.endswith("sam"):
        filein=open(path+'/'+fi.replace('sam','fa'),'r')
        for line in filein:
            outfile.write(line)
        filein.close()
outfile.close()