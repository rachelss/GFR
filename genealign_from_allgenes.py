#! /usr/local/bin/python
import os
import subprocess
import sys
from StringIO import StringIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob

def read_resort(species,seqdict):
    print 'Reading data: '+species
    for seq_record in SeqIO.parse(f, "fasta"):
        gene=seq_record.id[:]
        if gene not in seqdict:
            seqdict[gene]=list()
        seq_record.id=species.split('.')[0]
        seqdict[gene].append(seq_record)    
    return seqdict

##########################################################
fafiles = glob.glob("*.fa")
fafiles.remove('ref_genes.fa')
seqdict=dict()
for f in fafiles:
    seqdict = read_resort(f,seqdict)        #gene:(SeqRecord(species,sequence))

for gene, records in seqdict.iteritems():
    SeqIO.write(records, gene+'.fa', "fasta")
