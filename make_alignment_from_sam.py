#! /usr/local/bin/python
#take sam file aligned to single contig and produce whole alignment
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio import AlignIO,SeqIO
from collections import Counter

#######################
def sep_cigar(cigar):
    d=list()
    code=list()
    nums=list()
    for i,c in enumerate(cigar):
        if c.isdigit():
            d.append(c)
        else:
            code.append(c)
            d="".join(d)
            nums.append(int(d))
            d=list()

    return nums, code


def adjust_seq(seq,cigar):
    seq = list(seq)
    newseq = []
    nums,code = sep_cigar(cigar)
    for j,num in enumerate(nums):
        if code[j] == 'S':
            if j == 0:
                for z in range(num):
                    b = seq.pop(0)
            else:
                for z in range(num):
                    b = seq.pop(0)
                    newseq.append(b)
        elif code[j] == 'I':
            for z in range(num):
                b = seq.pop(0)
        elif code[j] == 'D':
            for z in range(num):
                newseq.append('-')
        elif code[j] == 'M':
            for z in range(num):
                b = seq.pop(0)
                newseq.append(b)
    
    return newseq

def get_consensus(alignment,pos):
    pos_bases = []
    for seq in alignment:
        if seq[pos] is not '-':
            pos_bases.append(seq[pos])
    pos_bases_counts = Counter(pos_bases).most_common()
    if len(pos_bases_counts)>1:
        prop = float(pos_bases_counts[0][1]) / (float(pos_bases_counts[1][1])+float(pos_bases_counts[0][1]))
        if prop > 0.7:
            c_base = pos_bases_counts[0][0]
        else:
            c_base = 'N'
    elif len(pos_bases_counts)==1:
        c_base = pos_bases_counts[0][0]
    else:
        c_base = 'N'
    
    return c_base
######################
contig_read_mappings=sys.argv[1]
folder_name=contig_read_mappings.split('/')[0]
file_name=contig_read_mappings.split('/')[1]
node_name=file_name.split('.')[0]

totalseq=[]
#get alignment
samfile=open(contig_read_mappings,'r')    #open file
for line in samfile:      #go through file
    if line.startswith('@'):
        continue
    elif len(line.strip())>0:
        splitline=line.split()
        flag = bin(int(splitline[1]))
        flagl = flag.split('b')
        if len(flagl[1])>=3:
            flag = flagl[1][-3]     #check for mapping
            if flag == '1':
                continue
        if len(flagl[1])>=5:
            flag = flagl[1][-5]     #check for reverse mapping
        else:
            flag = '0'
            
        cigar = splitline[5]
        pos = int(splitline[3])
        seq=splitline[9]
            
        seq = adjust_seq(seq,cigar)
        
        while pos>1:
            seq.insert(0, '-')
            pos=pos-1
        totalseq.append("".join(seq))

numreads = len(totalseq)
maxlen=0
for i,j in enumerate(totalseq):
    if len(j)>maxlen:
        maxlen = len(j)

for i,j in enumerate(totalseq):
    plus = maxlen - len(j)
    totalseq[i] = j + ('-' * plus)

#outfile = open('out_contig_seq_alignment.txt','w')
#for i,j in enumerate(totalseq):
#    outfile.write(j+"\n")
#outfile.close()

consensus = []
for i in range(len(totalseq[0])):
    c_base = get_consensus(totalseq,i)
    consensus.append(c_base)
consensus = ''.join(consensus)
consensus_seq = Seq(consensus)
consensus_seq_r = SeqRecord(consensus_seq, id=node_name)
SeqIO.write(consensus_seq_r,folder_name+'/'+node_name+'.fa', "fasta")

print numreads