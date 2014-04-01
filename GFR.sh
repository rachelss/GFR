#!/bin/bash
#requires ref_genes.fa
#one folder per species containing pe reads as *1.fastq and *2.fastq
#do not have bam or sam files in these folders
#bash GFR.sh <number processors>

P=$1
FILELIST=($(ls */*1.fastq))
declare -a ALLFOLDERLIST=()
for F in "${FILELIST[@]}"; do ALLFOLDERLIST+=("$( dirname "${F}" )"); done       #array of directories containing FORMAT files
FOLDERLIST=( $(echo "${ALLFOLDERLIST[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )  #sorted unique list of folders with paired FORMAT files as array

bowtie2-build ref_genes.fa ref_genes    #build index for conserved contigs
rm */*bam
for FOLDER in "${FOLDERLIST[@]}"; do
    for FILE in "${FOLDER}"/*1.fastq; do
        bowtie2 -p "${P}" -N 1 --local -x ref_genes -1 ${FILE} -2 ${FILE/1.fastq/2.fastq} > >(tee ${FILE/1.fastq/}_stdout.log) 2> >(tee ${FILE/1.fastq/}_stderr.log >&2) | samtools view -Su -F 4 - | samtools sort - ${FILE/1.fastq/}  #output sorted bam file w/o unaligned reads - lots per folder
    done
done   
    
parallel -j $P "samtools merge {}/merged.bam {}/*.bam" ::: "${FOLDERLIST[@]}"       #merge bam files
parallel -j $P "samtools view {}/merged.bam > {}/merged.sam" ::: "${FOLDERLIST[@]}"     #convert to sam file
parallel -j $P "python separate_sam_by_contig.py {}/merged.sam" ::: "${FOLDERLIST[@]}"   #separate by gene -> lots of sam files per folder
rm */merged.sam
ls */*sam | parallel -j $P "python make_alignment_from_sam.py {}" #take sam file, do dumb consensus and make new conserved contig file -> lots of fa files per folder
parallel -j $P "cat {}/*fa > {}.fa"  ::: "${FOLDERLIST[@]}"     #concatenate fa files -> all species fa files in main folder
python genealign_from_allgenes.py  #make files for each gene containing seq for each species
for F in ${FOLDERLIST}; do rm ${F}.fa
    done      #delete species fasta files
for F in *fa; do
    if [ "${F}" != 'ref_genes.fa' ]; then   #align files of each gene
        mafft $F
    fi
done