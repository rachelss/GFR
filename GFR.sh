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
    for FILE in *1.fastq; do
        bowtie2 -p "${P}" -N 1 --local -x ref_genes -1 ${FILE} -2 ${FILE/1.fastq/2.fastq} > >(tee ${FILE/1.fastq/}_stdout.log) 2> >(tee ${FILE/1.fastq/}_stderr.log >&2) | samtools view -Su -F 4 - | samtools sort - ${FILE/1.fastq/}  #output sorted bam file w/o unaligned reads
    done
done   
    
parallel -j $P "samtools merge {}/merged.bam {}/*.bam" ::: "${FOLDERLIST[@]}"
parallel -j $P "samtools view {}/merged.bam > {}/merged.sam" ::: "${FOLDERLIST[@]}"
parallel -j $P "python separate_sam_by_contig.py {}/merged.sam" ::: "${FOLDERLIST[@]}"   #separate by gene
rm */merged.sam
ls */*sam | parallel -j $P "python make_alignment_from_sam.py {}" #take sam file, do dumb consensus and make new conserved contig file
python concatenate.py   #make files for each gene containing each species