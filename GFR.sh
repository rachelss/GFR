#!/bin/bash
#requires ref_genes.fa
#one folder per species containing pe reads as *R1*.fastq and *R2*.fastq
#do not have bam or sam files in these folders
#bash GFR.sh

help() {
    echo "
Default settings can be changed using the following flags:

    -p : use to specify the number of processors available
    -s : use to specify the steps to skip: 1 skips aligning to the genes
    Example command: bash GFR.sh -s 1"
    }

P=1
SKIP=0

while getopts p:s:h option
do
case "${option}"
    in
        p) P=${OPTARG};;
        s) SKIP=${OPTARG};;
        h) help; exit;;
        \? ) echo "Unknown option" >&2; exit 1;;
        esac
done

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"    #where is the main program file? 

FILELIST=($(ls */*R1*.fastq))
declare -a ALLFOLDERLIST=()
for F in "${FILELIST[@]}"; do ALLFOLDERLIST+=("$( dirname "${F}" )"); done       #array of directories containing FORMAT files
FOLDERLIST=( $(echo "${ALLFOLDERLIST[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ') )  #sorted unique list of folders with paired FORMAT files as array

if [ "$SKIP" -lt 1 ]; then
    bowtie2-build ref_genes.fa ref_genes    #build index for conserved contigs
    rm */*bam
    for FOLDER in "${FOLDERLIST[@]}"; do
        for FILE in "${FOLDER}"/*R1*.fastq; do
            bowtie2 -p "${P}" -N 1 --local -x ref_genes -1 ${FILE} -2 $( echo ${FILE}|sed 's/R1/R2/' ) > >(tee ${FILE/.fastq/}_stdout.log) 2> >(tee ${FILE/.fastq/}_stderr.log >&2) | samtools view -Su -F 4 - | samtools sort - ${FILE/.fastq/}  #output sorted bam file w/o unaligned reads - lots per folder
        done
    done   
fi

rm */merged.bam
parallel -j $P "samtools merge {}/merged.bam {}/*.bam" ::: "${FOLDERLIST[@]}"       #merge bam files
for FOLDER in "${FOLDERLIST[@]}"; do if [ ! -f ${FOLDER}/merged.bam ]; then                   #samtools won't merge one file - need to copy
    A=${FOLDER}/*bam
    ln $(echo "$A") ${FOLDER}/merged.bam
fi
done

parallel -j $P "samtools view {}/merged.bam > {}/merged.sam" ::: "${FOLDERLIST[@]}"     #convert to sam file
parallel -j $P "python ${DIR}/separate_sam_by_contig.py {}/merged.sam" ::: "${FOLDERLIST[@]}"   #separate by gene -> lots of sam files per folder
rm */merged.sam
rm */*fa
rm */*counts
mv ref_genes.fa ref_genes.fasta
rm *fa
for i in */*sam; do python ${DIR}/make_alignment_from_sam.py $i #take sam file, do dumb consensus and make new conserved contig file -> lots of fa files per folder
parallel -j $P "cat {}/*fa > {}.fa"  ::: "${FOLDERLIST[@]}"     #concatenate fa files -> all species fa files in main folder
python ${DIR}/genealign_from_allgenes.py  #make files for each gene containing seq for each species
mv ref_genes.fasta ref_genes.fa
for F in "${FOLDERLIST[@]}"; do rm ${F}.fa
    done      #delete species fasta files
for F in *fa; do
    if [ "${F}" != 'ref_genes.fa' ]; then   #align files of each gene
        mafft $F > ${F/.fa/_align.fa}
    fi
done
