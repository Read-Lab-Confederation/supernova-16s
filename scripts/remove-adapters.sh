#! /bin/bash
# Use BBDuk to remove adapters from input FASTQs and fastq-scan to for quality profiles
# /home/rpetit/projects/supernova-16s
# ls data/fastq/ | grep "R1" | sed 's/_L001_R1_001.fastq.gz//' | awk '{print $1" data/fastq/"$1"_L001_R1_001.fastq.gz data/fastq/"$1"_L001_R2_001.fastq.gz data/adapter-cleaned"}' | xargs -I {} sh -c 'scripts/remove-adapters.sh {}'
set -e
set -u

SAMPLE=$1
R1=$2
R2=$3
OUTDIR=$4

mkdir -p ${OUTDIR}

# Remove Adapters
STATS=${OUTDIR}/${SAMPLE}-bbduk.txt
FINAL_R1=${OUTDIR}/${SAMPLE}_R1.fastq.gz
FINAL_R2=${OUTDIR}/${SAMPLE}_R2.fastq.gz
 bbduk.sh -Xmx8g in=${R1} in2=${R2} out=${FINAL_R1} out2=${FINAL_R2} \
    ref=adapters k=23 ktrim=r mink=11 hdist=1 tpe=t tbo=t ordered=t tossjunk=t 2> ${STATS}
zcat ${FINAL_R1} | fastq-scan > ${OUTDIR}/${SAMPLE}_R1.json
zcat ${FINAL_R2} | fastq-scan > ${OUTDIR}/${SAMPLE}_R2.json
