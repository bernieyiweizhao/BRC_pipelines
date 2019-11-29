#!/bin/bash

fastq_path=$(readlink -f $1)
out_prefix=$2
ref=$3

mkdir -p $out_prefix

sample=$(echo $fastq_path | rev | cut -d "/" -f 1 | rev)
project=$(echo $fastq_path | rev | cut -d "/" -f 3 | rev)
lab=$(echo $fastq_path | rev | cut -d "/" -f 4 | rev)

echo $sample,$project,$lab

R1=$(ls $fastq_path/*/${sample}_S*_L*_R1_001.fastq.gz | paste -sd "," -)
R2=$(ls $fastq_path/*/${sample}_S*_L*_R2_001.fastq.gz | paste -sd "," -)

command="STAR --runThreadN 8 --genomeDir $ref --readFilesIn $R1 $R2 --outFileNamePrefix $out_prefix --readFilesCommand gunzip -c --outSAMtype BAM Unsorted"
$command

