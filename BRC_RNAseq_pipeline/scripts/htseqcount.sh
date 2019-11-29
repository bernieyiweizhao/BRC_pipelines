#/bin/bash

bam=$1
gtf=$2
count=$3

out_dir=$(echo $count | rev | cut -d "/" -f 2- |rev)
mkdir -p $out_dir
htseq-count -f bam -s reverse $bam $gtf --additional-attr=gene_name > $count
