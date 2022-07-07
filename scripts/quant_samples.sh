#!/bin/bash
for fn in /media/studentsgh129/keldadrive/MSc/data/raw_data/*/;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i /media/studentsgh129/keldadrive/MSc/data/grch38_cdna_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
		 -p 8 \
         --recoverOrphans \
         -o quants/${samp}_quant\ \
         --seqBias \
         --gcBias \
         --validateMappings


done 