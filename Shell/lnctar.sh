#!/bin/bash

outdir=/media/data7/hxy/Mus/DE

cd ~/software/LncTar

echo "--- LncTar: `date` ---"
perl LncTar.pl -p 1 -l ${outdir}/degCorLnc.oneline.fa -m ${outdir}/degCormRNA.transcript.oneline.fa -d -0.15 -s F -o ${outdir}/lnc_mRNA_interaction.txt
echo "--- Finished ---"
