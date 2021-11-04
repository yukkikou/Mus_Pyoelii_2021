#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda activate circ 
conda info --envs
 
workdir=/media/data7/hxy/Mus/yoelii/bln/Fa
#outdir=/media/data7/hxy/Mus/yoelii/bln/rls_falci
#dbdir=/media/data7/hxy/Mus/yoelii/bln/DB/Falci
#db=Falciparum
outdir=/media/data7/hxy/Mus/yoelii/bln2
dbdir=/media/data7/hxy/Mus/yoelii/bln/DB
thread=20
evalue=1e-10

cd $dbdir

for id in `ls $dbdir` 
do
	mkdir -p $outdir/rls_${id}
	cd $outdir/rls_${id}
	echo "*** $id START:`date` ***"

	for file in `ls ${workdir}/*.faa`
	do
		echo "* $(basename $file ".faa") *"
		#blastp -query ${file} -db $id -num_threads $thread \
		#-evalue $evalue \
		#-out ${outdir}/rls_${id}/$(basename $file ".faa").rlt -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send evalue length pident bitscore'
		Rscript ${outdir}/homo_identity.R $(basename $file ".faa").rlt
	done
	cat *.tsv >> Final_${id}.tsv
done

