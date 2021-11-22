#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda activate circ 
conda info --envs
 
workdir=/media/data7/hxy/Mus/yoelii/bln3/Fa
outdir=/media/data7/hxy/Mus/yoelii/bln3
dbdir=/media/data7/hxy/Mus/yoelii/bln/DB
thread=40
evalue=1e-10

cd $dbdir

for id in `ls $dbdir` 
do
	mkdir -p $outdir/rls_${id}
	cd $dbdir/${id}
	echo "*** $id START:`date` ***"

	for file in `ls ${workdir}/*.faa`
	do
		echo "* $(basename $file ".faa") *"
		blastp -query ${file} -db $id -num_threads $thread \
		-evalue $evalue \
		-out ${outdir}/rls_${id}/$(basename $file ".faa").rlt -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send evalue length pident bitscore'
		cd $outdir/rls_${id}
		Rscript ${outdir}/homo_identity.R $outdir/rls_${id}/$(basename $file ".faa").rlt
		cd $dbdir/${id}
	done
	cat $outdir/rls_${id}/*.tsv >> $outdir/rls_${id}/Final_${id}.tsv
done

