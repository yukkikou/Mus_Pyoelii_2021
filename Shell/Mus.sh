#!/bin/bash
#SBATCH -A ylab
#SBATCH -J MusPipeLine
#SBATCH -p compute
#SBATCH -D /media/work/hxy/Mus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=100
#SBATCH -o /home/hxy/log/Mus_miRNA-%A_%a.out
#SBATCH -e /home/hxy/log/Mus_miRNA-%A_%a.out
#SBATCH --mem 240G

#######################################################
gtf=/media/share/node13/disk3/ylab/hxy/annotation/Mus_musculus.GRCm39.103.gtf
genome=/media/share/node13/disk3/ylab/hxy/genome/Mus_musculus.GRCm39.dna.toplevel.fa
cds=/media/share/node13/disk3/ylab/hxy/genome/Mus_musculus.GRCm39.cds.all.fa
cdna=/media/share/node13/disk3/ylab/hxy/genome/Mus_musculus.GRCm39.cdna.all.fa
ncrna=/media/share/node13/disk3/ylab/hxy/genome/Mus_musculus.GRCm39.ncrna.fa
fastq=/media/share/node13/disk1/ylab/hxy/RNA/C328/cleandata
outdir=/media/share/node13/disk3/ylab/hxy/C328
workdir=/media/work/hxy/Mus
tmpdir=${workdir}/_tmp
prefix=Mus
thread=100
list=/media/share/node13/disk3/ylab/hxy/C328/sample.fa.list

Pfam_db=/media/share/node13/disk3/ylab/hxy/Pfam_db

mature_miRNA=/media/share/node13/disk3/ylab/hxy/miRNA/mature.fa
hairpin_miRNA=/media/share/node13/disk3/ylab/hxy/miRNA/hairpin.fa
miRNA_gff=/media/share/node13/disk3/ylab/hxy/miRNA/mmu.gff3
miRNA_list=/media/share/node13/disk3/ylab/hxy/C328/miRNA/miRNA.list
mature_mmu_miRNA=/media/share/node13/disk3/ylab/hxy/miRNA/mature.mmu.fa
hairpin_mmu_miRNA=/media/share/node13/disk3/ylab/hxy/miRNA/hairpin.mmu.fa

########################################################
echo ---Time is--- `date`;
echo WorkDirectory is `pwd`;
echo This job runs on $SLURM_CPUS_PER_TASK CPUs;
echo SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID;
echo SLURM_CPUS_ON_NODE is $SLURM_CPUS_ON_NODE of CPUs on the allocated node;
echo SLURM_JOB_ACCOUNT is $SLURM_JOB_ACCOUNT;
echo SLURM_JOB_NODELIST is $SLURM_JOB_NODELIST;
echo SLURM_JOB_NUM_NODES is $SLURM_JOB_NUM_NODES;
echo SLURM_JOB_PARTITION is $SLURM_JOB_PARTITION;
echo SLURM_SUBMIT_DIR is $SLURM_SUBMIT_DIR;
echo SLURM_SUBMIT_HOST is $SLURM_SUBMIT_HOST;
echo SLURM_TASKS_PER_NODE is $SLURM_TASKS_PER_NODE;
echo SLURMD_NODENAME is $SLURMD_NODENAME;
echo "This is job #${SLURM_ARRAY_JOB_ID}, with parameter ${input[$SLURM_ARRAY_TASK_ID]}";
echo "There are ${SLURM_ARRAY_TASK_COUNT} task(s) in the array.";
echo "  Max index is ${SLURM_ARRAY_TASK_MAX}";
echo "  Min index is ${SLURM_ARRAY_TASK_MIN}";
echo now files: `ls ./`;
echo -----------------------
#########################################################
#conda activate
#conda activate circ
#conda install -y hisat2
#conda install -y stringtie
#conda install -y gffcompare

#conda activate miRNA
#conda install bowtie
#conda install fastqc
#conda install fastx_toolkit

mkdir -p ${workdir}
mkdir -p ${tmpdir}
mkdir -p ${workdir}/stringtie
mkdir -p ${workdir}/ballgown
mkdir -p ${workdir}/mRNA
mkdir -p ${workdir}/miRNA
mkdir -p ${workdir}/miRNA/QC
mkdir -p ${workdir}/miRNA/QC/filter
mkdir -p ${workdir}/circ/CIRIquant

#mkdir -p ${outdir}/QC
#mkdir -p ${outdir}/index
#mkdir -p ${outdir}/hisat2
#mkdir -p ${outdir}/stringtie
#mkdir -p ${outdir}/mRNA
#mkdir -p ${outdir}/miRNA
#mkdir -p ${outdir}/miRNA/QC
#mkdir -p ${outdir}/miRNA/QC/filter
#mkdir -p ${outdir}/circRNA/CIRIquant

hisat2_index=${outdir}/index/hisat_index
bwa_index=${outdir}/index/bwa_index
bowtie2_index=${outdir}/index/bowtie2_index
mature_miRNA_bowtie_index=${outdir}/index/mature_miRNA_bowtie_index
hairpin_miRNA_bowtie_index=${outdir}/index/hairpin_miRNA_bowtie_index


###################################################################################################################
cd ${workdir}

source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda info --envs

##########################################################mRNA(C328)###############################################
###fq
echo '--- fastqc ---'
function getdir(){
    for element in `ls $1`
      do
        dir_or_file=$1"/"$element
    if [ -d $dir_or_file ]
      then
        getdir $dir_or_file
      else
        echo $dir_or_file     
    fi
done
}
getdir $fastq | grep gz |  while read id; do fastqc -o QC/ -t 6 $id ; done
cp QC/* ${outdir}/QC
getdir $fastq | grep gz > ${outdir}/sample.fa.list
#change sample.fa.list format!

###mapping
echo '--- building index ---'
hisat2-build -p $thread $genome ${outdir}/index/hisat_index 

echo '--- mapping with hisat2 ---'
cat $list | while read id
do
        arr=(${id})
        pre=${arr[0]}
        fq1=${arr[1]}
        fq2=${arr[2]}
	hisat2 -x ${hisat2_index} -p ${thread} --sensitive --no-discordant --no-mixed -t \
	--summary-file ${outdir}/hisat2/${pre}.hisat2sam.txt \
	-1 ${fq1} -2 ${fq2} -S ${outdir}/hisat2/${pre}.sam
	echo ---------------${pre}-----------------
	echo `date`
done

###assemble
ls ${outdir}/hisat2/*.sam | while read id; do ( samtools view -bF 12 --threads ${thread} -o ${outdir}/hisat2/$(basename ${id} ".sam").bam ${id} ); done
ls ${outdir}/hisat2/*.bam | while read id; do ( samtools sort -O bam -@ ${thread} -o ${outdir}/hisat2/$(basename ${id} ".bam").sort.bam ${id} ); done
echo '--- s2b ---'
echo `date`
rm ${outdir}/hisat2/*.sam
echo '--- sam removal ---'
ls ${outdir}/hisat2/*.sort.bam | xargs -i samtools index {}
echo '--- index ---'
echo `date`
ls ${outdir}/hisat2/*.sort.bam | while read id; do ( samtools flagstat -@ ${thread} $id > ${outdir}/hisat2/$(basename ${id} ".sort.bam").flagstat ); done
echo '--- flagstat ---'
echo `date` 

ls ${outdir}/hisat2/*.sort.bam | while read id
do
	echo "--- stringtie start: `date` ---" && \
	time stringtie -f 0.3 -j 3 -c 5 -g 100 -s 10000 -p ${thread} \
	-G $gtf -o ${workdir}/stringtie/$(basename ${id} ".sort.bam").gtf ${id}
	echo "--- stringtie end: `date` ---"
done

ls ${outdir}/stringtie/*.gtf > ${outdir}/stringtie/gtf.list
echo "--- stringtie merge: `date` ---" && \
	time stringtie --merge -p ${thread} -G $gtf -l ${prefix} -o ${workdir}/stringtie/merged.gtf ${outdir}/stringtie/gtf.list

###gffcompare
genome_sm=/media/share/node13/disk3/ylab/hxy/genome/Mus_musculus.GRCm39.dna_sm.toplevel.fa

gffcompare -r $gtf -s ${genome_sm} -o ${workdir}/stringtie/gffcompare ${outdir}/stringtie/merged.gtf
prepDE.py -g ${workdir}/mRNA/gene_count_matrix.csv -t ${workdir}/mRNA/transcript_count_matrix.csv -i ${outdir}/stringtie/gtf.samlist

cp ${workdir}/stringtie/* ${outdir}/stringtie/
cp ${workdir}/mRNA/* ${outdir}/mRNA/

######################################################miRNA(C355)############################################
####fq
conda activate miRNA
conda info --envs

echo '--- miRNA part ---'
cat ${miRNA_list} | while read id
do
        arr=(${id})
        pre=${arr[0]}
        fq=${arr[1]}
        echo "--- ${pre} ---"
	echo `date`
	echo $fq
	fastqc -o ${workdir}/miRNA/QC --nogroup -t ${thread} $fq
	zcat $fq | fastq_quality_filter -v -q 20 -p 80 -i - -o ${pre}_tmp
	fastx_trimmer -v -m 15 -i ${pre}_tmp -o -z ${workdir}/miRNA/${pre}.clean.fq.gz
	fastqc -o ${outdir}/miRNA/QC/filter --nogroup -t ${thread} ${outdir}/miRNA/${pre}.clean.fq.gz
	rm ${pre}_tmp
done

####two way to quantity
####method1: miRBase
perl -alne '{if(/^>/){if(/Mus musculus/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' ${hairpin_miRNA} > $(dirname ${hairpin_miRNA})/hairpin.mmu.fa
perl -alne '{if(/^>/){if(/Mus musculus/){$tmp=1}else{$tmp=0}};next if $tmp!=1;s/U/T/g if !/>/;print }' ${mature_miRNA} >  $(dirname ${hairpin_miRNA})/mature.mmu.fa

mature_mmu_miRNA=$(dirname ${hairpin_miRNA})/mature.mmu.fa
hairpin_mmu_miRNA=$(dirname ${hairpin_miRNA})/hairpin.mmu.fa
mature_miRNA_bowtie_index=${outdir}/index/mature_miRNA_bowtie_index
hairpin_miRNA_bowtie_index=${outdir}/index/hairpin_miRNA_bowtie_index

###bowtie
bowtie-build --threads ${thread} ${mature_mmu_miRNA} ${outdir}/index/mature_miRNA_bowtie_index
bowtie-build --threads ${thread} ${hairpin_mmu_miRNA} ${outdir}/index/hairpin_miRNA_bowtie_index
echo '--- index done ---'


cat ${miRNA_list} | while read id
do
        arr=(${id})
        pre=${arr[0]}
        fq=${arr[1]}
        clean=${arr[2]}
        echo "--- ${pre}: `date` ---"
        bowtie -n 0 -m1 -f --best --strata -x ${mature_miRNA_bowtie_index} ${clean} -S ${workdir}/miRNA/${pre}.mature.sam
        bowtie -n 0 -m1 -f --best --strata -x ${hairpin_miRNA_bowtie_index} ${clean} -S ${workdir}/miRNA/${pre}.hairpin.sam
done
echo '--- assignment end ---'

ls *.sam | while read id ;do (samtools sort -O bam -@ 5  -o $(basename ${id} ".sam").bam   ${id}); done
rm *.sam

ls *.bam | xargs -i samtools index {}
ls *.bam | while read id ;do (samtools idxstats ${id} > ${id}.txt ); done
#samtools view  matrue.bam | cut -f 3 | sort | uniq  -c
echo '--- miRNA end ---'

#########################################################circRNA####################################################
#CIRI2
#bwa index -p ${outdir}/index/bwa_index ${genome}
#echo 'mapping with bwa'
#cat $list | while read id
#do
#       arr=(${id})
#       pre=${arr[0]}
#       fq1=${arr[1]}
#       fq2=${arr[2]}
#       bwa mem -t ${thread} ${bwa_index} ${fq1} ${fq2} > ${pre}.sam
#	perl ${outdir}/circRNA/CIRI2.pl -I ${pre}.sam -O ${outdir}/circRNA/${pre}.circ.txt -F $genome -G ${outdir}/circRNA/${pre}.log -high -T ${thread}
#	#rm ${pre}.sam
#	echo ---------------${pre}-----------------
#       	echo `date`
#done

#CIRIquant
#need more test!!!
#awk '{print $1}' ${outdir}/sample.fa.list | xargs -n1 mkdir ${outdir}/circRNA/CIRIquant/

#echo '----CIRIquant------'
#cat $list | while read id
#do
#       	arr=(${id})
#       	pre=${arr[0]}
#       	fq1=${arr[1]}
#       	fq2=${arr[2]}		
#	CIRIquant -t 20 -v -1 ${fq1} -2 ${fq2} --config ${outdir}/circRNA/CIRIquant/pm1.yml -p ${pre} --circ ${outdir}/circRNA/${pre}.circ.txt --tool CIRI2 -o ${outdir}/circRNA/CIRIquant/${pre}
#done
##step 1
#prep_CIRIquant -i ${outdir}/circRNA/CIRIquant/sample.lst \
#                 --lib ${outdir}/circRNA/CIRIquant/library_info.csv \
#                 --circ ${outdir}/circRNA/CIRIquant/circRNA_info.csv \
#                 --bsj ${outdir}/circRNA/CIRIquant/circRNA_bsj.csv \
#                 --ratio ${outdir}/circRNA/CIRIquant/circRNA_ratio.csv
#echo '---step 1 end------'
#These count matrices (CSV files) can then be imported into R for use by DESeq2 and edgeR
#(using the DESeqDataSetFromMatrix and DGEList functions, respectively).
##step 2
#The output of StringTie should locate under output_dir/gene/prefix_out.gtf
#use prepDE.py from stringTie (conda env circ) to generate the gene count matrix for normalization.
#prepDE.py -i ${outdir}/circRNA/CIRIquant/sample_gene.lst -g ${outdir}/circRNA/CIRIquant/gene_count_matrix.csv -t ${outdir}/circRNA/CIRIquant/transcript_count_matrix.csv -v
#echo '---step 2 end------'
#use gene_count_matrix.csv generated under current working directory for further analysis.
##step 3
#R packages: edgeR(BioManager) optparse statmod
#CIRI_DE_replicate --lib  ${outdir}/circRNA/CIRIquant/library_info.csv \
#            --bsj  ${outdir}/circRNA/CIRIquant/circRNA_bsj.csv \
#            --gene ${outdir}/circRNA/CIRIquant/gene_count_matrix.csv \
#            --out  ${outdir}/circRNA/CIRIquant/circRNA_de.tsv
#echo '---step 3 end------'
#the output results is unfiltered
#echo '----circquant end-------'

#rm -r QC/
echo '----all-----end----'
echo 'now files:'
ls ./*
