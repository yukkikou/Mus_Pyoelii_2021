#!/bin/bash
#SBATCH -A ylab
#SBATCH -J combination
#SBATCH -p single-ib
#SBATCH -D /media/work/hxy/Mus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=35
#SBATCH -o /home/hxy/log/Py_com-%A_%a.out
#SBATCH -e /home/hxy/log/Py_com-%A_%a.out
#SBATCH --mem 80G

#######################################################
gtf=/media/share/node13/disk3/ylab/hxy/Mus/annotation/combination/para_mus.gtf
genome=/media/share/node13/disk3/ylab/hxy/Mus/annotation/combination/para_mus_genome.fa
fastq=/media/share/node13/disk1/ylab/hxy/RNA/C328/cleandata
outdir=/media/share/node13/disk3/ylab/hxy/Mus/yoelii
workdir=/media/work/hxy/Mus/yoelii
tmpdir=${workdir}/_tmp
prefix=Pyoelii
thread=35
list=/media/share/node13/disk3/ylab/hxy/Mus/C328/sample.fa.list

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
echo "Max index is ${SLURM_ARRAY_TASK_MAX}";
echo "Min index is ${SLURM_ARRAY_TASK_MIN}";
echo now files: `ls ./`;
echo -----------------------
#########################################################

mkdir -p ${workdir}
mkdir -p ${tmpdir}

#mkdir -p ${workdir}/hisat2
#mkdir -p ${workdir}/stringtie
#mkdir -p ${workdir}/mRNA
#mkdir -p ${workdir}/circ/CIRIquant

mkdir -p ${outdir}/hisat2
mkdir -p ${outdir}/stringtie
mkdir -p ${outdir}/mRNA
mkdir -p ${outdir}/index
 
hisat2_index=${outdir}/index/hisat_index

###################################################################################################################
cd ${workdir}

source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda activate circ
conda info --envs

##########################################################mRNA###############################################
#mapping
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
	echo "--- ${pre} end: `date` ---"
done

#assemble
ls ${outdir}/hisat2/*.sam | while read id; do ( samtools view -bF 12 --threads ${thread} -o ${outdir}/hisat2/$(basename ${id} ".sam").bam ${id} ); done
ls ${outdir}/hisat2/*.bam | while read id; do ( samtools sort -O bam -@ ${thread} -o ${outdir}/hisat2/$(basename ${id} ".bam").sort.bam ${id} ); done
echo "--- sam to bam: `date` ---"

rm ${outdir}/hisat2/*.sam
echo "--- sam removal: `date` ---"

ls ${workdir}/hisat2/*.sort.bam | xargs -i samtools index {}
echo "--- index: `date` ---"

#ls ${workdir}/hisat2/*.sort.bam | while read id; do ( samtools flagstat -@ ${thread} $id > ${workdir}/hisat2/$(basename ${id} ".sort.bam").flagstat ); done
#echo "--- flagstat: `date` ---"

ls ${outdir}/hisat2/*.sort.bam | while read id
do
	echo "--- stringtie start: `date` ---" && \
	time stringtie -e -p ${thread} \
	-G $gtf -o ${outdir}/stringtie/$(basename ${id} ".sort.bam").gtf ${id}
	echo "--- stringtie end: `date` ---"
done

ls ${outdir}/stringtie/*.gtf > ${outdir}/stringtie/gtf.list
echo "--- stringtie merge: `date` ---" && \
	time stringtie --merge -p ${thread} -G $gtf -l ${prefix} -o ${outdir}/stringtie/merged.gtf ${outdir}/stringtie/gtf.list

# prepDE.py3 is in ${outdir}
python ${outdir}/prepDE.py3 -g ${outdir}/mRNA/gene_count_matrix.csv -t ${outdir}/mRNA/transcript_count_matrix.csv -i ${outdir}/stringtie/gtf.samlist

echo "--- Finished: `date` ---"
echo 'now files:'
ls ./*
