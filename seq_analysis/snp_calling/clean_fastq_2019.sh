#!/bin/bash
#SBATCH --job-name clean_fastq_2019
#SBATCH --output /shared/projects/most_kmer/vcf_bgi_2021/config/slurm-%x.log
#SBATCH --error /shared/projects/most_kmer/vcf_bgi_2021/config/slurm-%x.log
#SBATCH --cpus-per-task=48
#SBATCH --partition=fast

module load fastp/0.23.1
module load cutadapt/3.1

fastq_2019_dir=/shared/projects/elai_most/data/fastq/raw_fastq_2019
read1_list_2019=$(find $fastq_2019_dir -wholename "**/*1.fastq.gz")
fastq_2021_dir=/shared/projects/elai_most/data/fastq/raw_fastq_bgi_2021
read1_list_2021=$(find $fastq_2021_dir -wholename "**/*1.fq.gz")
#read1_list=/shared/projects/elai_most/data/fastq/raw_fastq_2019/C047-MERGED/C047-MERGED_R1.fastq.gz
fastp_dir=/shared/projects/elai_most/data/fastq/fastq_2019_noadapter
outdir=/shared/projects/elai_most/data/fastq/fastq_2019_clean
mkdir -p $fastp_dir
mkdir -p $out_dir

for i1 in $read1_list_2019, $read1_list_2021
do
	fastq_dir=$(dirname $(dirname $i1))
	i2=$(echo ${i1/'1.f'/'2.f'}) # read 2
	o1=$(echo ${i1/$fastq_dir/$fastp_dir}) # output read 1 of fastp
	o2=$(echo ${i2/$fastq_dir/$fastp_dir}) # output read 2 of fastp
	mkdir -p $(dirname $o1) # create outdir for fastp
		
	fo1=$(echo ${i1/$fastq_dir/$outdir}) # output read 1 of cutadapt
	fo2=$(echo ${i2/$fastq_dir/$outdir}) # output read 2 of cutadapt
	mkdir -p $(dirname $fo1) # create outdir for cutadapt
	
	if [[ -f $fo1 && -f $fo2 ]] # if reads are TR samples
	then # skip cleaning
		echo -e "already cleaned $i1 and $i2"
	else # clean other samples
		# echo -e "cutadapt -j 0 -m 35 -q 20,20 -o $fo1 -p $fo2 $i1 $i2"
		# fastp -i $i1 -o $o1 -I $i2 -O $o2 --detect_adapter_for_pe --disable_quality_filtering --disable_length_filtering -w $SLURM_CPUS_PER_TASK -h $(dirname $o1)/fastp.html -j $(dirname $o1)/fastp.json # only trimming adapters
		cutadapt -j 0 -m 35 -q 20,20 -o $fo1 -p $fo2 $i1 $i2 # filter read ends by Q20
	fi
done
