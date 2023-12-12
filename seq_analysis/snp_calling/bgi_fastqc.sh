#!/bin/bash
#SBATCH --job-name bgi_fastqc
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --output /shared/projects/most_kmer/vcf_bgi_2021/config/slurm-%x_%j.log
#SBATCH --error /shared/projects/most_kmer/vcf_bgi_2021/config/slurm-%x_%j.log

module load singularity
module load conda
module load snakemake/6.12.3

cd /shared/projects/most_kmer/scripts/gbp_variantcalling

snakemake fastqc --use-conda -p --verbose --rerun-incomplete \
--jobs 60 --latency-wait 30 --nolock --no-temp \
--configfile /shared/projects/most_kmer/vcf_bgi_2021/config/config.yaml \
--cluster "python3 slurm_wrapper.py" \
--cluster-status "python3 slurm_status.py" 
