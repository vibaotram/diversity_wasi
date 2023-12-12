#!/bin/bash
#SBATCH --job-name bgi_gbp_2021
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --output /shared/projects/afrob_seq/wgs_snp/vcf_bgi_2021/config/slurm-%x_%j.log
#SBATCH --error /shared/projects/afrob_seq/wgs_snp/vcf_bgi_2021/config/slurm-%x_%j.log

module load singularity
module load conda
module load snakemake/6.12.3

# source activate /shared/projects/vietcaf/conda/env/snakemake

cd /shared/projects/afrob_seq/wgs_snp/scripts/gbp_variantcalling

snakemake gbp_variantcalling --use-conda -p --verbose --rerun-incomplete \
--jobs 60 --latency-wait 30 --nolock --keep-going --notemp \
--configfile /shared/projects/afrob_seq/wgs_snp/vcf_bgi_2021/config/config.yaml \
--cluster "python3 slurm_wrapper.py" \
--cluster-status "python3 slurm_status.py" \
--forcerun GenotypeGVCFs # use this option when adding new sequences to an existing pipeline
