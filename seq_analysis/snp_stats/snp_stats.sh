#!/bin/bash
#SBATCH --job-name snp_stats
#SBATCH --partition=long
#SBATCH --cpus-per-task=2
#SBATCH --output /shared/projects/elai_most/diversity_wasi/seq_analysis/snp_calling/slurm-%x_%j.log
#SBATCH --error /shared/projects/elai_most/diversity_wasi/seq_analysis/snp_calling/slurm-%x_%j.log


module load vcftools/0.1.16

cd /shared/projects/elai_most/diversity_wasi/seq_analysis/snp_stats

vcftools \
--gzvcf /shared/projects/afrob_seq/wgs_snp/vcf_bgi_2021/gbp_output/filtered_vcf_by_chrom/all_final_quality_SNPbiallelic_depth_missing_singletons_singletons.vcf.gz \
--weir-fst-pop vn_hybrids.txt \
--weir-fst-pop donor_AG.txt \
--fst-window-size 5000 \
--out hybrid_AG


# vcftools \
# --gzvcf all_final_quality_SNPbiallelic_depth_missing_singletons_singletons.vcf.gz \
# --TajimaD 5000 \
# --out pop_thin