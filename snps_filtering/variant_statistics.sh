#!/bin/bash
#SBATCH --output /shared/projects/most_kmer/vcf_bgi_2021/gbp_output/vcf_by_chrom/slurm-%x_%j.log
#SBATCH --error /shared/projects/most_kmer/vcf_bgi_2021/gbp_output/vcf_by_chrom/slurm-%x_%j.log

module load bcftools/1.9
module load vcftools/0.1.16

mydir=/shared/projects/most_kmer/vcf_bgi_2021/gbp_output/vcf_by_chrom
cd $mydir

vcf_output=all_chrom_raw.vcf

vcf_inputs=vcf_inputs=$(echo $(ls *.vcf))

bcftools concat -o $vcf_output $vcf_inputs

vcftools --vcf $vcf_output --site-quality --site-mean-depth --missing-site --SNPdensity 10


