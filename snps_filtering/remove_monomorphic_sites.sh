#!/bin/bash
#SBATCH --output /data3/projects/vietcaf/baotram/slurm-%x_%j.log
#SBATCH --error /data3/projects/vietcaf/baotram/slurm-%x_%j.log

#module load bioinfo/gatk/4.1.4.1
module load bioinfo/vcftools/0.1.16

mydir=/data3/projects/vietcaf/baotram
scratch=$(mktemp -dp /scratch)

vcf_input=vietcaf_general_DP_missing_indv_singletons_Filtered.recode.vcf
vcf_biallelic=vietcaf_general_DP_missing_indv_singletons_Filtered_biallelic

cd $scratch

# transfer input files to scratch
rsync -vauP $mydir/$vcf_input  ./
#rsync -vaurP /data3/projects/vietcaf/reference ./

# remove monomorphic sites

vcftools --vcf $vcf_input --non-ref-ac 1 --max-non-ref-af 0.999 --out $vcf_biallelic --recode

# transfer output files to nas
rsync -vauP $vcf_biallelic.recode.vcf $mydir/

cd ../
rm -rf $scratch
