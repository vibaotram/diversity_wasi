#!/bin/bash
#SBATCH --output /data3/projects/vietcaf/baotram/slurm-%x_%j.log
#SBATCH --error /data3/projects/vietcaf/baotram/slurm-%x_%j.log

#module load bioinfo/gatk/4.1.4.1
module load bioinfo/vcftools/0.1.16

mydir=/data3/projects/vietcaf/baotram
scratch=$(mktemp -dp /scratch)

vcf_input=vietcaf_general_DP_missing_indv_singletons_Filtered_biallelic.recode.vcf
vcf_out=vietcaf_final_biallelic_random_10perc.vcf
cd $scratch

# transfer input files to scratch
rsync -vauP $mydir/$vcf_input  ./
#rsync -vaurP /data3/projects/vietcaf/reference ./

# filter random SNPs

gatk SelectVariants \
     -R $mydir/reference/CC1.8_v2_pseudomolecule_cat.fa \
     -V $vcf_input \
     --select-random-fraction 0.1 \
     -O $vcf_out

# transfer output files to nas
rsync -vauP $vcf_biallelic.recode.vcf $mydir/

cd ../
rm -rf $scratch
