#!/bin/bash
#SBATCH --output /data3/projects/vietcaf/baotram/vcf_by_chrom/slurm-%x_%j.log
#SBATCH --error /data3/projects/vietcaf/baotram/vcf_by_chrom/slurm-%x_%j.log
#SBATCH --exclude=node19

module load bioinfo/vcftools/0.1.16

mydir=/data3/projects/vietcaf/baotram
scratch=$(mktemp -dp /scratch)

vcf_file=vietcaf_general_DP_missing_indv_singletons_Filtered_biallelic.recode.vcf

cd $scratch

# transfer input file to scratch
rsync -vauP $mydir/$vcf_file ./
#rsync -vaurP /data3/projects/vietcaf/reference ./

# extract snps by chromosome
for i in {01..11}
do 
	vcftools --vcf $vcf_file --chr CC1.8.Chr$i --out vietcaf_final_chr$i --recode
	#gatk SelectVariants \
	#-R reference/CC1.8_v2_pseudomolecule_cat.fa \
	#-V $vcf_file \
	#-L CC1.8.Chr$i \
	#-O vietcaf_final_chr$i.vcf
done


# transfer output files to nas
rsync -vauP vietcaf_final_chr*.recode.vcf $mydir/vcf_by_chrom

cd ../
rm -rf $scratch
