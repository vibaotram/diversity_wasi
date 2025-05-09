# diversity_wasi
Scripts and some data available for [paper]()

## genotyping
- snp_panels: SNP data before and after filtering
- lgc_results: raw data from LGC
- raw_genotype: genotype data of all individuals (including the ones not on the paper)
- genetic_diveristy.R: all the genotyping analyses


## seq_analysis
- snp_calling: scripts and config files to run the [snp calling workflow](https://github.com/vibaotram/gbp_variantcalling) (the workflow is on another repository)
- snp_filtering: scripts for filtering snps after calling
- snp_stats: scripts to perform some statistic analyses on the snp data
- elai: config files and scripts to run [elai workflow](https://github.com/vibaotram/snakelai) (the workflow on another repository), and scripts to merge elai results from different runs
