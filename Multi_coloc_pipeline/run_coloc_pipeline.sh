#!/bin/bash

#########

## I think a superior approach would be to have the bash arguments as $1 = number of GWAs, $2 file containing conditionally independent lead SNPs of primary GWAS $3 path of primary GWAS, $4, $5 etc ... Secondary GWASs

## It would also be best to use this R optarg interface that Vid used in a previous project. Maybe for next time.

# $1 - Path to directory
# 
# 

######## Extract regions from secondary GWAS ########
# Extract_trait_regions will output file formatted_trait_regions
# The R file needs to be manually changed
Rscript Extract_trait_regions.R

# It will also extract the regions from the summary statistics files and output the results to GWAS data, this will speed up SUSIE by reducing time for loading full sum stats.
# Maybe this should be moved to the begining for efficency.

######## Run Hypercoloc ########

# Rscript Hyprcoloc_base.R

######## OBTAIN LD #############
#  Extact SNPs outputted by Hyprcoloc, vcftools script will output SNP.recode file to LD_region folder

convcf.sh $DIRECTORY

wait

# Conputes LD using plink for extracted SNPs by region

plinkld.sh $DIRECTORY + LD_regions

####### OBTAIN SNP ESTIMATES ####### 

# Run cat on files with prefix 'res' to obatin trait_regions.txt
cat $1"/res_*" > trait_regions.txt


Rscript suse_coloc_base.R

#This will out the results of SUSIE to results

####### OPTIONAL ONLY NECESSARY IF YOU NEED COLOCALISED SNPs TO BE AT THE SAME REGION #######

We need to go back to our conditionally independent SNPs from Mahajan and tag which regions colocalised.
We then take all SNPs in the colocalised region and extract from the traits in which they where colocalised

file6 = Extract_snps_from_sumstats