####### Guide for colocalisation #######

# BEFORE RUNNING PIPELINE ALL GWAS DATA NEEDS TO BE FORMATTED TO THE FOLLOWING HEADER NAMES
# chr
# position
# p 
# beta 
# se 
# eaf 
# rsid
# effect_allele
# other_allele
# N
# trait

####### DATA PREPROCESSING STEP #######

file0 = Extract_trait_regions.R

Extract_trait_regions will output file formatted_trait_regions
It will also extract the regions from the summary statistics files and output the results to GWAS data, this will speed up SUSIE by reducing time for loading full sum stats

####### STEP 1 ####### HYPRCOLOC

file1 = Hyprcoloc_base.R

We start by running the Hyprcoloc base script to identify regionals with SNPs in close proximity
The script will determine in regional prob > 0.8, if so, the SNPs in that region will be written to mult_coloc/regions/chr and the region will be the file name
Also a file with the prefix "res" will be produced containing the region colocalised and traits

####### STEP 2 ####### OBTAIN LD

file2 = convcf.sh
Extact SNPs outputted by Hyprcoloc, vcftools script will output SNP.recode file to LD_region folder
file3 = plinkld.sh
Conputes LD using plink for extracted SNPs by region

####### STEP 3 ####### OBTAIN SNP ESTIMATES

Run cat on files with prefix 'res' to obatin trait_regions.txt

file5 = suse_coloc_base.R
This will out the results of SUSIE to results


####### STEP 5  Extracting SNPs ####### 

####### OPTIONAL ONLY NECESSARY IF YOU NEED COLOCALISED SNPs TO BE AT THE SAME REGION #######

We need to go back to our conditionally independent SNPs from Mahajan and tag which regions colocalised.
We then take all SNPs in the colocalised region and extract from the traits in which they where colocalised

file6 = Extract_snps_from_sumstats