### THESE SCRIPTS ARE MADE TO WORK IN MY ENVIRONMENT ###
### MODIFY THEM ACCORDING TO YOUR OWN DETAILS ###

Steps how the script works.

0) Make sure you have files:
-Create-PGS.R
-Extract-dosages-vcf.sh
-TBBMD-SNPs.csv (your own SNPs for creating the PRS)

0a) Run it by invoking:
bash Script-Vid.sh TBBMD-SNPs.csv

1) Extracts SNPs from the .vcf files and creates recoded .vcf files per
chromosome in "OUTPUT" folder
1a) Extracts SNPs from the "TBBMD-SNPs.csv"
1b) If you wish to filter for certain SNPs i.e. exclude them, delete them from
the "TBBMD-SNPs.csv" file

2) Waits until all the extraction is done
2a) Looks at "qstat" information for your MSN

3) Extracts the dosage information for each variant

4) Joins all the vcf files into a single one

5) Runs the R script to extract polygenic scores
5a) First harmonizes the dosage information with the GWAS data i.e. changes
dosages so that the dosage information is aligned with effect allele information.
5b) Recodes the GWAS data to alleles so that they are effect increasing i.e. if
variant "A>G" has effect -0.5, it recodes it into variant "G>A" with effect
0.5, while at the same time changes the dosages "0 becomes 2, 1 remains 1, 2
becomes 0 etc."
5c) Extracts number of trait increasing alleles (unweighted approach) and trait
increasing alleles multiplied by their respective weights (effects; weighted
approach)

6) Polygenic scores are written to the (trait)_SNPS-(study)-DOSAGES.txt file
