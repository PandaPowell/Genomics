#!bin/bash

# Make sure reference genome is in bed format if not use
# wget -qO- ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz | gunzip -c - | vcf2bed --sort-tmpdir=${PWD} --max-mem=2G - > hg19.dbSNP150.bed

# Next convert your file of chr:pos to a bed format and sort using awk, requires chr, pos in columns 1 & 2, no headers

awk -vOFS="\t" '{ print $1, $2-1, $2, $3, $4,$5,$6,$7,$8,$9,$10; }' positions.txt | sort-bed - > positions.bed

# Next you map the positions between the two bed files

bedmap --echo --echo-map-id --delim '\t' positions.bed hg19.dbSNP150.common.bed > answer.bed
