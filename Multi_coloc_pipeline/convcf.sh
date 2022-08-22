#!/bin/bash

INPUTDIR=$1
OUTDIR=$1"/LD_regions"

mkdir $OUTDIR

# Extract SNPs from VCF files
#### important !!!! for chrm1 add noIBD #######

echo "Extracting SNPs from .vcf files"

for CHR in 1
do
  cd $INPUTDIR/chr_"$CHR"
  
  for FILE in *
  do
  
    /opt/sge/bin/lx-amd64/qsub -cwd -V -b y -q long.q vcftools --gzvcf "/home/055442/reference_panels/HRCr1.1/VCF/HRC.r1-1.EGA.GRCh37.chr"$CHR".haplotypes.noIBD.vcf.gz" \
                              --positions $FILE \
                              --maf 0.01 \
                              --recode \
                              --out $OUTDIR/$FILE
    echo "Job submitted"
  done
done


for CHR in {2..22}
do
  cd $INPUTDIR/chr_"$CHR"
  
  for FILE in *
  do
  
    /opt/sge/bin/lx-amd64/qsub -cwd -V -b y -q long.q vcftools --gzvcf "/home/055442/reference_panels/HRCr1.1/VCF/HRC.r1-1.EGA.GRCh37.chr"$CHR".haplotypes.vcf.gz" \
                              --positions $FILE \
                              --maf 0.01 \
                              --recode \
                              --out $OUTDIR/$FILE
    echo "Job submitted"
  done
done