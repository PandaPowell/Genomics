#!/bin/bash

#########
# $1 - File with the 81 SNPs
# 1st column should be CHR, 2nd POS
# Necessary column names CHR,POS,A1,A2,beta, comma seperated
########

INPUT=$1 # SNP file needs to be in current directory, i.e cant have a directory name before the file name
INPUTNAME="${INPUT%.*}" 
RSDIR="/home/055442/pPRS/GENR4/"
OUTDIR=$RSDIR"/"${INPUTNAME}
GENDIR="/home/055442/GENR4/Imputed/1000G_PhaseIIIv5/GenR_Kids/"
OUT="GENR4"

mkdir $RSDIR
mkdir $OUTDIR

echo "Extracting SNPs from .vcf files"

for i in {1..22}
do
    awk -v chr=$i '{if ($1==chr) print $1"\t"$2}' $INPUT > $OUTDIR/chr${i}.snps
    /opt/sge/bin/lx-amd64/qsub -cwd -V -b y -q long.q vcftools --gzvcf $GENDIR/GenRKids_IDC_Chr${i}.dose.vcf.gz --positions $OUTDIR/chr${i}.snps --recode --out $OUTDIR/chr${i}
done

echo sleep 5

# Get last job ID
JOBID=$(qstat | tail -n 1 | cut -d ' ' -f 2)

echo "Waiting for the extraction to finish, last jobID ="$JOBID

while true
do
   if ! qstat | cut -d ' ' -f 2 | grep -q $JOBID
   then
       break
   fi
done


echo "Extraction done"
echo "Filling in the information"
echo "##########################"

for i in {1..22}
do
    echo "Filling in the information for chromosome ${i}"
    bcftools query -l $OUTDIR/chr${i}.recode.vcf > $OUTDIR/chr${i}.personid
    printf "SNP"$'\t'"CHR"$'\t'"POS"$'\t'"REF"$'\t'"ALT"$'\t' > $OUTDIR/chr${i}.dosages.txt
    tr "\n" "\t" < $OUTDIR/chr${i}.personid >> $OUTDIR/chr${i}.dosages.txt
    printf "\n" >> $OUTDIR/chr${i}.dosages.txt
    bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' $OUTDIR/chr${i}.recode.vcf >> $OUTDIR/chr${i}.dosages.txt
done

echo "Dosages extracted for the ALT allele."

head -n 1 $OUTDIR/chr1.dosages.txt > $OUTDIR/${INPUTNAME}-$OUT-DOSAGES.txt
tail -n +2 -q $OUTDIR/chr*.dosages.txt >> $OUTDIR/${INPUTNAME}-$OUT-DOSAGES.txt
sed -i 's/\([0-9]*\)\t$/\1/g' $OUTDIR/${INPUTNAME}-$OUT-DOSAGES.txt

#### RELEVANT ONLY IF EXTRACTING GENOTYPE DATA, NOT DOSAGE DATA!
#echo "Changing genotype data to 0, 1, 2"
#sed -i 's/\(0|0\)/0/g; s/\(0|1\)\|\(1|0\)/1/g; s/\(1|1\)/2/g' ./LIPO-DG-DOSAGES-GENR4.txt

MIS=$(comm -13 <(cut -f2,3 $OUTDIR/${INPUTNAME}-$OUT-DOSAGES.txt | tail -n +2 | tr "\t" ":" | sort) <(cut -f1,2 -d " " $INPUT | tail -n +2 | tr " " ":" | sort))

echo "Missing SNPs: ${MIS}"
echo $MIS > $OUTDIR/MISSING-SNPs-${INPUTNAME}.txt

# Remove individual chromosome files to make it neat
for i in {1..22}
do
  rm $OUTDIR/chr${i}.*
done
  

#### DO IT INDEPENDENTLY
#echo "Running the script to extract PGS"
echo "Run the script to extract PGS"

#Rscript Create-PGS.R ./TBBMD-DOSAGES-GENR4.txt $INPUT
