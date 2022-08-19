#!/bin/bash

#########
# $1 - File with the 81 SNPs
# 1st column should be CHR, 2nd POS
# Necessary column names CHR,POS,A1,A2,beta
########

INPUT=$1
INPUTNAME="${INPUT%.*}" 
RSDIR="/home/055442/pPRS/RS1/"
OUTDIR="/home/055442/pPRS/RS1/"${INPUTNAME}
GENDIR="/home/055442/RS1/Imputed/HRC_release1.1"

echo $OUTDIR

mkdir $RSDIR
mkdir $OUTDIR
mkdir ./DOSAGES

echo "Extracting SNPs from .vcf files"

for i in {1..22}
do
    awk -v chr=$i '{if ($1==chr) print $1":"$2"-"$2}' $INPUT > $OUTDIR/chr${i}.snps
    cat $OUTDIR/chr${i}.snps | tr "\n" " " | sed 's/ $/\n/g' > $OUTDIR/snps
    mv $OUTDIR/snps $OUTDIR/chr${i}.snps

    if [ -s $OUTDIR/chr${i}.snps ]
    then
        tabix -h $GENDIR/chr${i}.dose.vcf.gz $(cat $OUTDIR/chr${i}.snps) > $OUTDIR/chr${i}.recode.vcf
    else
        rm $OUTDIR/chr${i}.snps
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

head -n 1 $OUTDIR/chr1.dosages.txt > ./${INPUTNAME}-RS1-DOSAGES.txt
tail -n +2 -q $OUTDIR/chr*.dosages.txt >> ./${INPUTNAME}-RS1-DOSAGES.txt
sed -i 's/\([0-9]*\)\t$/\1/g' ./${INPUTNAME}-RS1-DOSAGES.txt

#### RELEVANT ONLY IF EXTRACTING GENOTYPE DATA, NOT DOSAGE DATA!
#echo "Changing genotype data to 0, 1, 2"
#sed -i 's/\(0|0\)/0/g; s/\(0|1\)\|\(1|0\)/1/g; s/\(1|1\)/2/g' ./LIPO-DG-DOSAGES-GENR4.txt

MIS=$(comm -13 <(cut -f2,3 ./${INPUTNAME}-RS1-DOSAGES.txt | tail -n +2 | tr "\t" ":" | sort) <(cut -f1,2 -d "," $INPUT | tail -n +2 | tr "\t" ":" | sort))

echo "Missing SNPs: ${MIS}"
echo $MIS > ./MISSING-SNPs-RS1-${INPUTNAME}.txt

#### DO IT INDEPENDENTLY
#echo "Running the script to extract PGS"
echo "Run the script to extract PGS"

#Rscript Create-PGS.R ./TBBMD-DOSAGES-GENR4.txt $INPUT
