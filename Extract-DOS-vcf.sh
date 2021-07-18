#!/bin/bash

#########
# $1 - File with the 81 SNPs
########

INPUT=$1
OUTDIR="/home/055442/Preeclampsia/OUTPUT_GENR4"
GENDIR="/data/GENR4/Imputed/1000G_PhaseIIIv5/GenR_Moms"

mkdir $OUTDIR

echo "Extracting SNPs from .vcf files"

for i in {1..22}
do
    awk -v chr=$i '{FS=","; if ($1==chr) print $1"\t"$2}' $INPUT > $OUTDIR/chr${i}.snps
    /opt/sge/bin/lx-amd64/qsub -cwd -V -b y -q long.q vcftools --gzvcf $GENDIR/GenRMoms_IDC_Chr${i}.dose.vcf.gz --positions $OUTDIR/chr${i}.snps --recode --out $OUTDIR/chr${i}
done

echo "Waiting for the extraction to finish"

#### MAKE SURE NO INTERACTION QUEUE IS OPEN
while true
do
    a=$(qstat -u 055442 | wc -l)
    if [[ $a -eq 0 ]]
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

head -n 1 $OUTDIR/chr1.dosages.txt > ./TBBMD-DOSAGES-GENR4.txt
tail -n +2 -q $OUTDIR/chr*.dosages.txt >> ./TBBMD-DOSAGES-GENR4.txt
sed -i 's/\([0-9]*\)\t$/\1/g' ./TBBMD-DOSAGES-GENR4.txt

#### RELEVANT ONLY IF EXTRACTING GENOTYPE DATA, NOT DOSAGE DATA!
#echo "Changing genotype data to 0, 1, 2"
#sed -i 's/\(0|0\)/0/g; s/\(0|1\)\|\(1|0\)/1/g; s/\(1|1\)/2/g' ./TBBMD-DOSAGES-GENR4.txt

MIS=$(comm -13 <(cut -f2,3 ./TBBMD-DOSAGES-GENR4.txt | tail -n +2 | tr "\t" ":" | sort) <(cut -f1,2 -d "," $INPUT | tail -n +2 | tr "," ":" | sort))

echo "Missing SNPs: ${MIS}"
echo $MIS > ./TBBMD-MISSING-SNPs-GENR4.txt

#### DO IT INDEPENDENTLY
#echo "Running the script to extract PGS"
echo "Run the script to extract PGS"

#Rscript Create-PGS.R ./TBBMD-DOSAGES-GENR4.txt $INPUT
