#!/bin/bash

# Make folders for each region

cd $1

for FILE in *.vcf
do
mkdir "${FILE%.*}" && mv "$FILE" "${FILE%.*}"
done

# Calculate LD for each region

for d in */
do
  cd $1"/$d"
  
  echo $d
  
  for file in *.recode.vcf
  do
  
  fname="${file}"
  nname="${fname%%.*}"
  
  echo "$nname"
  echo "$fname"
  
  /opt/tools/bin/plink19 --noweb \
        --vcf "$fname" \
        --vcf-idspace-to _ \
        --const-fid \
        --allow-extra-chr 0 \
        --split-x b37 no-fail \
        --make-bed \
        --keep-allele-order \
        --out $nname

  echo "Calculating LD for SNPs in region"
  
  /opt/tools/bin/plink19 --bfile  $nname \
          --r square \
          --keep-allele-order \
          --out $nname
  done

done