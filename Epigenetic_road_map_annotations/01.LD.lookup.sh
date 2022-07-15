#!/bin/bash
# Get LD for SNPs r >=0.8, 1000G based
# get token from https://ldlink.nci.nih.gov/?tab=apiaccess
# CGB, V1.0
# $1 = input file 

while read -r line
do
	call="https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=${line}&pop=EUR&r2_d=r2&token=d8fcf4e19b43"
	curl -k -X GET $call | awk '$7 >= 0.8 {print $1,$2}' | sed -e "s/:/ /g" | sed -e "s/ /\t/g" > test.$line.txt
	echo $'SNP\tCHR\tPOS' > LD.$line.txt
	tail -n +2 -q test.$line.txt >> LD.$line.txt
	rm test.$line.txt
done < $1
exit
