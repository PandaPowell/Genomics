
#######################################################################
### Script for mapping cloest genes to a particular genetic variant ###
#######################################################################

# SNP file format required #
# Three columns tab seperated, chromsome, SNP start genomic position, SNP end genomic position #
# example, chr1    110222901       110222901

SNPS=$1
fname="${SNPS}"
nname="${fname%%.*}"

mkdir "./closest_gene_"$nname

cd "./closest_gene_"$nname

cat $SNPS > $nname"_snps.bed"

bedtools sort -i $nname"_snps.bed" > $nname"_sorted_snps.bed"

(echo -e "Chr\SNP_Start\SNP_End\Chr\tx_Start\tx_End\Gene_name\score\strand\cd_Start\cd_End"; # this line here just creates headers for the new file
  bedtools closest -a $nname"_sorted_snps.bed" -b ../NCBI37.3.gene.loc.sorted.bed) > $nname"annotated.tsv"
  
echo "Mapping complete"