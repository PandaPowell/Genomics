library(LDlinkR)
library(data.table)
library(foreign)
library(tidyverse)

####################################################
#### Script for finding LD proxies using LDlink ####
####################################################
# Created by Sam Ghatan
# Last updated 28th November 2022
# Requires snps with rsnumbers & sumstats

setDTthreads(10)

args = commandArgs(trailingOnly = T)

if ( length(args) == 0 ){
  stop("Supply list of variants as rsid IDs & summary statistics")
} else {
  SNP = args[1]
  SUM <- args[2]
}

cat("Reading in SUM STATS and SNP files\n")
SNPS <- fread(file = SNP, header = F, stringsAsFactors = FALSE)

if(sum(grepl("rs", SNPS$V1) == 0)){
   stop("Supply list of variants as rsid IDs & summary statistics")
}

SS = fread(file = SUM, header = T, stringsAsFactors = FALSE)
SS$SNP = paste(SS$chromosome,SS$base_pair_location, sep=":")

cat(nrow(SNPS),"SNPs that require LD proxies \n")

cat(nrow(SNPS),"Finding LD proxies \n")
LDproxy_batch(SNPS$V1, pop = "CEU", r2d = "r2", token = "d8fcf4e19b43", append = T)

proxy = fread("combined_query_snp_list_grch37.txt", fill=F)
proxy$SNP = gsub("chr","",proxy$Coord)

proxy2 = proxy[,c(1,3:9,13)] %>% filter(R2 > 0.60) %>% filter(SNP %in% SS$SNP) %>%
         distinct(RS_Number, .keep_all=T)

n = which(proxy2$RS_Number %in% SNPS$V1)
proxy2$query_snp = ""

for(i in 1:length(n)){
  if(i==length(n)){
  proxy2$query_snp[n[i]:nrow(proxy2)] = proxy2$RS_Number[n[i]]
  } else {
  proxy2$query_snp[n[i]:n[i+1]] = proxy2$RS_Number[n[i]]
  }
}

proxy3 = proxy2 %>% filter(!Alleles %in% c("(A/T)", "(C/G)", "(G/C)", "(T/A)")) %>%
         arrange(desc(R2)) %>% distinct(query_snp, .keep_all=T)

ex_proxies = SS %>% filter(SNP %in% proxy3$SNP) %>%
  left_join(proxy3, "SNP")
  
cat("Writing extracting proxy SNPs\n")

fwrite(ex_proxies, "proxies.txt", row.names = F, quote = F, sep=",")
