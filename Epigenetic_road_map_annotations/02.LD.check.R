#----------------------------------------------#
#           check LD Files                     #
#----------------------------------------------#
### Dependancies
library("data.table")
setDTthreads(threads = 1)

### Args
args <- commandArgs(trailingOnly=TRUE)
list.snps    <- args[1]

### Script
files.snp  <- list.files(path = list.snps, pattern="LD.*", full.names=TRUE)
rsmid      <- fread("zcat ~/GWAS/GWAS_Files/HRCv1.1_snplist_MAF0.001_noInDel.txt.gz", header= TRUE, stringsAsFactors = FALSE, data.table = FALSE)

### Loop-di-check ###
for (i in 1:length(files.snp)){
  tmp <- read.table(files.snp[i], header=T,strings=F)
  if((is.na(tmp[1,1]))==TRUE){
    x <-files.snp[i]
    y <- sapply(strsplit((sapply(strsplit(x, "LD."), "[[",3)),".txt"),"[[",1)
    z <- rsmid[(which(rsmid$SNP == y)),]
    z$CHR <- paste0("chr",z$CHR)
    write.table(z, x, quote = FALSE, row.names = FALSE, sep = "\t")
  } else {
    print(paste0(files.snp[i],"is good, no changes made"))
  }
}

#########
q()
