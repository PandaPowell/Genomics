library(data.table)
library(foreign)
library(ggplot2)

# FORMAT only names that matter, Chr, Pos, EA, EAF, Beta, 

setDTthreads(10)
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Supply DOSAGES file as well as SNPs file")
} else {
  DOS <- args[1]
  SNP <- args[2]
}

cat("Reading in DOSAGE and SNP/GWAS files\n")
DOSAGES <- fread(file = DOS, header = TRUE, stringsAsFactors = FALSE)
DOSAGES[, SNP := paste(CHR, POS, sep = ":")]
SNPS <- fread(file = SNP, header = TRUE, stringsAsFactors = FALSE)
SNPS[, SNP := paste(Chr, Pos, sep = ":")]

# Remove duplicates
SNPS = SNPS[!duplicated(SNPS$SNP)]
DOSAGES = DOSAGES[!duplicated(DOSAGES$SNP)]

cat("Harmonizing files\n")
cat("#################\n")
SNPS <- SNPS[order(SNPS$SNP)]
DOSAGES <- DOSAGES[order(DOSAGES$SNP)]

MIS <- SNPS[!which(SNPS$SNP %in% DOSAGES$SNP)]$SNP
cat("SNPs missing in DOSAGE data:", MIS, "\n")

SNPS <- SNPS[which(SNPS$SNP %in% DOSAGES$SNP)]
DOSAGES <- DOSAGES[which(DOSAGES$SNP %in% SNPS$SNP)]

cat("Working with", dim(DOSAGES)[1], "SNPs\n")

# cat("Changing SNP/GWAS data to reflect alleles with an increasing effect")
# COMPL <- list("A", "C", "G", "T")
# names(COMPL) <- c("T", "G", "C", "A")
# SNPS[beta < 0, A1 := as.vector(unlist(COMPL[SNPS[beta < 0, A1]]))]
# SNPS[beta < 0, beta := beta * (-1)]

cat("Matching DOSAGES with SNP/GWAS effect alleles\n")

if (!identical(SNPS$SNP, DOSAGES$SNP)) {
  stop("ERROR, DOSAGE and SNP/GWAS data not aligned")
}

DROPPED <- c()

# ALT allele should match effect allele
for (i in 1:nrow(SNPS)) {
  if (SNPS[i, "EA"] != DOSAGES[i, "ALT"]) {
    if (SNPS[i, "EA"] == DOSAGES[i, "REF"]) {
      DOSAGES[i, 6:ncol(DOSAGES)] <- -1 * (DOSAGES[i, 6:ncol(DOSAGES)] - 2)
      DOSAGES[i, "REF"] <- DOSAGES[i, "ALT"]
      DOSAGES[i, "ALT"] <- SNPS[i, "EA"]
    } else {
      cat("Problem with phasing in SNP", DOSAGES[i, SNP], "; Dropping it.\n")
      DROPPED <- c(DROPPED, i)
    }
  }
}

if (!is.null(DROPPED)){
  SNPS = SNPS[-DROPPED, ]
  DOSAGES = DOSAGES[-DROPPED,]
}

cat("Making sure DOSAGES match SNP/GWAS efect allales.\n")

if (!identical(DOSAGES$ALT, SNPS$EA)) {
  stop("ERROR, DOSAGE and SNP/GWAS data not harmonized")
} else {
  cat("All seems OK!")
}

cat("Changing data to align with effect increasing alleles.\n")

for (i in 1:nrow(SNPS)) {
  if (SNPS[i, Beta] < 0) {
    SNPS[i, Beta := -1 * Beta]
    SNPS[i, EAF := 1 - EAF]
    DOSAGES[i, ALT := REF]
    DOSAGES[i, ]$REF <- SNPS[i, ]$EA
    SNPS[i, ]$EA <- DOSAGES[i, ]$ALT
    DOSAGES[i, 6:ncol(DOSAGES)] <- -1 * (DOSAGES[i, 6:ncol(DOSAGES)] - 2)
  }
}

if (any(SNPS$Beta < 0)) {
  stop("ERROR, some effect sizes still negative.\n")
} else {
  cat("All effect sizes seem OK!\n")
}

if (!identical(DOSAGES$ALT, SNPS$EA)) {
  stop("ERROR, DOSAGE and SNP/GWAS data not harmonized")
} else {
  cat("All alleles seem harmonized, extracting PGS!\n")
}

PRS <- DOSAGES[, 6:ncol(DOSAGES)]
WPRS <- DOSAGES[, 6:ncol(DOSAGES)]
WPRS <- as.data.frame(as.matrix(PRS) * SNPS$Beta)

final_df <- data.frame(PGS = colSums(PRS), wPGS = colSums(WPRS))

study = strsplit(tail(strsplit(DOS, "/")[[1]],1),"-")[[1]][2]
trait = tail(strsplit(sapply(strsplit(SNP, "_SNPs"), "[",1), "/")[[1]], n=1)

cat("Writing PGS into file \n")
fwrite(x = final_df, file = paste("./PGS-",trait,"-",study,".txt", sep = ""), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

