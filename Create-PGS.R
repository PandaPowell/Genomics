library(data.table)
library(foreign)
library(ggplot2)

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
SNPS <- fread(file = SNP, header = TRUE, sep = ",", stringsAsFactors = FALSE)
SNPS[, SNP := paste(CHR, BP, sep = ":")]

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

cat("Checking for duplicates in DOSAGE file\n")

if (length(DOSAGES) != length(SNPS)) {
  DOSAGES <- unique(DOSAGES, by = "SNP")
  cat("duplicates removed working with", dim(DOSAGES)[1], "SNPs\n")
}

cat("Matching DOSAGES with SNP/GWAS effect alleles\n")

if (!identical(SNPS$SNP, DOSAGES$SNP)) {
  stop("ERROR, DOSAGE and SNP/GWAS data not aligned")
}

DROPPED <- c()

for (i in 1:nrow(SNPS)) {
  if (SNPS[i, "A1"] != DOSAGES[i, "ALT"]) {
    if (SNPS[i, "A1"] == DOSAGES[i, "REF"]) {
      DOSAGES[i, 6:ncol(DOSAGES)] <- -1 * (DOSAGES[i, 6:ncol(DOSAGES)] - 2)
      DOSAGES[i, "REF"] <- DOSAGES[i, "ALT"]
      DOSAGES[i, "ALT"] <- SNPS[i, "A1"]
    } else {
      cat("Problem with phasing in SNP", DOSAGES[i, SNP], "; Dropping it.\n")
      DROPPED <- c(DROPPED, i)
    }
  }
}

cat("Making sure DOSAGES match SNP/GWAS efect allales.\n")

if (!identical(DOSAGES$ALT, SNPS$A1)) {
  stop("ERROR, DOSAGE and SNP/GWAS data not harmonized")
} else {
  cat("All seems OK!")
}

cat("Changing data to align with effect increasing alleles.\n")

for (i in 1:nrow(SNPS)) {
  if (SNPS[i, beta] < 0) {
    SNPS[i, beta := -1 * beta]
    SNPS[i, EAF := 1 - EAF]
    DOSAGES[i, ALT := REF]
    DOSAGES[i, ]$REF <- SNPS[i, ]$A1
    SNPS[i, ]$A1 <- DOSAGES[i, ]$ALT
    DOSAGES[i, 6:ncol(DOSAGES)] <- -1 * (DOSAGES[i, 6:ncol(DOSAGES)] - 2)
  }
}

if (any(SNPS$beta < 0)) {
  stop("ERROR, some effect sizes still negative.\n")
} else {
  cat("All effect sizes seem OK!\n")
}

if (!identical(DOSAGES$ALT, SNPS$A1)) {
  stop("ERROR, DOSAGE and SNP/GWAS data not harmonized")
} else {
  cat("All alleles seem harmonized, extracting PGS!\n")
}

PRS <- DOSAGES[, 6:ncol(DOSAGES)]
WPRS <- DOSAGES[, 6:ncol(DOSAGES)]
WPRS <- as.data.frame(as.matrix(PRS) * SNPS$beta)

final_df <- data.frame(PGS = colSums(PRS), wPGS = colSums(WPRS))

cat("Writing PGS into file \"PGS-GENR4.txt\".\n")
write.table(x = final_df, file = "~/PGS-GENR4.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
