library(data.table)
library(tidytable)
library(foreign)

args = commandArgs(trailingOnly = TRUE)

# Argument 1 - Pathway to directory
# Argument 2 - lead SNP and region file path


if (length(args) == 0) {
  stop("Supply PATH file as well as SNPs file")
} else {
  DIR <- args[1]
  REG <- args[2]
}

setwd(DIR)

cat("Reading in SNP/REGIONS files\n")

# Load primary lead variant loci and regions
maha_sig <- fread(file = REG, header = TRUE, stringsAsFactors = FALSE)

# define +/-500kp around lead varaint
sig_regions <- maha_sig %>% mutate.(lower = pos - 500000, upper = pos + 500000,
                                    region = paste(paste(chr,lower, sep = ":"),upper, sep = "-"))
sig_regions <- as.data.table(sig_regions)
sig_regions$trait = "T2D"

# regions to extract
all_regions <- pull.(sig_regions,"region")

cat("Reading in GWAS file\n")

gwas <- fread("/home/055442/hyprcoloc/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz") %>%
  select.(SNP = SNPID,
         rsid = SNP,
         chr = CHR,
         position = BP,
         effect_allele = ALLELE1,
         other_allele = ALLELE0,
         eaf = A1FREQ,
         beta = logOR,
         se = logOR.SE,
         p = P.NI,
         n = N)

cat("Extracting regions from GWAS file\n")

regions <- data.frame()
for (i in 1:22){
  print(paste("Extracting SNPs in regions in CHR  ",i))
  reg_loop <- gwas[chr == i & inrange(position, sig_regions$lower[sig_regions$chr == i], 
                                             sig_regions$upper[sig_regions$chr == i])]
  regions <- rbind(regions, reg_loop)
}

head(regions)

regions$trait <- "Fracture"

cat("Writing file\n")

fwrite(regions, file = "./Fracture.txt", quote = F, row.names = F, sep = "\t")
