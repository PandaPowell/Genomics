libraries <- c("hyprcoloc",
               "data.table",
               "tidytable",
               "tidyverse",
               "foreign")

setwd("/home/055442/ebmd_t2d_coloc/")

#browseVignettes("hyprcoloc")
invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

# pre-subsetted lead variant loci
maha_sig1 <- fread("/home/055442/hyprcoloc/one_signal_loci.txt")
maha_sig2 = fread("/home/055442/hyprcoloc/multi_signal_loci.txt")
maha_sig = bind_rows.(maha_sig1, maha_sig2)

# define +/-500kp around lead varaint
sig_regions <- maha_sig %>% mutate.(lower = pos - 500000, upper = pos + 500000,
                                    region = paste(paste(chr,lower, sep = ":"),upper, sep = "-"))
sig_regions <- as.data.table(sig_regions)
sig_regions$trait = "T2D"
# regions to extract
all_regions <- pull.(sig_regions,"region")

################### Extract SNPs summary statistics ##########################################

maha <- fread("unzip -p /home/055442/hyprcoloc/Mahajan.NatGenet2018b.T2D.European.zip", quote = " ")
# Add for loop here and replace 1 with i to extract all regions
t2d_regions <- data.frame()
for (i in 1:22){
  reg_loop <- maha[Chr == i & inrange(Pos, sig_regions$lower[sig_regions$chr == i], 
                                      sig_regions$upper[sig_regions$chr == i])]
  t2d_regions <- rbind(t2d_regions, reg_loop)
}
t2d_data <- t2d_regions %>% filter.(EAF >= 0.01 & EAF <= 0.99) %>% mutate.(trait = "T2D")
cat("Summary statistc one regions extracted")
#########################   Load GWAS summary statistic Data #########################

bmd <- fread("/home/055442/hyprcoloc/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt") %>% 
  mutate.(SNP = paste(CHR, BP, sep = ":"), trait = "eBMD") %>% 
  rename.(chr = CHR, position = BP, p = P.NI, beta = BETA, se = SE, eaf = EAF, rsid = RSID, effect_allele = EA, other_allele = NEA) %>%
  dplyr::mutate_at(c("chr", "position", "beta", "se", "p", "eaf"), as.numeric) %>%
  filter.(eaf <= 0.99 & eaf >= 0.01)

############ Load SNPS we are interested in ##################################
# lower=39535928
# upper = 40535928
# CHR=1
# trait="eBMD"
# binary=T

run_hyprcoloc <- function(lower, upper, CHR, trait, binary){
  
  # Extract ranges from GWAS data 
  
  cat("loading data for CHR", CHR, "\n")
  
  t2d_d2 <- t2d_data[Chr == CHR & inrange(Pos, lower, upper)] %>%
    select.(SNP, beta = Beta, se = SE, trait)
  
  bmd2 <- bmd[chr == CHR & inrange(position, lower,upper)] %>%
    select.(SNP, beta, se, trait)
  
  cat("SNPs in region",length(t2d_d2$SNP), "\n")
  
  # Combine datasets andremove duplicate SNPs
  
  locus <- bind_rows.(t2d_d2,bmd2) %>%
    mutate.(unique = paste(SNP,trait, sep = ";")) %>%
    distinct.(unique, .keep_all = T) %>% select.(-unique)
  
  # Create Matrix of beta values
  
  b <-  locus %>% select.(-se) %>%
    pivot_wider.(names_from = "trait", values_from = "beta", values_fill = NA) %>% 
    drop_na.() %>% tibble::column_to_rownames("SNP")
  
  b2 <- locus %>% select.(-se) %>%
    pivot_wider.(names_from = "trait", values_from = "beta", values_fill = NA) %>%
    filter.(is.na(trait) == F)
  
  # Check missing SNPs in relation to T2D
  nas <- sapply(b2, function(x) sum(is.na(x)))
  cat("GWAS with the most missing snps in comparsion to ",trait,print(sort(nas)), "\n")
  
  b_matrix <- as.matrix(b)
  cat("Number of overlapping SNPs in region",dim(b)[1], "\n")
  
  # Create Matrix of se values
  se <- locus %>% select.(-beta) %>%
    pivot_wider.(names_from = "trait", values_from = "se", values_fill = NA) %>% 
    drop_na.() %>% tibble::column_to_rownames("SNP")
  se_matrix <- as.matrix(se)
  
  # Assign Column and Row names
  trait_name <- colnames(b)
  snp_list <- rownames(b)
  
  # code binary traits
  if(binary == T){
    binary.traits = c(rep(0,grep(pattern = trait, trait_name)-1),1,rep(0,length(trait_name)-grep(pattern = trait, trait_name)))
  } else {
    binary.traits = rep(0,length(trait_name))
  }

  # Run regional hyprcoloc with prob >0.8
  results <- hyprcoloc(b_matrix, se_matrix, trait.names=trait_name, snp.id=snp_list, binary.outcomes = binary.traits, prior.c = 0.02,
                       reg.thresh = 0.8, bb.selection = "reg.only")
  pos_prob <- as.data.frame(results$results)
  pos_prob <- pos_prob[is.na(pos_prob$posterior_prob) == F,]
  
  if(sum(grepl(trait, pos_prob$traits, fixed = T)) >= 1) {
    cat(trait," region colocalised, writing SNPs\n")
    # Output overlapping SNPs as a dataframe
    overlapping_snps <- data.frame(stringr::str_split_fixed(t2d_d2$SNP, ":", n = 2))
    range <- paste(lower, upper, sep = "-")
    
    # Create df with range and traits
    df2 <- data.frame(region = paste(CHR, range, sep = ":"), traits = pos_prob$traits)
    print(df2)
    
    fwrite(df2, paste("res", paste(range, ".txt", sep = ""), sep = "_"), quote = F, row.names = F, col.names = F)
    
    dir.create(paste("./chr", CHR, sep = "_"), showWarnings = F)
    fwrite(overlapping_snps, paste( paste("chr_", CHR, sep = ""), range,sep = "/"), quote = F, row.names = F, col.names = F, sep = "\t")
  } else {
    cat("No traits colocalise with T2D in this region \n")
  }
  
  return(pos_prob)
}

cat("Performing colocalization \n")

for (CHR in sort(unique(sig_regions$chr)) ){ 
  
  # It is necessary to re-create this dataframe of regions everytime 
  sig_regions2 <- maha_sig %>%
    mutate.(lower = pos - 500000, upper = pos + 500000,region = paste(paste(chr,lower, sep = ":"),upper, sep = "-")) %>%
    filter.(chr == CHR)
  
  yoo <- lapply(1:nrow(sig_regions2), function(x) 
    run_hyprcoloc(sig_regions2$lower[x], sig_regions2$upper[x], CHR, c("T2D"), F)
  )
  
  if (is.null(yoo[[1]])){
    stop("Something went wrong with hyprcoloc")
  }
  
  cat("Hyprcoloc complete, writing results \n")
  
  capture.output(yoo, file=paste(CHR,"txt", sep = "."))
  
  trait = unique(sig_regions$trait)
  
  loj <- c()
  trait_names <- c()
  for (i in c(1:length(yoo))){
    if (any( grep(trait,yoo[[i]], fixed = T) )){
      roj <- yoo[[i]]$candidate_snp[grep(trait, yoo[[i]]$traits)]
      loj <- c(loj,roj)
      temp_names <- yoo[[i]]$traits[grep(trait, yoo[[i]]$traits)]
      temp_names <- stringr::str_split(temp_names, pattern = ", ")
      trait_names <- c(trait_names, temp_names)
    }
  }
  
  cat("Extracting candidate snps from colocalised traits \n")
  
  coloc_ass <- data.frame()
  for (i in 1:length(loj)){
    bmd_can_snps <- bmd %>% filter.(SNP %in% loj[i], trait %in% unlist(trait_names[i])) %>%
      select.(rsid, SNP, effect_allele, other_allele, eaf,beta, se, p, trait)
    risk_associations <- bmd_can_snps
    coloc_ass <- bind_rows.(coloc_ass, risk_associations)
  }
  
  cat("writing table \n")
  
  t2d_can_snps <- t2d_data %>% filter.(SNP %in% loj) %>% left_join.(coloc_ass, by = "SNP")
  write.table(t2d_can_snps, file=paste(paste("coloc_ass.chr", CHR, sep = ""),"txt", sep = "."), quote = F, row.names = F, sep = ",")
}


