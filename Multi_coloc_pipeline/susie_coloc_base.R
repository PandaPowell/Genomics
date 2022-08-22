libraries <- c("coloc",
               "susieR",
               "data.table",
               "tidytable",
               "tidyverse",
               "foreign",
               "purrr")

invisible(suppressMessages(lapply(libraries, require, character.only = TRUE)))

setDTthreads(10)

args = commandArgs(trailingOnly = TRUE)

####### Define the following variables #####
wd = "/home/055442/ebmd_t2d_coloc/" # Working directory
lead_loci = "/home/055442/ebmd_t2d_coloc/maha_lead_loci.txt" # path of base GWAS lead loci
base_trait = "T2D" # base GWAS trait
base_gwas = "/home/055442/hyprcoloc/Mahajan.NatGenet2018b.T2D.European.zip"
hypr_res = "trait_regions.txt"
############################################

setwd(wd)

# Load primary lead variant loci and regions
cat("Loading primary regions\n")
maha_sig <- fread(file = lead_loci, header = TRUE, stringsAsFactors = FALSE)

# define +/-500kp around lead varaint
sig_regions <- maha_sig %>% mutate.(lower = pos - 500000, 
                                    upper = pos + 500000,
                                    region = paste(paste(chr,lower, sep = ":"),upper, sep = "-"))

sig_regions <- as.data.table(sig_regions)
sig_regions$trait = base_trait

# regions to extract
all_regions <- pull.(sig_regions,"region")

################### Load Primary GWAS summary statistics ##########################################
cat("Loading primary summary statistcs\n")

maha <- fread("unzip -p /home/055442/hyprcoloc/Mahajan.NatGenet2018b.T2D.European.zip", quote = " ")
# Add for loop here and replace 1 with i to extract all regions
t2d_regions <- data.frame()
for (i in 1:22){
  reg_loop <- maha[Chr == i & inrange(Pos, sig_regions$lower[sig_regions$chr == i], 
                                      sig_regions$upper[sig_regions$chr == i])]
  t2d_regions <- rbind(t2d_regions, reg_loop)
}
t2d_data <- t2d_regions %>% filter.(EAF >= 0.01 & EAF <= 0.99) %>% mutate.(trait = base_trait)
cat("Primary summary statistc regions extracted\n")
  
# load dataframe of the traits that acheive regional colocalisation > 0.8 at each regions
tr <- fread(hypr_res, fill = T, header = F)

tr$chr <- as.integer(unlist(lapply(1:nrow(tr), function(x) strsplit(tr$V1[x], split = ":")[[1]][1] )))
tr$lower <- as.integer(unlist(lapply(1:nrow(tr), function(x) strsplit( strsplit(tr$V1[x], split = "-")[[1]][1], ":")[[1]][2] )))
tr$upper <- as.integer(unlist(lapply(1:nrow(tr), function(x) strsplit(tr$V1[x], split = "-")[[1]][2] )))
tr$range <- unlist(lapply(1:nrow(tr), function(x) strsplit(tr$V1[x], split = ":")[[1]][2] ))

tr <- tr %>% rename.(chr_range = V1) %>%
  tidyr::unite(c_traits, V2:V3, sep = ",", na.rm = T) %>%
  mutate.(lower = if_else(lower < 0 | is.na(lower) == T,0,lower))
tr$upper[tr$chr_range == "2:-77856-922144"] = 922144
tr$range[tr$chr_range == "2:-77856-922144"] = "0-922144"
  #%>% filter(grepl("T2D",c_traits) == T)
head(tr)

## We need something here that fixes the regions so that we can't have any negative regions but it will be better to put this in a pre-step

## Sample size of GWAS needs to be included in the summary statisics otherwise it has to be input manually

### Start of SUSIE COLOC function ###
run_coloc_susie = function(x){

  cat("Loading LD and BIM files\n")
  
  # Load LD and bim file data
  LDfilename <- paste("./LD_regions/", tr$range[x], ".recode/", tr$range[x],sep = "") #
  BIMfilename <- paste("./LD_regions/", tr$range[x], ".recode/", tr$range[x], sep = "") #
  LD <- fread(paste(LDfilename, "ld", sep = "." ))
  BIM <- fread( paste(BIMfilename, "bim", sep = ".") )
  colnames(BIM) = c("chr", "rsid", "dk", "pos", "alt", "ref")
  BIM$SNP <- paste(BIM$chr, BIM$pos, sep = ":")
  
  # Remove duplicates from LD matrix
  colnames(LD) <- BIM$SNP
  LD$SNP <- BIM$SNP
  LD <- LD[!duplicated(LD$SNP),]
  LD <- LD[, !duplicated(colnames(LD)), with=F ]
  
  # Remove duplicates from BIM file
  BIM  = distinct.(BIM, "SNP", .keep_all = T)
  
  # define CHR & Traits
  CHR = tr$chr[x] # 
  coloc_traits <- stringr::str_split(c(tr[x,2]), ",") #
  coloc_traits = coloc_traits[[1]]
  coloc_traits = coloc_traits[coloc_traits != ""]
  
  coloc_tr_snps = data.frame()
  
  # filter T2D data
  d2 <- t2d_data[SNP %in% LD$SNP] %>% left_join.(BIM, "SNP") %>%
    mutate.(EA.x = ifelse.(EA == ref & NEA == alt, EA, NEA),
            NEA.x = ifelse.(EA == ref & NEA == alt, NEA, EA),
            EAF.x = ifelse.(EA == ref & NEA == alt, EAF, 1-EAF),
            Beta.x = ifelse(EA == ref & NEA == alt, Beta,Beta*-1)) %>%
    filter.(EA.x == ref & NEA.x == alt) %>%
    select.(SNP, chr = Chr, position = Pos, EA = EA.x, NEA = NEA.x, EAF = EAF.x, Beta = Beta.x, se = SE, p =Pvalue, n = Neff) %>%
    distinct.(SNP, .keep_all = T) %>%
    mutate.(varbeta = se^2, z = Beta/se) %>%
    arrange.(SNP)
  
  # Loop through traits that coloclaise
  
    for (name in c("eBMD")){ # coloc_traits[coloc_traits != base_trait] , A hacking way to run this quick is just enter the trait names here
    
    cat("Colocalising region",tr$range[x], "and trait", name, "GWAS \n")
    
    # Load colocalised data, allign with BIM
    d <- fread(paste(name, ".txt", sep = "")) %>% 
      mutate.(SNP = paste(chr,position, sep=":"), p = ifelse.(p == 0, as.numeric(2e-308), as.numeric(p))) %>%
      filter.(SNP %in% d2$SNP) %>% left_join.(BIM, "SNP") %>%
      mutate.(EA.x = ifelse.(effect_allele == ref & other_allele == alt, effect_allele, other_allele),
              NEA.x = ifelse.(effect_allele == ref & other_allele == alt, other_allele, effect_allele),
              EAF.x = ifelse.(effect_allele == ref & other_allele == alt, eaf, 1-eaf),
              Beta.x = ifelse(effect_allele == ref & other_allele == alt, beta,beta*-1)) %>%
      filter.(EA.x == ref & NEA.x == alt) %>%
      select.(SNP, chr, position,effect_allele = EA.x, other_allele = NEA.x, eaf = EAF.x, beta = Beta.x, se, p,n)
    
    d1 = d %>% distinct.(SNP, .keep_all = T) %>%
      mutate.(varbeta = se^2) %>%
      rename.(N = n, MAF = eaf) %>%
      arrange.(SNP) %>% mutate.(z = beta/se)
    
    # Format colocalised trait 
    b = c(d1$beta)
    names(b) = d1$SNP
    vb = c(d1$varbeta)
    names(vb) = d1$SNP
    maf = c(d1$MAF)
    names(maf) = c(d1$SNP)
    
    # filter out bad SNPs from LD matrix
    LD = LD %>% filter.(SNP %in% d1$SNP)
    LD <- LD[, c(d1$SNP,"SNP"), with=F ]
    
    # remove SNP column from LD dataframe
    LD_coloc = as.matrix( LD[, 1:(ncol(LD)-1)] )
    dimnames(LD_coloc)[[1]] <- dimnames(LD_coloc)[[2]]
    
    
    D1 <- list(d1$SNP, d1$position, b, vb, maf, d1$N[1], "quant", LD_coloc)
    names(D1) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type", "LD")
    check_dataset(D1, req= "LD", warn.minp = 1e-05)
    check_alignment(D1)
    plot_dataset(D1, main = name)
    
    # filter out bad snps from T2DM data
    d2 = d2 %>% filter.(SNP %in% d1$SNP)
    # Format T2D data
    b2 = c(d2$Beta)
    names(b2) = d2$SNP
    vb2 = c(d2$varbeta)
    names(vb2) = d2$SNP
    maf2 = c(d2$EAF)
    names(maf2) = c(d2$SNP)
    
    D2 = list(d2$SNP, d2$position, b2, vb2, maf2, d2$n[1], "cc", LD_coloc)
    names(D2) = c("snp", "position", "beta", "varbeta", "MAF", "N", "type", "LD")
    check_dataset(D2, req= "LD",warn.minp = 1e-05)
    check_alignment(D2)
    plot_dataset(D2, main = "T2D")
    
    # Set coverage low so we capture as many SNPs as possible, these can be filtered out later
    cat("running susie for T2D \n")
    S2 = try(runsusie(D2, coverage = 0.1))
    print(summary(S2))
    
    cat("running SUSIE for", name, "\n")
    S1 = try(runsusie(D1, coverage = 0.1))
    print(summary(S1))
    
    if(!inherits(S2, "try-error") | !inherits(S2, "try-error") ){
      test_res = coloc.susie(S1,S2)
    } else{
      test_res=NA
    }
    
    
    if(is.na(test_res) == F){
      tr2 = test_res$summary[PP.H4.abf >= 0.6]
    } else {tr2 = data.table()}
    
    # If susie cant identify credible sets or doesnt detect coloclaisation then run coloc
    if(is.na(test_res) |  nrow(tr2) == 0){
      
      cat("no snps identified, running COLOC \n")
      res = coloc.abf(D1,D2)
      
      
        coloc_snps = subset(res$results,SNP.PP.H4 == max(SNP.PP.H4))
        colnames(coloc_snps)[1] = "SNP...1"
        coloc_snps$PP.H4.abf = res$summary[6]
        df = filter.(d1, SNP %in% coloc_snps$SNP...1) %>%
          mutate.(trait = name)
        
        df2 = filter.(d2, SNP %in% coloc_snps$SNP...1) %>%
          rename.(effect_allele = EA, other_allele = NEA , MAF = EAF, beta = Beta, N = n) %>% 
          mutate.(trait = "T2D")
        
        EX_SNP = bind_cols.(df,df2, .name_repair = "unique") %>%
          left_join.(coloc_snps, "SNP...1") %>%
          mutate.(region = tr$chr_range[x])
        
        coloc_tr_snps = bind_rows.(EX_SNP, coloc_tr_snps)
        
        print(coloc_tr_snps)

      
    } else {
  
    cat("running colocalisation \n")
    res = coloc.susie(S1,S2)
    
    # Check to make sure we have results to work with

    colnames(res$summary)[2] = "SNP...1"
    print(res$summary)
    print(res$summary[PP.H4.abf > 0.6])
    
    #sensitivity(res,"H4 > 0.9",row=3,dataset1=D1,dataset2=D2)
    
    
   coloc_snps =  res$summary
   
   df = filter.(d1, SNP %in% coloc_snps$SNP...1) %>%
     mutate.(trait = name)
   
   df2 = filter.(d2, SNP %in% coloc_snps$hit2) %>%
     rename.(effect_allele = EA, other_allele = NEA , MAF = EAF, beta = Beta, N = n) %>% 
     mutate.(trait = "T2D")
   
   EX_SNP = bind_cols.(df,df2, .name_repair = "unique") %>%
     left_join.(coloc_snps, "SNP...1") %>%
     mutate.(region = tr$chr_range[x])
   
   coloc_tr_snps = bind_rows.(EX_SNP, coloc_tr_snps) }
    
  }
  return(coloc_tr_snps)
}

### In case of an error return NA ### 
run_coloc_susie2 = purrr::possibly(run_coloc_susie, otherwise = NA)

results = data.frame() # 2,38,43
for (i in c(1:nrow(tr))) {
  temp_res = run_coloc_susie2(i)
  results = bind_rows.(results,temp_res)
}

cat("Writing final results")

dir.create("./susie_coloc_results", showWarnings = F)

fwrite(results, "./susie_coloc_results/eBMD_new_results.txt", quote = F, row.names = F)
 
