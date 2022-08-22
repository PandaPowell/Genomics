library(data.table)
library(tidytable)
library(tidyverse)

setwd("/home/055442/ebmd_t2d_coloc/susie_coloc_results/")

comb = fread("/home/055442/ebmd_t2d_coloc/susie_coloc_results/eBMD_new_results.txt") %>% filter(PP.H4.abf > 0.59)
comb_old = fread("/home/055442/ebmd_t2d_coloc/susie_coloc_results/eBMD_results_combined.txt")

maha = fread("/home/055442/multi_coloc/Mahajan_leadSNPs.csv") %>% filter.(region %in% comb$region)

cc  = comb %>% select(SNP...1,gene,trait...13) %>%
  group_by(gene) %>% mutate(coloc_traits = list(trait...13)) %>% distinct(gene, .keep_all = T) %>%
  select(-trait...13)

maha2 = fread("/home/055442/multi_coloc/Mahajan_leadSNPs.csv") %>% filter.(`Nearest gene` %in% cc$gene) %>%
  rename(gene = `Nearest gene`) %>% left_join(cc, "gene")
colnames(maha2)[3] = "position"

tr = maha2 %>%
  mutate(position = gsub(",","",position),SNP = paste(Chromosome, position,sep=":"))

# Load T2D data
t2d_data <- fread("/home/055442/multi_coloc/ALL_data/All.Multi.t2d_regions.txt") %>%
  mutate.(trait = "T2D")

# filter T2D data
d2 <- t2d_data[SNP %in% tr$SNP] %>%
  select.(SNP, chr = Chr, position = Pos, EA , NEA , EAF , Beta , se = SE, p =Pvalue, n = Neff) %>%
  distinct.(SNP, .keep_all = T) %>%
  mutate.(varbeta = se^2, z = Beta/se) %>%
  arrange.(SNP)

gly_traits = c("T2D")

# DF of traits and there respective GWAS sample sizes
traits = c("ALT", "GGT", "HDL_cholesterol", "LDL_cholesterol", "triglycerides", "VAT", "AFR_males", "AFR_females",
                         "LFR_males", "LFR_females", "TFR_males", "TFR_females", "BMI", "WHR", "Proinsulin_levels", "Homeostasis_model_assessment_of_insulin_resistance",
                         "Modified_Stumvoll_Insulin_Sensitivity_Index", "Leptin", "Adiponectin","HbA1c", "Fasting_glucose", "Fasting_insulin","2hGlu")



#coloc_traits <- unlist(str_split(unlist((tr[x,'coloc_traits'])), ",")) #
#coloc_traits = coloc_traits[coloc_traits != "T2D"]
#coloc_traits = coloc_traits[coloc_traits != ""]

results = data.frame()
for (name in traits){
  cat("extracting",name,"/n")
  
  # filter T2D snps that 
  coloc_snps = tr[grepl(name ,tr$coloc_traits)]
  
  # Load colocalised data, allign with BIM
  d <- fread(paste("/home/055442/multi_coloc/GWAS_data/", name, ".txt", sep = "")) %>% #
    mutate.(SNP = paste(chr,position, sep=":"), p = ifelse.(p == 0, 2e-308, p)) %>%
    filter(SNP %in% coloc_snps$SNP) %>%
    select.(SNP, chr, position,effect_allele , other_allele , eaf , beta , se, p)
  
  d1 = d %>% distinct.(SNP, .keep_all = T) %>%
    mutate.(varbeta = se^2) %>%
    rename( MAF = eaf) %>%
    arrange.(SNP) %>% mutate(z = beta/se) %>% mutate(trait =name)
  
  results = bind_rows.(results,d1)
  
}

res2 = d2 %>% left_join(results, "SNP") %>% left_join.(tr[,c(1,16,17,20)], by = "SNP")

fwrite(res2, "extracted_mult_coloclaised_snps.txt", quote = F, row.names = F)

  
