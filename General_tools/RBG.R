library(readxl)
library(tidyverse)
library(TwoSampleMR)
library(haven)
library(ggplot2)

setwd("~/pPRS")

T2D <- read_sav("DM_data followup study .sav") %>% rename(IID = ergoid) %>%
  mutate(IID = as.character(IID)) %>% mutate(Inci_DM_2015 = as.factor(Inci_DM_2015))

frac <- read.table("simple_frac_pheno.txt", header = T) %>%
  mutate(IID = as.character(IID))

prs <- read.table("PGS-RS1.txt", header = T, row.names = NULL) %>% 
  rename(IID = row.names) %>%
  left_join(frac, "IID") %>%
  left_join(T2D, "IID") %>% 
  arrange(wPGS,desc()) 


fit <- glm(anyfr ~ PRS_, data = prs2, family = "binomial")
summary(fit)
exp(coef(fit))
exp(confint(fit))

get_quantile <- function(x, num.quant, quant.ref){
  quant <- as.numeric(cut(x,
                          breaks = unique(quantile(
                            x, probs = seq(0, 1, 1 / num.quant)
                          )),
                          include.lowest = T))
  if(is.na(quant.ref) | is.null(quant.ref)){
    quant.ref <- ceiling(num.quant / 2)
  }
  quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
  return(quant)
}

#quant.cutoff <- "1,5,10,20,40,60,80,90,95,99,100"

set_uneven_quant <- function(quant.cutoff, ref.cutoff, num.quant, prs, quant.index){
  
  quant <- get_quantile(prs, num.quant, 1)
  quant <- factor(quant,1:num.quant)
  quant.cut <- sort(as.numeric(strsplit(quant.cutoff, split=",")[[1]]))
  if((!is.null(ref.cutoff) & sum(ref.cutoff == quant.cut)==0) | num.quant < max(quant.cut)){
    stop(
      "Invalid quant-break. quant-break must be smaller than total number of quantiles and quant-ref must be one of the quant-break"
    )
  }
  prev.name <- 0
  ref.level <- NULL
  for(i in 1:length(quant.cut)){
    up.bound <- max(which(suppressWarnings(as.numeric(levels(quant))) <= quant.cut[i]))
    cur.name <- levels(quant)[up.bound]
    name.level <- paste0("(",prev.name, ",", cur.name, "]")
    
    if(prev.name==0){
      name.level <- paste0("[",prev.name, ",", cur.name, "]")
      
    }
    if(!is.null(ref.cutoff)) {
      if (quant.cut[i] == ref.cutoff) {
        ref.level <- name.level
        
      }
    }
    range <- i:up.bound
    levels(quant)[range] <- rep(name.level, length(range))
    prev.name <- cur.name
  }
  if(is.null(ref.level)){
    writeLines("=======================================")
    writeLines("Warning: Cannot find required reference level, will use the middle level as reference")
    writeLines("=======================================")
    ref.level <- levels(quant)[ceiling(length(quant.cut)/2)]
  }
  ref.index <- which(levels(quant)==ref.level)
  quant.index <- c(ref.index,c(1:length(levels(quant)))[-ref.index])
  quant<- relevel(quant, ref=ref.level)
  return(list(quant, quant.index))
}

# Plottings ---------------------------------------------------------------

# Standard Theme for all plots
theme_sam <- NULL
theme_sam <- theme_classic()+theme(axis.title=element_text(face="bold", size=18),
                                     axis.text=element_text(size=14),
                                     legend.title=element_text(face="bold", size=18),
                                     legend.text=element_text(size=14),
                                     axis.text.x=element_text(angle=45, hjust=1))

call_quantile <-
  function(pheno.merge,
           covariance, 
           prefix,
           num_quant, 
           quant.index,
           pheno.as.quant, # is outcome continous T or F
           use.residual, 
           binary, # is outcome binary T or F
           extract, 
           device,
           uneven) {
    if (!pheno.as.quant) {
      family <- gaussian
      if (binary) {
        family <- binomial
      }
      independent.variables <-
        c("quantile", "Pheno", colnames(covariance[, !colnames(covariance) %in% c("FID", "IID")]))
      reg <-
        summary(glm(Pheno ~ ., family, data = pheno.merge[, independent.variables]))
      coef.quantiles <- (reg$coefficients[1:num_quant, 1])
      ci <- (1.96 * reg$coefficients[1:num_quant, 2])
      
      ci.quantiles.u <-
        coef.quantiles + ci
      ci.quantiles.l <-
        coef.quantiles - ci
      if (binary) {
        ci.quantiles.u <- exp(ci.quantiles.u)
        ci.quantiles.l <- exp(ci.quantiles.l)
        coef.quantiles <- exp(coef.quantiles)
      }
      
      coef.quantiles[1] <-
        ifelse(binary , 1, 0)
      ci.quantiles.u[1] <-
        ifelse(binary , 1, 0)
      ci.quantiles.l[1] <-
        ifelse(binary , 1, 0)
      quantiles.for.table <- factor(levels(pheno.merge$quantile), levels(pheno.merge$quantile))
      quantiles.df <-
        data.frame(
          Coef = coef.quantiles,
          CI.U = ci.quantiles.u,
          CI.L = ci.quantiles.l,
          DEC = quantiles.for.table
        )
      quantiles.df$Group = 0
      if (!is.null(extract)) {
        # Last element should be the cases
        quantiles.df$Group[nrow(quantiles.df)] <- 1
      }
      quantiles.df$Group <-
        factor(quantiles.df$Group, levels = c(0, 1))
      quantiles.df <- quantiles.df[order(quant.index),]
      quantiles.df$DEC <- factor(quantiles.df$DEC, levels=levels(quantiles.df$DEC)[order(quant.index)])
      row.names(quantiles.df) <- quantiles.df$DEC
      quantiles.df <- cbind(Quantile = rownames(quantiles.df), quantiles.df) 
      quant.out <- quantiles.df[,-5]
      sample.size <- as.data.frame(table(pheno.merge$quantile))
      colnames(sample.size) <- c("Quantile", "N")
      quant.out$Order <- 1:nrow(quant.out)
      quant.out <- merge(quant.out,sample.size )
      if(binary & !use.residual){
        colnames(quant.out)[2] <- "OR"
      }
      
      quant.out <-
        quant.out[order(quant.out$Order), !colnames(quant.out) %in% "Order"]
      write.table(
        quant.out,
        paste(prefix, "_QUANTILES_", Sys.Date(), ".txt", sep = ""),
        sep = "\t",
        quote = F,
        row.names = F
      )
      
      plot.quant(quantiles.df,
                   num_quant,
                   binary,
                   extract,
                   prefix,
                   uneven,
                   device)

    } else{
      # TODO: Maybe also change this to regression? Though might be problematic if we have binary pheno without cov
      pheno.sum <-
        data.frame(
          mean = numeric(num_quant),
          quantile = factor(levels(pheno.merge$quantile)[order(quant.index)],levels=levels(pheno.merge$quantile)[order(quant.index)]),
          UCI = numeric(num_quant),
          LCI = numeric(num_quant)
        )
      for (i in 1:length(levels(pheno.sum$quantile))) {
        
        cur.prs <-
          pheno.merge$PRS[as.character(pheno.merge$quantile) %in% as.character(levels(pheno.sum$quantile))[i]]
        pheno.sum$mean[i] <- mean(cur.prs, na.rm = T)
        pheno.sum$UCI[i] <-
          pheno.sum$mean[i] + sd(cur.prs, na.rm = T)
        pheno.sum$LCI[i] <-
          pheno.sum$mean[i] - sd(cur.prs, na.rm = T)
      }
      pheno.sum$Group = 0
      if (!is.null(extract)) {
        pheno.sum$Group[num_quant] = 1
      }
      pheno.sum$Group <-
        factor(pheno.sum$Group, levels = c(0, 1))
      write.table(
        pheno.sum,
        paste(prefix, "_PHENO_QUANTILES_", Sys.Date(), ".txt", sep = ""),
        sep = "\t",
        quote = F,
        row.names = F
      )

        plot.pheno.quant(pheno.sum,
                         use.residual,
                         num_quant,
                         extract,
                         prefix,
                         uneven,
                         device)
    }
    
  }

ch <- set_uneven_quant("1,5,10,20,40,60,80,90,95,99,100", 60, 100, prs$wPGS, NULL)





