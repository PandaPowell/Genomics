#----------------------------------------------------------------#
#                     Tissue specificity                         #
#                 GO GWAS data [indepent SNV]                    #
#----------------------------------------------------------------#
# Version 1.0
#----------------------------------------------------------------#
#                        Dependancies                            #
#----------------------------------------------------------------#
library(ggplot2)

#----------------------------------------------------------------#
#                        Data.preparation                        #
#----------------------------------------------------------------#
enrichment <- read.table("Enrichment_results_FINEMAP.txt", header=T)
sig.results <- t(enrichment)

# annotation_data 
ann                <- read.table("EID_identifier_ROADMAP.txt", header=T, sep= "\t", stringsAsFactors = T)
res                <- as.data.frame(sig.result)
res$EID            <- row.names(res)
res.ann            <- merge(ann, res, by ="EID", all.x = FALSE, all.y = FALSE)
res.ann$Total      <- -log10(res.ann$`12.Total`)
res.ann$fig        <- "lightgrey"
res.ann$fig        <- ifelse(res.ann$Total < 7.9, "#999999",res.ann$fig)
res.ann$fig        <- ifelse(res.ann$EID == "E026"|res.ann$EID == "E049" | res.ann$EID == "E025"|
                               res.ann$EID == "E023"| res.ann$EID == "E129", "#D55E00", res.ann$fig)
res.ann$fig        <- ifelse(res.ann$EID == "E052"|res.ann$EID == "E108" | res.ann$EID == "E107"|
                               res.ann$EID == "E089"| res.ann$EID == "E090" | res.ann$EID == "E083" |
                               res.ann$EID == "E095", "#0072B2", res.ann$fig)
res.ann$fig        <- ifelse(res.ann$EID == "E097"|res.ann$EID == "E080" | 
                               res.ann$EID == "E091", "#56B4E9", res.ann$fig)
res.ann$fig        <- ifelse(res.ann$EID == "E053", "#000000", res.ann$fig)

# order on 
res.ann.order       <- res.ann[order(res.ann$Total,decreasing = TRUE), ]
res.ann.order$order <- seq.int(nrow(res.ann.order))

#----------------------------------------------------------------#
#                        Data.preparation                        #
#----------------------------------------------------------------#
df <- res.ann.order
l <- as.vector(df$Mnemonic)
g <- as.vector(df$fig)

ggplot(data=df, aes(x = order, y = Total, group = GROUP)) +
  geom_bar(stat="identity",  color="black", fill = g) +
  scale_x_continuous(breaks=0:126, labels = l) +
  geom_hline(yintercept = 7.9, size = 1, linetype = "dashed") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) 




