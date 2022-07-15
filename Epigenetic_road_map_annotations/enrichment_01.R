#---------------------------------------------------#
#            Enrichment - analysis                  # 
#---------------------------------------------------#

#------------------------------------------------------------#
#  Libraries and dependencies                                #
#------------------------------------------------------------#
library (gplots)
library (RColorBrewer)
library(ComplexHeatmap)

#------------------------------------------------------------#
#                                 #
#------------------------------------------------------------#

counts <- read.table("All_counts_FINEMAP_Oct_2019.txt", header = TRUE, row.names = 1, sep = "\t")
enrichment <- counts[1:12]

for(i in 1:12){
  for(j in 1:127){
    prop.result <- prop.test(x = c(counts[i,j], counts[13,i]), n = c(counts[i,128], 533), alternative = "greater")
    enrichment[i,j] <-prop.result$p.value
  }
}

write.table(enrichment, "Enrichment_results_FINEMAP.txt", row.names =  TRUE, quote = FALSE, sep = "\t")


enrichment <- read.table("Enrichment_results_FINEMAP.txt", header=T)
# multiple testing correction: 3.55e-04
# select all significant total results
tresults <- t(enrichment)
sig.result <- tresults[tresults[,12]<= 0.00000005,]
sig.result <- tresults

# annotation_data 
ann <- read.table("EID_identifier_ROADMAP.txt", header=T, sep= "\t", stringsAsFactors = T)
res <- as.data.frame(sig.result)
res$EID <- row.names(res)

res.ann <- merge(ann, res, by ="EID", all.x = FALSE, all.y = FALSE)
res.ann.order <- res.ann[order(res.ann$`12.Total`), ]
row.names(res.ann.order) <- res.ann.order$Epigenome

# plot results 
results.plot =  t(-log10(res.ann.order[,7:18]))
names.plot <- t(res.ann.order)

#------------------------------------------------------------#
#  Customizing and plotting the heatmap                      #
#------------------------------------------------------------#
library(circlize)
my_palette = colorRamp2(c(0, 5, 7.5, 15 , 30), c("White","White", "#df65b0", "#dd1c77", "#980043"))

mat <- results.plot
ann <- as.vector(res.ann.order$GROUP)
col_fun <- list(bar = c("Adipose" = "#AF5B39", "Blood & T-cell" = '#55A354', "Brain" = '#C5912B', "Digestive" = '#C58DAA',
                        "ENCODE2012" = '#000000', "Epithelial" = '#FF9D0C', "ES-deriv" = '#4178AE', "ESC" = '#924965',
                        "Heart" = '#D56F80', "HSC & B-cell" = '#678C69', "IMR90" = '#678C69', "iPSC" = '#69608A', 
                        "Mesench" = '#B65C73', "Muscle" = '#C2655D',"Myosat" = '#E67326', "Neurosph" = '#FFD924', 
                        "Other" = '#999999',"Sm. Muscle" = '#F182BC',"Thymus" = '#DAB92E'))

ha = HeatmapAnnotation(bar = ann, col = col_fun)
ht = Heatmap(mat,
        col = my_palette,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_dend_side = "right", 
        column_names_side = "top",
        bottom_annotation = ha
          )
draw(ht)









