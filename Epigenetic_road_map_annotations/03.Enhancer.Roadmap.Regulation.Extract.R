#------------------------------------------------------------#
# ENCODE/Roadmap chromatin states extraction script          #
#------------------------------------------------------------#

# --Version 3.0 :  October 2019

### Pre-script Run ###
# SNP data, contains header: SNP, CHR, BP, (rest header/colums does not matter)
# mnemonics data, NO HEADER, colums: CHR, BP_START, BP_ENd, chromain state

#------------------------------------------------------------#
# Arguments and load files                                   #
#------------------------------------------------------------#
args <- commandArgs(trailingOnly=TRUE)

data.15     <- args[1] # Path to folfer with CHROMM_15.txt.gz
list.SNPs   <- args[2] # PatH to folder with SNPS
out.name    <- args[3] # output name.prefix
ncores      <- as.numeric(args[4]) # number of cores allowed to run

### DO NOT CHANGE CODE BELOW, WILL AFFECT Script function  ###
#------------------------------------------------------------#
#                 Dependancies                               #
#------------------------------------------------------------#
require(data.table)
setDTthreads(threads = ncores)

require(parallel)
require(doParallel)
registerDoParallel(cores = ncores)

#------------------------------------------------------------#
#                 Safety Check                               #
#------------------------------------------------------------#
# safety check if more cores are used than specified
# Script will be killed

safety <- as.numeric(getDoParWorkers())
if(safety < ncores){
  print("ERROR: more cores used than requested")
  q()
} else{
  print(paste0("will use ", ncores, " cores for analysis."))
}

#------------------------------------------------------------#
#  Functions                                                 #
#------------------------------------------------------------#

### ChromHMM mining Function parallel ##

StateFinder <- function(roadmap, snp){
  foreach(j = 1:nrow(snp), .combine = rbind) %dopar% {
    roadmap[snp[j,2] == roadmap[,1] & (roadmap[,2] < snp[j,3] & roadmap[,3] >= snp[j ,3]),]
  }
}
  
#---------------------------------------------------------#
# ENCODE/Roadmap chromatin states (loop-di-loop)             #                                    #
#------------------------------------------------------------#

### SNP lists detection ###
snp <- fread(list.SNPs, header = TRUE, stringsAsFactors = F, data.table = FALSE)

### CHROMM HMM detection ###
files             <- paste0("zcat ", (list.files(path = data.15, pattern="*mnemonics.bed", full.names=TRUE)))
naming            <- lapply(strsplit((list.files(path = data.15, full.names = FALSE)), split = "_"), function(x){x[1]})
result.lst        <- vector("list", length(files))
names(result.lst) <- naming

### Loop-di-Loop  parallel !!! ###
for(i in 1:length(files)) {
  roadmap.file  <- fread(cmd = files[i], header = FALSE, 
                        stringsAsFactors = FALSE, data.table = FALSE)
  chrom.states  <- StateFinder(roadmap.file, snp)
  chromm.table  <- as.data.frame(table(chrom.states[,4]))
  if(is.null(chromm.table)){chromm.table <- as.data.frame(matrix(0, nrow = 2, ncol=2, dimnames = list(c("1", "2"), c("Var1", "Freq"))))}
  Enhan.count   <- chromm.table[chromm.table[,1] == "6_EnhG" | chromm.table[,1] == "7_Enh" | chromm.table[,1] == "12_EnhBiv" , ]
  output.loop   <- as.data.frame(matrix((sum(Enhan.count$Freq)), nrow=1, ncol=1, dimnames = list ("1", naming[i])))
  
  # save out.put in named list
  result.lst[[i]]  <- output.loop
}

### Creating data Set ###
states <- as.data.frame(rbind(unlist(result.lst)))

#------------------------------------------------------------#
# saving and naming output file                              #
#------------------------------------------------------------#

out.file <- paste0(out.name, "_15State_enhancer_counts.txt")
write.table(states, out.file, quote = FALSE, row.names = FALSE)

#------------------------------------------------------------#
# Exit Script                                                #                                    #
#------------------------------------------------------------#
warnings()
q()
