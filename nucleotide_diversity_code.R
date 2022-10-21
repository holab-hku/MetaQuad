library(vcfR)

args <- commandArgs(trailingOnly = TRUE)


vcf <- read.vcfR(args[1], verbose = FALSE )

#vcf <- read.vcfR("cellSNP.cells_ARG.vcf.gz", verbose = FALSE )

data <- as.data.frame(vcf@gt)

row.names(data) <- paste(vcf@fix[,1],vcf@fix[,2],vcf@fix[,4],vcf@fix[,5],sep = "_")
data <- data[,-1]



passed_mutation <- read.csv(args[2])
#passed_mutation <- read.csv("BIC_params_ARG.csv")

passed_mutation <- passed_mutation[which(passed_mutation$PASS_cutoff == TRUE & passed_mutation$PASS_MINSAMPLES == TRUE),]

data <- data[which(row.names(data) %in% passed_mutation$variant_name),]


cellSNP_informative_AD <- data


for (i in 1:ncol(cellSNP_informative_AD)) {
  cellSNP_informative_AD[,i] <- read.table(text = as.character(cellSNP_informative_AD[,i]), sep = ":")$V2
}

cellSNP_informative_AD[cellSNP_informative_AD=="."] <- 0

cellSNP_informative_AD <- as.data.frame(apply(cellSNP_informative_AD, 2, as.numeric))


cellSNP_informative_DP <- data

for (i in 1:ncol(cellSNP_informative_DP)) {
  cellSNP_informative_DP[,i] <- read.table(text = as.character(cellSNP_informative_DP[,i]), sep = ":")$V3
}


cellSNP_informative_DP <- as.data.frame(apply(cellSNP_informative_DP, 2, as.numeric))


row.names(cellSNP_informative_DP) <- row.names(data)
row.names(cellSNP_informative_AD) <- row.names(data)


cellSNP_informative_AP <- cellSNP_informative_AD/cellSNP_informative_DP



cellSNP_informative_nd <- cellSNP_informative_AP*2-2*cellSNP_informative_AP*cellSNP_informative_AP

cellSNP_informative_nd$genome <- gsub("_[0-9]*_._.$","",row.names(cellSNP_informative_nd))

cellSNP_informative_nd[is.na(cellSNP_informative_nd)] <- 0
cellSNP_informative_nd <- aggregate(.~genome,cellSNP_informative_nd,sum)

row.names(cellSNP_informative_nd) <- cellSNP_informative_nd$genome
cellSNP_informative_nd <- cellSNP_informative_nd[,-1]
colnames(cellSNP_informative_nd) <- gsub("\\..*$","",colnames(cellSNP_informative_nd))


#####normalization



genome_length <- read.csv(args[3],row.names = 1)
#genome_length <- read.csv("genome_length_ARG.csv",row.names = 1)

colnames(genome_length) <- gsub('_coverage.tab','', colnames(genome_length))
colnames(genome_length) <- gsub('\\.','-', colnames(genome_length))

cellSNP_informative_nd_1 <- cellSNP_informative_nd

for (gene in row.names(cellSNP_informative_nd_1)) {
  for (sample in colnames(cellSNP_informative_nd_1)) {
    cellSNP_informative_nd_1[gene,sample] <- cellSNP_informative_nd_1[gene,sample]/genome_length[gene,sample]
  }
}



write.csv(cellSNP_informative_nd_1,"Nucleotide_diversity.csv")









