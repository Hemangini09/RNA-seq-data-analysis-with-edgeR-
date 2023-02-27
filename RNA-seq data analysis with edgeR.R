# install required packages
BiocManager::install("edgeR")

# import important libraries
library(edgeR)
library(dplyr)

# read count matrix file
exprs <- get(load("dataFilt.RData"))

# Load clinical data 
clinical_data <- read.csv("Clinical.csv")
clinical_data <- data.frame(clinical_data, row.names = T)

all(rownames(clinical_data) %in% colnames(exprs))
all(rownames(clinical_data) == colnames(exprs))


# build a count matrix
count_matrix <- as.matrix(exprs)
sample_info <- clinical_data$braf_status

# group samples
group = factor(sample_info)
group = relevel(group, ref="WT")

# Create DGE object
dim(dataFilt)
dge <- DGEList(counts = count_matrix, group = group)
dim(dataFilt)


# Calculate normalization factor
dge <- calcNormFactors(dge)

# Filter lowly- expressed genes
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge <- dge[-drop,] 
head(dge$counts)


# Estimate dispersion
dge <- estimateDisp(dge)

# Testing for DEG
test <- exactTest(dge)          # If it is not working try test  <- exactTest(dge,pair=c(1,2))
top_degs <- topTags(test,n="Inf")
top_degs <- top_degs$table

# Total number of differentially expressed genes at FDR < 0.05 and threshold value 1
summary(decideTests(test, lfc = 1))


# save degs file
write.csv(top_degs,file=" top_degs.csv")
save(top_degs, file="top_degs.RData")

