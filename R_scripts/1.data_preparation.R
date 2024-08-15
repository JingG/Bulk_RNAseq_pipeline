## script to take in the count matrix and sample table into DESeq2 object.

## output ##
# 1. dds.RDS
# 2. TPM with annotation 


options(stringsAsFactors=FALSE)
options(digits=10)
library(knitr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# source("/Users/jing/nBox//Projects/bin/CommonFunctionR/plot_pdf_png.r")

### file path ###
project_dir <- "/path/to/project/"
path_input <- paste0(project_dir,"inputFiles/")
path_output <- paste0(project_dir,"Rspace/results/")
path_deposit <- paste0(project_dir,"Rspace/dataDeposit/")

### read in data ###

sampleInfo <- read.csv(file = paste0(path_input, "sampleInfo.csv"))
rownames(sampleInfo) <- sampleInfo$Enrolment.ID
df <- read.table(file = paste0(path_input, "combined_counts.txt"), header = T, row.names = 1)
sampleInfo$Classification <- as.factor(sampleInfo$Classification)
sampleInfo$Diagnosis <- as.factor(sampleInfo$Diagnosis)

##### filter genes 
## filte ron 

check_columns_greater_than_zero <- function(matrix) {
  apply(matrix, 1, function(x) sum(x > 0) > 10)
}

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = sampleInfo,
                              design = ~ Classification + Diagnosis)

dds <- dds[ check_columns_greater_than_zero(counts(dds)), ] #more than 10 samples has expression 

dds <- estimateSizeFactors(dds) 
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

saveRDS(dds, file = paste0(path_deposit, "inital_dds.RDS"))

############################# Annotation ##################################

ensembl.genes <- readRDS("/Users/jing/nBox/Projects/bin/CommonFunctionR/ensembl_genes_gtf110.RDS")

ensembl.genes <- as.data.frame(ensembl.genes)

############################### TPM #######################################
library(GenomicFeatures)

count.df <- as.data.frame(dds@assays@data$counts)

hs.gtf.db <- makeTxDbFromGFF( "~/nBox/Projects/bin/reference_genome/Homo_sapiens.GRCh38.110.chr.gtf", format="gtf" )
ensembl.genes = genes(hs.gtf.db)
exonsByGene <- GenomicFeatures::exonsBy(hs.gtf.db, by="gene" )
exonsByGene <- reduce(exonsByGene)
hs.gene.length = sum(width(exonsByGene))
saveRDS(hs.gene.length, file = "~/nBox/Projects/bin/reference_genome/Homo_sapiens.GRCh38.110.geneLength_by_ENSID.RDS")

# the gene shall be no change if using the same gtf version 
count.df <- count.df[rownames(count.df) %in% names(hs.gene.length),]

hs.gene.length2 <- hs.gene.length[rownames(count.df)]
identical(names(hs.gene.length2), rownames(count.df))

############# calculate TPM #############
genes_rpk = count.df/((hs.gene.length2[rownames(count.df)]) / 1000)
#genes_rpk[is.na(genes_rpk)] <- 0  #for corrected annotated with concordant gtf, there shall be no NA -- use this to check the gtf version 
PM_samples_scalingfactor = colSums(genes_rpk) / 1e6
TPM = sweep(genes_rpk, 2, PM_samples_scalingfactor, "/")
# write.xlsx(x = TPM, file = paste0(path_output, "TPM.xlsx"), col.names = T, row.names = T, append = F, showNA = T)
TPM.anno <- cbind(ensembl.genes[rownames(TPM),], TPM)

saveRDS(TPM.anno, file = paste0(path_deposit, "TPM.RDS"))





