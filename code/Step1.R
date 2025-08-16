## Code related to "Genome-wide changes in lncRNA, alternative splicing, and cortical patterning in autism"
## This script takes count data and filters for expressed genes using the criteria described in the Supplementary Information.
## It then uses conditional quantile normalization to perform gene length and GC normalization across samples and adjusts the data.
## Note that, for ease of sharing the code, the data are normalized without additional samples in the set (resulting in differences) at the GC normalization step and bootstrapping/permutation steps are not performed. Results are compared with those reported in the manuscript to demonstrate there are only minor differences and the final results are robust to these modifications on the data.

## Set current R session to be in the appropriate working directory

##### Set libraries and defaults and load provided files
rm(list=ls())
options(stringsAsFactors=FALSE)
file_path<-"C:/Users/ariad/OneDrive/Desktop/Proyect WCQN/data/provided/ASD_RNAseq_ExpressionData.Rdata"
#"C:/Users/ariad/OneDrive/Desktop/Transcriptomics/Genome-wide-changes-in-lncRNA-alternative-splicing-and-cortical-patterning-in-autism-master/data/provided/ASD_RNAseq_ExpressionData.Rdata"
#"C:/Users/palom/Desktop/Transcriptomics/Genome-wide-changes-in-lncRNA-alternative-splicing-and-cortical-patterning-in-autism-master/data/provided/ASD_RNAseq_ExpressionData.Rdata"
load(file=file_path)

## Check the contents
summary(datMeta)

#superior temporal gyrus (STG, also known as Brodmann’s area (BA) 41/42), prefrontal cortex (BA9) and cerebellar vermis
unique(datMeta$RegionID)

#cheecking same BrainID samples
sample_counts <- table(datMeta$BrainID)
summary(sample_counts)

###The genes to exclude have been selected as those that satisfied at least one of the following conditions:

#Split by sections before normalizing
keepC <- datMeta[,"RegionID"] == "vermis"
keepF <- datMeta[,"RegionID"] == "ba9"
keepT <- datMeta[,"RegionID"] == "ba41-42-22"


# 1. Cufflinks filtering: lower bound > 0 in 80% of all samples
datExpr <- datExpr.Cufflinks[,keepC]
passvec <- apply(datLB.gene[,keepC]>0,1,sum)>(0.8*sum(keepC)) ## OK in 80% of samples
datExpr.Cufflinks.C <- datExpr[passvec,]


datMeta.Cufflinks.C <- datMeta[keepC, ]
datMeta.Cufflinks.C <- datMeta.Cufflinks.C[order(datMeta.Cufflinks.C$RIN, decreasing = TRUE), ]
datMeta.Cufflinks.C <- datMeta.Cufflinks.C[!duplicated(datMeta.Cufflinks.C$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.Cufflinks.C), rownames(datMeta.Cufflinks.C))
datExpr.Cufflinks.C <- datExpr.Cufflinks.C[, commonSamples]
datMeta.Cufflinks.C <- datMeta.Cufflinks.C[commonSamples, ]
dim(datMeta.Cufflinks.C)  # Check dimensions to ensure correct subsetting
dim(datExpr.Cufflinks.C)  # Check dimensions to ensure correct subsetting


## FRONTAL CUFFLINKS
datExpr <- datExpr.Cufflinks[,keepF]
passvec <- apply(datLB.gene[,keepF]>0,1,sum)>(0.8*sum(keepF)) ## OK in 80% of samples
datExpr.Cufflinks.F <- datExpr[passvec,]

# Filter metadata
datMeta.Cufflinks.F <- datMeta[keepF, ]
datMeta.Cufflinks.F <- datMeta.Cufflinks.F[order(datMeta.Cufflinks.F$RIN, decreasing = TRUE), ]
datMeta.Cufflinks.F <- datMeta.Cufflinks.F[!duplicated(datMeta.Cufflinks.F$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.Cufflinks.F), rownames(datMeta.Cufflinks.F))
datExpr.Cufflinks.F <- datExpr.Cufflinks.F[, commonSamples]
datMeta.Cufflinks.F <- datMeta.Cufflinks.F[commonSamples, ]
dim(datMeta.Cufflinks.F)  # Check dimensions to ensure correct subsetting
dim(datExpr.Cufflinks.F)  # Check dimensions to ensure correct subsetting

## TEMPORAL CUFFLINKS
datExpr <- datExpr.Cufflinks[,keepT]
passvec <- apply(datLB.gene[,keepT]>0,1,sum)>(0.8*sum(keepT)) ## OK in 80% of samples
datExpr.Cufflinks.T <- datExpr[passvec,]

# Filter metadata
datMeta.Cufflinks.T <- datMeta[keepT, ]
datMeta.Cufflinks.T <- datMeta.Cufflinks.T[order(datMeta.Cufflinks.T$RIN, decreasing = TRUE), ]
datMeta.Cufflinks.T <- datMeta.Cufflinks.T[!duplicated(datMeta.Cufflinks.T$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.Cufflinks.T), rownames(datMeta.Cufflinks.T))
datExpr.Cufflinks.T <- datExpr.Cufflinks.T[, commonSamples]
datMeta.Cufflinks.T <- datMeta.Cufflinks.T[commonSamples, ]
dim(datMeta.Cufflinks.T)  # Check dimensions to ensure correct subsetting
dim(datExpr.Cufflinks.T)  # Check dimensions to ensure correct subsetting

pass.cufflinks.genes <- union(rownames(datExpr.Cufflinks.C), union(rownames(datExpr.Cufflinks.F), rownames(datExpr.Cufflinks.T)))    

# 2. HTSeq union exon filtering: counts > 10 in 80% of all samples
datExpr <- datExpr.HTSC.unionexon[,keepC]
passvec <- apply(datExpr,1,quantile,0.8) > 10 ## OK in 80% of samples
datExpr.HTSC.unionexon.C <- datExpr[passvec,]

# Filter metadata
datMeta.HTSC.unionexon.C <- datMeta[keepC, ]
datMeta.HTSC.unionexon.C <- datMeta.HTSC.unionexon.C[order(datMeta.HTSC.unionexon.C$RIN, decreasing = TRUE), ]
datMeta.HTSC.unionexon.C <- datMeta.HTSC.unionexon.C[!duplicated(datMeta.HTSC.unionexon.C$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.HTSC.unionexon.C), rownames(datMeta.HTSC.unionexon.C))
datExpr.HTSC.unionexon.C <- datExpr.HTSC.unionexon.C[, commonSamples]
datMeta.HTSC.unionexon.C <- datMeta.HTSC.unionexon.C[commonSamples, ]
dim(datMeta.HTSC.unionexon.C)  # Check dimensions to ensure correct subsetting
dim(datExpr.HTSC.unionexon.C)  # Check dimensions to ensure correct subsetting

## FRONTAL UNIONEXON
datExpr <- datExpr.HTSC.unionexon[,keepF]
passvec <- apply(datExpr,1,quantile,0.8) > 10 ## OK in 80% of samples
datExpr.HTSC.unionexon.F <- datExpr[passvec,]

# Filter metadata
datMeta.HTSC.unionexon.F <- datMeta[keepF, ]
datMeta.HTSC.unionexon.F <- datMeta.HTSC.unionexon.F[order(datMeta.HTSC.unionexon.F$RIN, decreasing = TRUE), ]
datMeta.HTSC.unionexon.F <- datMeta.HTSC.unionexon.F[!duplicated(datMeta.HTSC.unionexon.F$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.HTSC.unionexon.F), rownames(datMeta.HTSC.unionexon.F))
datExpr.HTSC.unionexon.F <- datExpr.HTSC.unionexon.F[, commonSamples]
datMeta.HTSC.unionexon.F <- datMeta.HTSC.unionexon.F[commonSamples, ]
dim(datMeta.HTSC.unionexon.F)  # Check dimensions to ensure correct subsetting
dim(datExpr.HTSC.unionexon.F)  # Check dimensions to ensure correct subsetting

## TEMPORAL UNIONEXON
datExpr <- datExpr.HTSC.unionexon[,keepT]
passvec <- apply(datExpr,1,quantile,0.8) > 10
datExpr.HTSC.unionexon.T<- datExpr[passvec,]

# Filter metadata
datMeta.HTSC.unionexon.T <- datMeta[keepT, ]
datMeta.HTSC.unionexon.T <- datMeta.HTSC.unionexon.T[order(datMeta.HTSC.unionexon.T$RIN, decreasing = TRUE), ]  
datMeta.HTSC.unionexon.T <- datMeta.HTSC.unionexon.T[!duplicated(datMeta.HTSC.unionexon.T$BrainID), ] # Keep only the first occurrence of each BrainID  

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.HTSC.unionexon.T), rownames(datMeta.HTSC.unionexon.T))
datExpr.HTSC.unionexon.T <- datExpr.HTSC.unionexon.T[, commonSamples]
datMeta.HTSC.unionexon.T <- datMeta.HTSC.unionexon.T[commonSamples, ]
dim(datMeta.HTSC.unionexon.T)  # Check dimensions to ensure correct subsetting
dim(datExpr.HTSC.unionexon.T)  # Check dimensions to ensure correct subsetting

pass.unionexon.genes <- union(rownames(datExpr.HTSC.unionexon.C), union(rownames(datExpr.HTSC.unionexon.F), rownames(datExpr.HTSC.unionexon.T)))

# 3. HTSeq whole gene filtering: counts > 10 in 80% of all samples
# (Remember to clean the suffix for whole gene names first)
rownames(datExpr.HTSC.wholegene) <- substr(rownames(datExpr.HTSC.wholegene), 1, 15)

datExpr <- datExpr.HTSC.wholegene[,keepC]
passvec <- apply(datExpr,1,quantile,0.8) > 10 ## OK in 80% of samples
datExpr.HTSC.wholegene.C <- datExpr[passvec,]

# Filter metadata
datMeta.HTSC.wholegene.C <- datMeta[keepC, ]
datMeta.HTSC.wholegene.C <- datMeta.HTSC.wholegene.C[order(datMeta.HTSC.wholegene.C$RIN, decreasing = TRUE), ]  
datMeta.HTSC.wholegene.C <- datMeta.HTSC.wholegene.C[!duplicated(datMeta.HTSC.wholegene.C$BrainID), ] # Keep only the first occurrence of each BrainID
# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.HTSC.wholegene.C), rownames(datMeta.HTSC.wholegene.C))
datExpr.HTSC.wholegene.C <- datExpr.HTSC.wholegene.C[, commonSamples]
datMeta.HTSC.wholegene.C <- datMeta.HTSC.wholegene.C[commonSamples, ]
dim(datMeta.HTSC.wholegene.C)  # Check dimensions to ensure correct subsetting
dim(datExpr.HTSC.wholegene.C)  # Check dimensions to ensure correct subsetting

## FRONTAL WHOLEGENE
datExpr <- datExpr.HTSC.wholegene[,keepF]
passvec <- apply(datExpr,1,quantile,0.8) > 10 ## OK in 80% of samples
datExpr.HTSC.wholegene.F <- datExpr[passvec,]

# Filter metadata
datMeta.HTSC.wholegene.F <- datMeta[keepF, ]
datMeta.HTSC.wholegene.F <- datMeta.HTSC.wholegene.F[order(datMeta.HTSC.wholegene.F$RIN, decreasing = TRUE), ]
datMeta.HTSC.wholegene.F <- datMeta.HTSC.wholegene.F[!duplicated(datMeta.HTSC.wholegene.F$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.HTSC.wholegene.F), rownames(datMeta.HTSC.wholegene.F))
datExpr.HTSC.wholegene.F <- datExpr.HTSC.wholegene.F[, commonSamples]
datMeta.HTSC.wholegene.F <- datMeta.HTSC.wholegene.F[commonSamples, ]
dim(datMeta.HTSC.wholegene.F)  # Check dimensions to ensure correct subsetting
dim(datExpr.HTSC.wholegene.F)  # Check dimensions to ensure correct subsetting


datExpr <- datExpr.HTSC.wholegene[,keepT]
passvec <- apply(datExpr,1,quantile,0.8) > 10
datExpr.HTSC.wholegene.T <- datExpr[passvec,]

# Filter metadata
datMeta.HTSC.wholegene.T <- datMeta[keepT, ]
datMeta.HTSC.wholegene.T <- datMeta.HTSC.wholegene.T[order(datMeta.HTSC.wholegene.T$RIN, decreasing = TRUE), ]    
datMeta.HTSC.wholegene.T <- datMeta.HTSC.wholegene.T[!duplicated(datMeta.HTSC.wholegene.T$BrainID), ] # Keep only the first occurrence of each BrainID

# Subset only common genes between datExpr.Cufflinks and datMeta, using snake_case
commonSamples <- intersect(colnames(datExpr.HTSC.wholegene.T), rownames(datMeta.HTSC.wholegene.T))
datExpr.HTSC.wholegene.T <- datExpr.HTSC.wholegene.T[, commonSamples]
datMeta.HTSC.wholegene.T <- datMeta.HTSC.wholegene.T[commonSamples, ]
dim(datMeta.HTSC.wholegene.T)  # Check dimensions to ensure correct subsetting
dim(datExpr.HTSC.wholegene.T)  # Check dimensions to ensure correct subsetting  

pass.whole.genes <- union(rownames(datExpr.HTSC.wholegene.C), union(rownames(datExpr.HTSC.wholegene.F), rownames(datExpr.HTSC.wholegene.T)))

#export_path <- "C:/Users/ariad/OneDrive/Desktop/Proyect WCQN/Raws"

#write.csv(datLB.gene, file = file.path(export_path,"raw_cufflinks.csv"), row.names = TRUE)
#write.csv(datExpr.HTSC.wholegene, file = file.path(export_path,"raw_htseq_wholegene.csv"), row.names = TRUE)
#write.csv(datExpr.HTSC.unionexon, file = file.path(export_path,"raw_htseq_unionexon.csv"), row.names = TRUE)

# 4. Genes that pass at least one of the filters
genes_passed <- Reduce(union, list(pass.cufflinks.genes, pass.unionexon.genes, pass.whole.genes))

length(genes_passed)  # Number of genes before applying gene length filter

# 5. Apply minimum gene length filter (>200bp)
# Load gene length annotation
# load("GC18unionAnno.Rdata") 

file_path<-"C:/Users/ariad/OneDrive/Desktop/Transcriptomics/Genome-wide-changes-in-lncRNA-alternative-splicing-and-cortical-patterning-in-autism-master/data/provided/GC18unionAnno.Rdata"
load(file_path) ## Modified version of Gencode v18 .gtf file to include union exon gene lengths

# Get only genes longer than 20bp
gc18unionAnno_200bp <- gc18unionAnno[gc18unionAnno[,1] > 20, ]
rownames(gc18unionAnno_200bp) <- substr(rownames(gc18unionAnno_200bp), 1, 15)

# Keep only genes present in both gene list and length annotation
genes_final <- intersect(rownames(gc18unionAnno_200bp), genes_passed)

length(genes_final)  # Final number of genes

#6. Let's apply the gene depth normalization to the expression data
# Function to calculate Counts Per Million
calculate_cpm <- function(counts) {
  total_counts <- sum(counts)      # Calculate total reads
  cpm <- (counts / total_counts) #* 1e6  # Normalize and scale by 1 million
  return(cpm)
}

datExpr.HTSC.unionexon.C <- apply(datExpr.HTSC.unionexon.C, 2, calculate_cpm)
datExpr.HTSC.unionexon.F <- apply(datExpr.HTSC.unionexon.F, 2, calculate_cpm)
datExpr.HTSC.unionexon.T <- apply(datExpr.HTSC.unionexon.T, 2, calculate_cpm)


# 7. Let's apply GeTMM (TMM does depth normalization, and lenght normalization) to the expression data

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("edgeR")

#GeTMM is a normalization method to implement using standard R tools, especially:
#edgeR (for TMM normalization)
#Raw counts
#Gene lengths

# Load necessary library
library(edgeR)#TMM normalization

##HERE STARTS THE GeTMM part
GeTMM <- function(counts) {

    gene_lengths <- gc18unionAnno_200bp
  # Ensure gene_lengths is a vector with names matching rownames of counts
  if (!all(rownames(counts) %in% rownames(gene_lengths))) {
    stop("Gene lengths must match the rownames of counts.")
  }
    # Find common genes
    genes_comunes <- intersect(rownames(counts), rownames(gene_lengths))

    # Subsetea ambas matrices para que tengan los mismos genes y en el mismo orden
    lengths_mat_aligned <- gene_lengths[genes_comunes,1]
    expr_mat_aligned    <- counts[genes_comunes, , drop=FALSE]

  # Calculate Reads Per Kilobase (RPK)
  rpk <- expr_mat_aligned / lengths_mat_aligned
  
  # Create DGEList object
  dge <- DGEList(counts = rpk)
  
  # Calculate normalization factors using TMM
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Compute GeTMM-normalized values
  lib_size_scaled <- colSums(rpk) * dge$samples$norm.factors
  geTMM <- t(t(rpk) / lib_size_scaled) * 1e6
  
  return(geTMM)
}


# 9. Log2 transformation (applies to all expression data)
log2Transform <- function(data) {
  return(log2(data + 1))#in case there are zeros in the data but it shouldnt happen because of the previous filtering
}


## For Union Exon
inter <- intersect(rownames(datExpr.HTSC.unionexon.C), genes_final)
datExpr.HTSC.unionexon.C.filtered <- na.omit(datExpr.HTSC.unionexon.C[inter, ])
inter <- intersect(rownames(datExpr.HTSC.unionexon.F), genes_final)
datExpr.HTSC.unionexon.F.filtered <- na.omit(datExpr.HTSC.unionexon.F[inter, ])
inter <- intersect(rownames(datExpr.HTSC.unionexon.T), genes_final)
datExpr.HTSC.unionexon.T.filtered <- na.omit(datExpr.HTSC.unionexon.T[inter, ])


datExpr.HTSC.unionexon.C.filtered <- GeTMM(datExpr.HTSC.unionexon.C.filtered)
datExpr.HTSC.unionexon.F.filtered <- GeTMM(datExpr.HTSC.unionexon.F.filtered)
datExpr.HTSC.unionexon.T.filtered <- GeTMM(datExpr.HTSC.unionexon.T.filtered)

datExpr.HTSC.unionexon.C.filtered <- log2Transform(datExpr.HTSC.unionexon.C.filtered)
datExpr.HTSC.unionexon.F.filtered <- log2Transform(datExpr.HTSC.unionexon.F.filtered)
datExpr.HTSC.unionexon.T.filtered <- log2Transform(datExpr.HTSC.unionexon.T.filtered)

#The metadata (datMeta) is filtered to include only the samples present in the normalized expression matrix.
datMeta.unionexon.C <- datMeta.HTSC.unionexon.C[match(colnames(datExpr.HTSC.unionexon.C.filtered),rownames(datMeta.HTSC.unionexon.C)),] 
datMeta.unionexon.F <- datMeta.HTSC.unionexon.F[match(colnames(datExpr.HTSC.unionexon.F.filtered),rownames(datMeta.HTSC.unionexon.F)),]
datMeta.unionexon.T <- datMeta.HTSC.unionexon.T[match(colnames(datExpr.HTSC.unionexon.T.filtered),rownames(datMeta.HTSC.unionexon.T)),] 

# Save the filtered and normalized data
output_path <- "C:/Users/ariad/OneDrive/Desktop/Proyecto/Exports"
write.csv(datExpr.HTSC.unionexon.C.filtered, file = file.path(output_path, "datExpr.HTSC.unionexon.C.filtered.csv"), row.names = TRUE)
write.csv(datExpr.HTSC.unionexon.F.filtered, file = file.path(output_path, "datExpr.HTSC.unionexon.F.filtered.csv"), row.names = TRUE)
write.csv(datExpr.HTSC.unionexon.T.filtered, file = file.path(output_path, "datExpr.HTSC.unionexon.T.filtered.csv"), row.names = TRUE)
write.csv(datMeta.unionexon.C, file = file.path(output_path, "datMeta.unionexon.C.csv"), row.names = TRUE)
write.csv(datMeta.unionexon.F, file = file.path(output_path, "datMeta.unionexon.F.csv"), row.names = TRUE)
write.csv(datMeta.unionexon.T, file = file.path(output_path, "datMeta.unionexon.T.csv"), row.names = TRUE)  
