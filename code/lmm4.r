library(lme4)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(multcomp)
library(stats)
library(lmerTest)

# Set your paths
expr_path <- "C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports"
output_path <- "C:/Users/ariad/OneDrive/Desktop/Proyecto/LMMResults"

# Get all relevant files
expr_files <- sort(list.files(expr_path, pattern="datExpr.*\\.csv$", full.names=TRUE))
meta_files <- sort(list.files(expr_path, pattern="datMeta.*\\.csv$", full.names=TRUE))

get_dataset_name <- function(file_path) {
  namee <- basename(file_path)
  namee <- gsub("datExpr.HTSC.unionexon.", "", namee)
  namee <- gsub(".filtered.csv", "", namee)
  if (namee %in% c("C", "CBL")) return("Vermis")
  if (namee %in% c("F", "FRT")) return("Frontal")
  if (namee %in% c("T", "TEM")) return("Temporal")
  return(namee)
}

encode_categorical_metadata <- function(meta) {
  # Binary encoding
  if ("Sex" %in% names(meta)) meta$Sex <- ifelse(meta$Sex == "M", 0, 1)
  if ("BrainBank" %in% names(meta)) meta$BrainBank <- ifelse(meta$BrainBank == "ATP", 0, 1)
  if ("ASD.CTL" %in% names(meta)) meta$Diagnosis <- ifelse(meta$ASD.CTL == "CTL", 0, 1)
  return(meta)
}

for (i in seq_along(expr_files)) {
  dataset_name <- get_dataset_name(expr_files[i])
  cat("\n\n====== Processing Dataset:", dataset_name, "======\n")
  
  expr_df <- fread(expr_files[i], data.table=FALSE)
  rownames(expr_df) <- expr_df[,1]
  expr_df <- expr_df[,-1]
  expr_df <- t(expr_df) # samples x genes
  
  meta <- fread(meta_files[i], data.table=FALSE)
  rownames(meta) <- meta[,1]
  meta <- meta[,-1]
  
  # Ensure alignment
  common_samples <- intersect(rownames(meta), rownames(expr_df))
  expr_df <- expr_df[common_samples, , drop=FALSE]
  meta <- meta[common_samples, , drop=FALSE]
  meta <- encode_categorical_metadata(meta)
  meta$Age <- scale(meta$Age)
  meta$Age2 <- meta$Age^2
  
  cat("Metadata shape:", dim(meta), "\n")
  cat("Expression shape:", dim(expr_df), "\n")
  cat("Number of samples:", length(common_samples), "\n")
  
  results <- list()
  
  for (gene in colnames(expr_df)) {
    df <- meta %>%
      select(Diagnosis, Age, Age2, BrainBank, Sex, SeqBatch) %>%
      mutate(Y = expr_df[,gene])
    
    if (nrow(df) < 20) next
    
    # Fit LMM
    fit <- tryCatch({
      lmer(Y ~ Diagnosis + Age + Age2 + BrainBank + Sex + (1|SeqBatch), data=df, REML=FALSE)
    }, error=function(e) NULL)
    
    if (is.null(fit)) {
      cat("Skipped", gene, "due to model error\n")
      next
    }
    
    coef <- summary(fit)$coefficients
    if (!"Diagnosis" %in% rownames(coef)) {
      cat("Skipped", gene, "due to missing Diagnosis effect\n")
      next
    }
    effect <- coef["Diagnosis", "Estimate"]
    pval <- coef["Diagnosis", "Pr(>|t|)"]
    std_effect <- effect / sd(df$Y, na.rm=TRUE)
    
    results[[gene]] <- data.frame(
      Gene=gene,
      EffectSize=effect,
      StdEffect=std_effect,
      PValue=pval
    )
  }
  
  results_df <- bind_rows(results)
  results_df <- na.omit(results_df)
  
  if (nrow(results_df) == 0) {
    cat("No valid models for dataset:", dataset_name, "\n")
    next
  }
  
  # FDR correction
  results_df$pBH <- p.adjust(results_df$PValue, method="BY")
  results_df$Significant <- (results_df$pBH < 0.05) & (abs(results_df$StdEffect) > 0.8)
  results_df$neglog10pBH <- -log10(results_df$pBH)
  
  # Save results
  write.csv(results_df, file=file.path(output_path, paste0(dataset_name, "_LMM.csv")), row.names=FALSE)
  write.csv(results_df[results_df$Significant,], file=file.path(output_path, paste0(dataset_name, "_LMM_SG.csv")), row.names=FALSE)
  cat("Results saved for", dataset_name, "\n")
  
  # Volcano plot
  p <- ggplot(results_df, aes(x=StdEffect, y=neglog10pBH, color=Significant)) +
    geom_point(alpha=0.7) +
    scale_color_manual(values=c("grey", "red")) +
    geom_vline(xintercept=c(-0.8, 0.8), color="blue", linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), color="green", linetype="dashed") +
    labs(title=paste("Volcano Plot:", dataset_name),
         x="Standardized Effect Size (β for Diagnosis)",
         y="-log10 Adjusted p-value (FDR)") +
    theme_minimal()
  ggsave(filename=file.path(output_path, paste0(dataset_name, "_VolcanoPlot.png")), plot=p, width=10, height=6, dpi=300)
}