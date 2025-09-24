# lmm_pipeline.R
# Mirrors your Python pipeline: fit one LMM per gene using lmer
# Requirements:
# install.packages(c("lme4","lmerTest","tidyverse","qvalue","broom.mixed","ggplot2"))

library(lme4)
library(lmerTest)       # provides p-values for lmer
library(tidyverse)
library(qvalue)
library(broom.mixed)    # tidy summaries from mixed models
library(ggplot2)

# Adjust these paths
expr_path <- "C:/Users/ariad/OneDrive/Desktop/Proyecto/DGEAnalysis/Exports"
output_path <- "C:/Users/ariad/OneDrive/Desktop/Proyecto/RLMMResults"

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# helper to map dataset name (mirrors your Python function)
get_dataset_name <- function(file_path) {
  namee <- basename(file_path) %>%
    str_replace("datExpr.HTSC.unionexon.", "") %>%
    str_replace(".filtered.csv", "")
  if (namee %in% c("C", "CBL")) return("Vermis")
  if (namee %in% c("F", "FRT")) return("Frontal")
  if (namee %in% c("T", "TEM")) return("Temporal")
  return(namee)
}

# find files
expr_files <- list.files(expr_path, pattern = "datExpr.*\\.csv$", full.names = TRUE) %>% sort()
meta_files <- list.files(expr_path, pattern = "datMeta.*\\.imp\\.csv$", full.names = TRUE) %>% sort()

if (length(expr_files) == 0) stop("No expression files found in path.")
if (length(meta_files) == 0) stop("No metadata files found in path.")
if (length(expr_files) != length(meta_files)) {
  message("Warning: unequal number of expression and metadata files. They will be zipped in order.")
}

for (i in seq_along(expr_files)) {
  expr_file <- expr_files[i]
  meta_file <- meta_files[i]
  dataset_name <- get_dataset_name(expr_file)
  cat("\n\n====== Processing Dataset:", dataset_name, "======\n")
  
  # read
  expr_df <- read.csv(expr_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- read.csv(meta_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  # transpose expression so rows = samples, cols = genes
  expr_df <- as.data.frame(t(expr_df))
  
  # align samples (use intersection)
  common_samples <- intersect(rownames(meta), rownames(expr_df))
  expr_df <- expr_df[common_samples, , drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  
  # Prepare covariates (mirror Python scaling: center / (2*sd))
  meta$Age_c <- (meta$Age - mean(meta$Age, na.rm = TRUE)) / (2 * sd(meta$Age, na.rm = TRUE))
  meta$Age2 <- meta$Age_c^2
  meta$Diagnosis <- meta$ASD.CTL
  
  cat("Metadata shape:", dim(meta), "Expression shape:", dim(expr_df), "Samples:", length(common_samples), "\n")
  
  results_list <- vector("list", ncol(expr_df))
  names(results_list) <- colnames(expr_df)
  first_plots <- 0
  
  for (g in colnames(expr_df)) {
    Y <- expr_df[[g]]
    df <- meta %>% select(Diagnosis, Age_c, Age2, BrainBank, Sex, SeqBatch, Seizures) %>% mutate(Y = Y)
    
    # require at least 20 samples like you did
    if (nrow(df) < 20) next
    
    # make sure categorical covariates are factors
    df$Diagnosis <- as.factor(df$Diagnosis)
    df$BrainBank <- as.factor(df$BrainBank)
    df$Sex <- as.factor(df$Sex)
    df$SeqBatch <- as.factor(df$SeqBatch)
    df$Seizures <- as.factor(df$Seizures)
    
    # try to fit
    fit <- tryCatch(
      lmer(Y ~ Diagnosis + Age_c + Age2 + BrainBank + Seizures + Sex + (1 | SeqBatch), data = df, REML = FALSE),
      error = function(e) e,
      warning = function(w) w
    )
    
    if (inherits(fit, "error")) {
      # skip and continue
      message("Skipped ", g, " due to error: ", conditionMessage(fit))
      next
    }
    
    # Extract tidy summary for Diagnosis row
    # broom.mixed::tidy with effects="fixed" returns a row per fixed effect
    tidy_sum <- tryCatch(tidy(fit, effects = "fixed", conf.int = FALSE), error = function(e) NULL)
    if (is.null(tidy_sum)) {
      message("Could not tidy model for ", g)
      next
    }
    
    # find Diagnosis row (exact name depends on factor coding; if Diagnosis has two levels, the row name will be Diagnosis<level>)
    # better to look for rows where term contains "Diagnosis"
    diagnosis_rows <- tidy_sum %>% filter(str_detect(term, "^Diagnosis"))
    if (nrow(diagnosis_rows) == 0) {
      # possibly Diagnosis is reference-coded differently; try to fetch contrast using emmeans or linearHypothesis (skipped here)
      message("Diagnosis effect not found for ", g, "; skipping.")
      next
    }
    
    # If multiple rows (multi-level), pick the first (you may want to adapt based on your encoding)
    row_diag <- diagnosis_rows[1, ]
    
    est <- row_diag$estimate
    # broom.mixed may not include p.value by default; lmerTest adds p-values to summary
    # fallback to summary from lmerTest
    s <- summary(fit)
    # attempt to extract the p-value for the Diagnosis term by searching the coef table
    coef_tab <- as.data.frame(s$coefficients)
    # Ensure row exists
    diag_rowname <- rownames(coef_tab)[str_detect(rownames(coef_tab), "^Diagnosis")]
    if (length(diag_rowname) == 0) {
      # fallback: if diagnosis is named exactly "Diagnosis" (rare), else skip
      message("No Diagnosis row in coef table for ", g)
      next
    }
    # If multiple, take the first
    diag_rowname <- diag_rowname[1]
    pval <- coef_tab[diag_rowname, "Pr(>|t|)"]
    est2 <- coef_tab[diag_rowname, "Estimate"]
    
    # ensure we use the Estimate from coef_tab (consistent) and p-value
    results_list[[g]] <- tibble(Gene = g, EffectSize = est2, PValue = pval)
    
    # plotting residuals for the first 5 genes (like your Python)
    if (first_plots < 5) {
      resid <- resid(fit)
      # KDE (use density) + overlay normal curve
      dens <- density(resid, na.rm = TRUE)
      png(filename = file.path(output_path, paste0(dataset_name, "_", g, "_resid_KDE.png")), width = 1200, height = 800, res = 150)
      plot(dens, main = paste0("KDE Plot of Residuals: ", g), xlab = "Residuals")
      x <- seq(min(resid, na.rm = TRUE), max(resid, na.rm = TRUE), length.out = 200)
      lines(x, dnorm(x, mean = mean(resid, na.rm = TRUE), sd = sd(resid, na.rm = TRUE)), lty = 2)
      dev.off()
      
      # QQ plot
      png(filename = file.path(output_path, paste0(dataset_name, "_", g, "_resid_QQ.png")), width = 800, height = 800, res = 150)
      qqnorm(resid); qqline(resid, col = "red")
      title(main = paste0("Q-Q Plot of Residuals: ", g))
      dev.off()
      first_plots <- first_plots + 1
    }
  } # end gene loop
  
  # combine results
  results_df <- bind_rows(results_list) %>% drop_na()
  if (nrow(results_df) == 0) {
    message("No valid models for dataset: ", dataset_name)
    next
  }
  
  # multiple testing corrections
  results_df <- results_df %>%
    mutate(
      pBH = p.adjust(PValue, method = "BH"),
      pBY = p.adjust(PValue, method = "BY"),
      qvalue = qvalue(PValue)$qvalues
    )
  
  # select significant genes similar to your Python: pBH < 0.05 and abs(EffectSize) > 0.8
  sig_genes <- results_df %>% filter(pBH < 0.05 & abs(EffectSize) > 0.8)
  
  # volcano plot columns
  results_df <- results_df %>% mutate(`-log10(pBH)` = -log10(pBH), Significant = (pBH < 0.05 & abs(EffectSize) > 0.8))
  
  # save
  write.csv(results_df, file = file.path(output_path, paste0(dataset_name, "_LMM.csv")), row.names = FALSE)
  write.csv(sig_genes, file = file.path(output_path, paste0(dataset_name, "_LMM_SG.csv")), row.names = FALSE)
  cat("Results saved for", dataset_name, "\n")
  
  # volcano plot
  p <- ggplot(results_df, aes(x = EffectSize, y = `-log10(pBH)`, color = Significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("grey", "red")) +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "green") +
    labs(title = paste("Volcano Plot:", dataset_name), x = "Effect Size (Diagnosis)", y = "-log10 Adjusted p-value (FDR)") +
    theme_minimal()
  
  ggsave(filename = file.path(output_path, paste0(dataset_name, "_VolcanoPlot.png")), plot = p, width = 10, height = 6, dpi = 300)
  
} # end dataset loop

cat("All done.\n")
