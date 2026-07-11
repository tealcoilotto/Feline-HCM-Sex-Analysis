# =============================================================================
# Feline Hypertrophic Cardiomyopathy (HCM) Sex-Specific Transcriptome Analysis
# Data: GSE275971 (Joshua et al. 2025)
# Reanalysis focusing on sex-specific transcriptional
# patterns not examined in the original publication
# =============================================================================

# 1. Setup

library(tidyverse)
library(GEOquery)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(ggrepel)
library(pheatmap)

# Install Bioconductor packages (uncomment to download)
# install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "clusterProfiler", "org.Hs.eg.db"))
# install.packages(c("pheatmap", "ggrepel", "patchwork"))

set.seed(12345)

# 2. Load Data

# Files contain pre-computed differential expression results 
# (logFC, FDR) from edgeR analysis performed by the original authors.
# Values are FPKM-normalized. The FDR is the BH-adjusted p-value. 
# No additional filtering or normalization was performed in this
# reanalysis.

# This analysis uses the differential expression statistics published by 
# Joshua et al. (2025). It is important to note that raw read counts were 
# not available in the supplementary files, so differential expression 
# was not recomputed. Instead, this project performs secondary comparative analyses 
# between male and female differential expression results.

# Analysis Scope:
# This reanalysis focuses exclusively on left ventricular (LV) tissue because
# hypertrophic cardiomyopathy (HCM) is characterized primarily by pathological
# remodeling of the left ventricle. The original study also included left atrial (LA)
# tissue. For this reanalysis, the LV was selected as the primary focus as it is
# most directly affected by HCM.
#
# As a future extension, these analyses could be repeated using the
# available LA datasets to determine whether the observed sex-specific
# transcriptional patterns are consistent across cardiac tissues, or if they
# represent chamber-specific responses.

# Download supplementary files from GEO (uncomment to dowload)
# getGEOSuppFiles("GSE275971")

# Create data path
data_path <- "GSE275971"

# Initial data exploration to confirm column structure prior to main analysis
# Uncomment to run, can be commented out after first run
# df <- read_excel(
#  "GSE275971//GSE275971_mRNA_Female_HCM_LV_vs_female_adult_healthy_LV.xlsx")
# 
# Check dataframe to ensure data is loading correctly
# head(df)
# dim(df)
# colnames(df)

# Load LV data, separating by sex
hcm_female <- read_excel(file.path(
  data_path, "GSE275971_mRNA_Female_HCM_LV_vs_female_adult_healthy_LV.xlsx"))
hcm_male <- read_excel(file.path(
  data_path, "GSE275971_mRNA_male_HCM_LV_vs_male_Adult_Healthy_LV.xlsx"))

# 3. Data Preprocessing

# Fix columns stored as characters instead of numerics
hcm_male <- hcm_male %>%
  mutate(PValue = as.numeric(PValue), FDR = as.numeric(FDR), LR = as.numeric(LR))

hcm_female <- hcm_female %>%
  mutate(PValue = as.numeric(PValue), FDR = as.numeric(FDR), LR = as.numeric(LR))

# Define sample columns for male and female data
male_cols <- c("X0457LV", "X0552LV", "MLVD1", "MLVD17", 
               "MLVD18", "MLVD19", "MLVD2", "LV15", 
               "LV22", "LV23B", "LV6", "LV9")

female_cols <- c("FLVD1", "FLVD3", "FLVD4", "FLVD5", 
                 "FLVD6", "FLVD7", "LV10", "LV21", 
                 "LV7", "LV8")

# Sample size summary
cat("\n=== Sample Sizes ===\n")

# Male samples
cat("Male HCM Samples:", sum(startsWith(colnames(hcm_male), "MLVD")), "\n")
cat("Male Healthy Samples:", sum(startsWith(colnames(hcm_male), "LV")) + 
      sum(startsWith(colnames(hcm_male), "X")), "\n")

# Female samples  
cat("Female HCM Samples:", sum(startsWith(colnames(hcm_female), "FLVD")), "\n")
cat("Female Healthy Samples:", sum(startsWith(colnames(hcm_female), "LV")), "\n")

# Check for duplicated gene names in raw data
# Results: 11 duplicates were found in male file, all corresponding to non-coding RNA 
# annotations (Metazoa_SRP, RNaseP_nuc, etc.) mapping to multiple locations
# in the cat genome. The female file does not have any duplicates. Therefore, we expect
# the shared genes of interest to be unaffected by duplicates, but will check again in 
# section 5 to confirm.
cat("Male Duplicated Gene Names:", sum(duplicated(hcm_male$GeneName)), "\n")
cat("Female Duplicated Gene Names:", sum(duplicated(hcm_female$GeneName)), "\n")

# 4. Differential Expression Summary

# Significant gene counts by sex (using an FDR cutoff of .05)
n_male_sig <- sum(hcm_male$FDR < 0.05)
n_female_sig <- sum(hcm_female$FDR < 0.05)

cat("\n=== Significant Genes (FDR < 0.05) ===\n")
cat("Male Significant Genes:", n_male_sig, "\n")
cat("Female Significant Genes:", n_female_sig, "\n")

# Validate against known HCM causal genes from literature (see README)
known_hcm_genes <- c("MYH7", "MYBPC3", "TNNT2", "TNNI3", 
                     "MYL2", "MYL3", "ACTC1", "TPM1")

cat("\n=== Known HCM Genes: Male ===\n")
hcm_male %>%
  filter(GeneName %in% known_hcm_genes) %>%
  dplyr::select(GeneName, logFC, FDR) %>%
  print()

cat("\n=== Known HCM Genes: Female ===\n")
hcm_female %>%
  filter(GeneName %in% known_hcm_genes) %>%
  dplyr::select(GeneName, logFC, FDR) %>%
  print()

# 5. Sex Comparisons

# Find genes that are significant in both male and female cats
# Note: genes are matched by GeneName (gene symbol). Both files originate 
# from the same GEO dataset and identical edgeR contrasts (HCM vs. healthy LV
# within each sex), ensuring logFC values are comparable between 
# sexes, also minimizing the risk of inconsistencies (e.g., GeneName).
# Genes appearing in only one sex are excluded from the shared analysis.

shared <- inner_join(filter(hcm_female, FDR < 0.05),
                     filter(hcm_male, FDR < 0.05),
                     by = "GeneName") %>% 
  mutate(same_direction = sign(logFC.x) == sign(logFC.y)) %>%
  rename(logFC_female = logFC.x,
         logFC_male = logFC.y,
         FDR_female = FDR.x,
         FDR_male = FDR.y)

# Confirm that duplicates do not impact the shared genes of interest
# Results: All male duplicates had FDR > 0.05 and were excluded by the FDR filter
# before the inner_join, which is confirmed by 0 duplicates in 'shared'.
cat("\nDuplicated Gene Names - Shared:", sum(duplicated(shared$GeneName)), "\n")

cat("\n=== Sex Comparisons ===\n")
cat("Male Significant Genes:", n_male_sig, "\n")
cat("Female Significant Genes:", n_female_sig, "\n")
cat("Shared Significant Genes:", nrow(shared), "\n")
cat("Shared Significant Genes as % of Male Genes:  ", 
    round(nrow(shared)/n_male_sig*100, 1), "%\n")
cat("Shared Significant Genes as % of Female Genes:", 
    round(nrow(shared)/n_female_sig*100, 1), "%\n")
cat("Male Only:", n_male_sig  - nrow(shared), "\n")
cat("Female Only:", n_female_sig - nrow(shared), "\n")

# Evaluate directional concordance within the shared significant genes
# TRUE indicates the same direction, FALSE indicates opposite direction
cat("\n=== Directional Concordance ===\n")
print(table(shared$same_direction))

# Examine sex-discordant genes
cat("\n=== Sex-Discordant Genes ===\n")
shared %>%
  filter(same_direction == FALSE) %>%
  dplyr::select(GeneName, logFC_female, logFC_male, FDR_female, FDR_male) %>%
  print()

# Examine MYOZ1 expression by sex (the one significant sex-discordant gene)
cat("\n=== MYOZ1 Sample Expression ===\n")

myoz1_male <- hcm_male %>%
  filter(GeneName == "MYOZ1") %>%
  dplyr::select(GeneName, all_of(male_cols))

myoz1_female <- hcm_female %>%
  filter(GeneName == "MYOZ1") %>%
  dplyr::select(GeneName, all_of(female_cols))

print(myoz1_male)
print(myoz1_female)

# Use Wilcoxon tests to evaluate MYOZ1 differences among disease-sex groups
male_hcm_myoz1 <- as.numeric(myoz1_male[c("MLVD1","MLVD17","MLVD18","MLVD19","MLVD2")])
female_hcm_myoz1 <- as.numeric(myoz1_female[c("FLVD1","FLVD3","FLVD4","FLVD5","FLVD6","FLVD7")])
male_healthy_myoz1 <- as.numeric(myoz1_male[c("X0457LV","X0552LV","LV15","LV22","LV23B","LV6","LV9")])
female_healthy_myoz1 <- as.numeric(myoz1_female[c("LV10","LV21","LV7","LV8")])

wilcoxon_hcm <- wilcox.test(male_hcm_myoz1, female_hcm_myoz1)
wilcoxon_healthy <- wilcox.test(male_healthy_myoz1, female_healthy_myoz1)

# Results: HCM Male vs. HCM Female MYOZ1 expression differs significantly (p = 0.017) with female HCM samples showing
# consistently lower MYOZ1 than male HCM samples. Healthy male and female cats did not exhibit significant differences
# in MYOZ1 expression (p = 0.927). However, this could be a result of small sample sizes with high variability within
# male healthy, female healthy, and male HCM samples, limiting our confidence in this finding (would need a larger cohort to be sure).
cat("HCM Male vs. Female MYOZ1 Wilcoxon p-value:", round(wilcoxon_hcm$p.value, 4), "\n")
cat("Healthy Male vs. Female MYOZ1 Wilcoxon p-value:", round(wilcoxon_healthy$p.value, 4), "\n")


# Top concordantly upregulated genes
cat("\n=== Top Shared Upregulated Genes ===\n")
shared %>%
  filter(same_direction == TRUE, logFC_female > 0) %>%
  arrange(desc(logFC_female)) %>%
  dplyr::select(GeneName, logFC_female, logFC_male, FDR_female, FDR_male) %>%
  head(20) %>%
  print()

# Fold change correlation between sexes
fc_cor <- cor(shared$logFC_female, shared$logFC_male, method = "pearson")
cat("\n=== Fold Change Correlation ===\n")
cat("Pearson r:", round(fc_cor, 3), "\n")

# Correlation CI and p-value
fc_cor_test <- cor.test(shared$logFC_female, shared$logFC_male, method = "pearson")
cat("95% CI:", round(fc_cor_test$conf.int[1], 3), "-", round(fc_cor_test$conf.int[2], 3), "\n")
cat("P-value:", formatC(fc_cor_test$p.value, format = "e", digits = 4), "\n")

# Spearman correlation (non-parametric, appropriate for fold change non-normality)
sp_cor <- cor.test(shared$logFC_female, shared$logFC_male, method = "spearman")
cat("Spearman rho:", round(sp_cor$estimate, 3), "\n")
cat("Spearman P-value:", formatC(sp_cor$p.value, format = "e", digits = 4), "\n")

# Highly significant genes of interest in both sexes (selected based on 
# statistical significance and known relevance to HCM and/or cardiac disease)
genes_of_interest <- c("SPP1", "MYBPC2", "MCEMP1", "ADAM8", "CXCL14", "NPPB")

cat("\n=== Key Genes: Male ===\n")
hcm_male %>%
  filter(GeneName %in% genes_of_interest) %>%
  dplyr::select(GeneName, logFC, FDR) %>%
  arrange(FDR) %>%
  print()

cat("\n=== Key Genes: Female ===\n")
hcm_female %>%
  filter(GeneName %in% genes_of_interest) %>%
  dplyr::select(GeneName, logFC, FDR) %>%
  arrange(FDR) %>%
  print()

# FLVD6 outlier investigation
# FLVD6 was identified in the direct examination of key HCM marker genes.
# Despite being classified as HCM, FLVD6 shows much lower expression of 
# cardiac stress markers (e.g., NPPB = 4.41 vs. mean 1057 in other female HCM samples),
# potentially suggesting earlier stage HCM. Some individual variability was observed 
# among other HCM samples, consistent with known heterogeneity in HCM progression.
cat("\n=== FLVD6 Outlier Investigation ===\n")
hcm_female %>%
  filter(GeneName %in% genes_of_interest) %>%
  dplyr::select(GeneName, FLVD1, FLVD3, FLVD4, FLVD5, FLVD6, FLVD7) %>%
  print()

# Calculate mean NPPB in other female HCM samples (excluding FLVD6)
nppb_values <- hcm_female %>%
  filter(GeneName == "NPPB") %>%
  dplyr::select(FLVD1, FLVD3, FLVD4, FLVD5, FLVD7) %>%
  unlist()

cat("Mean NPPB in other female HCM samples:", round(mean(nppb_values), 1), "\n")
cat("FLVD6 NPPB:", hcm_female %>% filter(GeneName == "NPPB") %>% pull(FLVD6), "\n")

# X0457LV and X0552LV outlier investigation
# Both samples are classified as male healthy controls but showed very
# elevated expressions of HCM marker genes relative to other healthy controls.
# X0457LV showed particularly high values: NPPB = 2424, SPP1 = 254, and CXCL14 = 148, 
# all substantially higher than other healthy males (NPPB range 10-39, SPP1 range 4-11).
# X0552LV similarly showed elevated values: NPPB = 398, SPP1 = 64.9.
# Both samples also showed elevated MYOZ1 expression (X0457LV = 13.7, 
# X0552LV = 7.80. vs mean 2.1 in other healthy male samples).
# 
# Taken together, these profiles strongly suggest subclinical HCM in both
# cats at the time of sampling, as their transcriptional signatures more closely
# resemble confirmed HCM than healthy tissue, despite being classified as healthy. 
# This likely explains their clustering with HCM samples in the correlation analysis. 
# 
# These samples are retained in the analysis but represent an important caveat,
# that the healthy control group may be impacted by undiagnosed HCM cases, which 
# could potentially reduce the evaluated magnitude of differential expression.
cat("\n=== X0457LV and X0552LV Outlier Investigation ===\n")
hcm_male %>%
  filter(GeneName %in% genes_of_interest) %>%
  dplyr::select(GeneName, X0457LV, X0552LV, LV15, LV22, LV23B, LV6, LV9) %>%
  print()

# Extract healthy male MYOZ1 values (excluding X0457LV and X0552LV)
healthy_male_myoz1_values <- hcm_male %>%
  filter(GeneName == "MYOZ1") %>%
  dplyr::select(LV15, LV22, LV23B, LV6, LV9) %>%
  unlist()

# Extract HCM male MYOZ1 values
hcm_male_myoz1_values <- hcm_male %>%
  filter(GeneName == "MYOZ1") %>%
  dplyr::select(MLVD1, MLVD17, MLVD18, MLVD19, MLVD2) %>%
  unlist()

# Extract healthy male SPP1 values (excluding X0457LV and X0552LV)
healthy_male_spp1_values <- hcm_male %>%
  filter(GeneName == "SPP1") %>%
  dplyr::select(LV15, LV22, LV23B, LV6, LV9) %>%
  unlist()

# Extract HCM male SPP1 values
hcm_male_spp1_values <- hcm_male %>%
  filter(GeneName == "SPP1") %>%
  dplyr::select(MLVD1, MLVD17, MLVD18, MLVD19, MLVD2) %>%
  unlist()

# Extract healthy male NPPB values (excluding X0457LV and X0552LV)
healthy_male_nppb_values <- hcm_male %>%
  filter(GeneName == "NPPB") %>%
  dplyr::select(LV15, LV22, LV23B, LV6, LV9) %>%
  unlist()

# Extract HCM male NPPB values
hcm_male_nppb_values <- hcm_male %>%
  filter(GeneName == "NPPB") %>%
  dplyr::select(MLVD1, MLVD17, MLVD18, MLVD19, MLVD2) %>%
  unlist()

# Calculate mean values and compare to outlier values
cat("Mean MYOZ1 in other healthy male samples:", round(mean(healthy_male_myoz1_values), 1), "\n")
cat("Mean MYOZ1 in hcm male samples:", round(mean(hcm_male_myoz1_values), 1), "\n")
cat("X0457LV MYOZ1:", hcm_male %>% filter(GeneName == "MYOZ1") %>% pull(X0457LV), "\n")
cat("X0552LV MYOZ1:", hcm_male %>% filter(GeneName == "MYOZ1") %>% pull(X0552LV), "\n")

cat("Mean SPP1 in other healthy male samples:", round(mean(healthy_male_spp1_values), 1), "\n")
cat("Mean SPP1 in hcm male samples:", round(mean(hcm_male_spp1_values), 1), "\n")
cat("X0457LV SPP1:", hcm_male %>% filter(GeneName == "SPP1") %>% pull(X0457LV), "\n")
cat("X0552LV SPP1:", hcm_male %>% filter(GeneName == "SPP1") %>% pull(X0552LV), "\n")

cat("Mean NPPB in other healthy male samples:", round(mean(healthy_male_nppb_values), 1), "\n")
cat("Mean NPPB in hcm male samples:", round(mean(hcm_male_nppb_values), 1), "\n")
cat("X0457LV NPPB:", hcm_male %>% filter(GeneName == "NPPB") %>% pull(X0457LV), "\n")
cat("X0552LV NPPB:", hcm_male %>% filter(GeneName == "NPPB") %>% pull(X0552LV), "\n")

# 6. Data Visualizations

# Plot 1: Volcano plots (side by side)

# Top genes to label for each sex
top_labels_male <- hcm_male %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(10)

top_labels_female <- hcm_female %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(10)

# Ensuring even scaling for plots
max_y <- max(-log10(hcm_male$FDR), -log10(hcm_female$FDR), na.rm = TRUE)

# Male volcano plot
p1 <- hcm_male %>%
  mutate(direction = case_when(
    FDR < 0.05 & logFC > 0 ~ "Up",
    FDR < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS")) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), color = direction)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_label_repel(data = top_labels_male,
                   aes(x = logFC, y = -log10(FDR), label = GeneName),
                   size = 3,
                   color = "black",
                   max.overlaps = 10) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  ylim(0, max_y) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = "Male HCM vs Healthy LV",
       x = "log2 Fold Change", y = "-log10 FDR")

# Female volcano plot
p2 <- hcm_female %>%
  mutate(direction = case_when(
    FDR < 0.05 & logFC > 0 ~ "Up",
    FDR < 0.05 & logFC < 0 ~ "Down",
    TRUE ~ "NS")) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), color = direction)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_label_repel(data = top_labels_female,
                   aes(x = logFC, y = -log10(FDR), label = GeneName),
                   size = 3,
                   color = "black",
                   max.overlaps = 10) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  ylim(0, max_y) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(title = "Female HCM vs Healthy LV",
       x = "log2 Fold Change", y = "-log10 FDR")

# Display side by side
p1 + p2

# Plot 2: Visualize directional concordance of shared genes
ggplot(shared, aes(x = logFC_female, y = logFC_male, color = same_direction)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "logFC in Female HCM",
       y = "logFC in Male HCM",
       title = "Shared Differentially Expressed Genes: Male vs Female HCM Response",
       color = "Same Direction") +
  theme_minimal()

# Plot 3: Fold change correlation plot
ggplot(shared, aes(x = logFC_female, y = logFC_male)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", 
           x = min(shared$logFC_female) + 2, 
           y = max(shared$logFC_male) - 2,
           label = paste0("r = ", round(fc_cor, 2)),
           size = 5, color = "red") +
  theme_minimal() +
  labs(x = "Female logFC", y = "Male logFC",
       title = "Fold Change Correlation: Male vs Female HCM")

# Plot 4: Heatmap of top 30 shared genes

# Get the top 30 shared genes by ranking normalized male and female FDR
# (normalized to account for discrepancies in FDR values by sex)
top_genes <- shared %>%
  mutate(rank_female = -log10(FDR_female) / max(-log10(FDR_female)),
         rank_male = -log10(FDR_male) / max(-log10(FDR_male)),
         rank_score = rank_female + rank_male) %>%
  arrange(desc(rank_score)) %>%
  head(30) %>%
  pull(GeneName)

# Create expression matrix
male_expr <- hcm_male %>%
  filter(GeneName %in% top_genes) %>%
  dplyr::select(GeneName, all_of(male_cols)) %>%
  column_to_rownames("GeneName")

female_expr <- hcm_female %>%
  filter(GeneName %in% top_genes) %>%
  dplyr::select(GeneName, all_of(female_cols)) %>%
  column_to_rownames("GeneName")

combined_expr <- cbind(male_expr, female_expr)

# Create annotation dataframe
col_annotation <- data.frame(Sex = c(rep("Male", length(male_cols)),
                                     rep("Female", length(female_cols))),
                             Disease = ifelse(c(male_cols, female_cols) %in% 
                                              c("MLVD1", "MLVD17", "MLVD18", "MLVD19", "MLVD2",
                                                "FLVD1", "FLVD3", "FLVD4", "FLVD5", "FLVD6", "FLVD7"),
                                              "HCM", "Healthy"))

rownames(col_annotation) <- colnames(combined_expr)

# Scale and plot
# Note: expression columns (FLVD1, MLVD1, etc.) contain FPKM-normalized 
# expression values from the original GEO supplementary files.
# These are row-scaled (z-scores) within pheatmap strictly for visualization.
# Default clustering: Euclidean distance, complete linkage
pheatmap(combined_expr,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = col_annotation,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Top 30 Shared Differentially Expressed Genes: Male vs Female HCM")


# Plot 5: Explore MYOZ1 expression by sex and disease status
myoz1_data <- data.frame(Expression = c(male_hcm_myoz1, male_healthy_myoz1, female_hcm_myoz1, female_healthy_myoz1),
                         Group = c(rep("Male HCM", 5), rep("Male Healthy", 7),
                                   rep("Female HCM", 6), rep("Female Healthy", 4)))

ggplot(myoz1_data, aes(x = Group, y = Expression, color = Group)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.5) +
  theme_minimal() +
  labs(title = "MYOZ1 Expression by Sex and Disease Status",
       x = NULL, y = "FPKM") +
  theme(legend.position = "none") +
  scale_x_discrete(limits = c("Male Healthy", "Male HCM",
                              "Female Healthy", "Female HCM"))

# 7. Pathway Enrichment

# Get the gene list for shared significant genes
gene_list <- shared$GeneName

# Convert gene symbols to Entrez IDs (required for clusterProfiler)
# Note: human ortholog annotation (org.Hs.eg.db) used as a proxy due to 
# incomplete feline annotation coverage. This may introduce mapping bias
# and results from clusterProfiler may be incomplete for cat-specific genes.
gene_ids <- bitr(gene_list, 
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# Mapping efficiency
cat("Mapped", nrow(gene_ids), "of", length(gene_list), 
    "shared genes to human Entrez IDs (", 
    round(nrow(gene_ids)/length(gene_list)*100, 1), "%).\n")

# Define background universe IDs as all genes tested in both sexes.
# The background universe helps mitigate the complexities of using human
# Entrez IDs as an approximation.
all_genes_tested <- union(hcm_male$GeneName, hcm_female$GeneName)

universe_ids <- bitr(all_genes_tested,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

# Run GO enrichment
go_results <- enrichGO(gene = gene_ids$ENTREZID,
                       universe = universe_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH", # Benjamini-Hochberg
                       pvalueCutoff = 0.05,
                       readable = TRUE)

# View top results
head(go_results, 20)

# Plot
dotplot(go_results, showCategory = 20) +
  ggtitle("GO Biological Process Enrichment: Shared HCM Differentially Expressed Genes") +
  theme(axis.text.y = element_text(size = 8))

# Upregulated genes only
# Note: using female logFC for direction (99.9% concordance with male logFC confirmed in Section 5)
up_genes <- shared %>% filter(logFC_female > 0) %>% pull(GeneName)

# Downregulated genes only  
down_genes <- shared %>% filter(logFC_female < 0) %>% pull(GeneName)

# Convert IDs and run enrichment separately for upregulated and downregulated genes
up_ids <- bitr(up_genes, fromType="SYMBOL", 
               toType="ENTREZID", OrgDb=org.Hs.eg.db)

down_ids <- bitr(down_genes, fromType="SYMBOL", 
                 toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Upregulated pathway enrichment
go_up <- enrichGO(gene = up_ids$ENTREZID,
                  universe = universe_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP", # Biological Process
                  pAdjustMethod = "BH", # Benjamini-Hochberg
                  pvalueCutoff = 0.05,
                  readable = TRUE)

# Downregulated pathway enrichment: Not plotted in README due to insignificance.
# There was no significant GO biological process enrichment detected after the
# background universe correction at the standard FDR < 0.05 threshold.
# However, with a relaxed threshold of FDR < 0.1, borderline enrichment is
# observed for cardiac structure terms (heart morphogenesis, contractile
# actin filament, etc.) suggesting partial loss of normal cardiac pathways.
# Eye development terms likely reflect GO annotation overlap between cardiac
# and eye development genes rather than true biological signals. Results should be
# interpreted under the context of a relaxed threshold. This contrasts with the 
# observed strongly enriched upregulated pathway analysis.
go_down <- enrichGO(gene = down_ids$ENTREZID,
                    universe = universe_ids$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP", # Biological Process
                    pAdjustMethod = "BH", # Benjamini-Hochberg
                    pvalueCutoff = 0.1, # relaxed since .05 had no significant results
                    readable = TRUE)

cat("\nDownregulated GO Pathways Detected (FDR < 0.1):", nrow(as.data.frame(go_down)), "\n")

# Plot
dotplot(go_up, showCategory=15) + ggtitle("GO Biological Process Enrichment: Upregulated Pathways") + theme(axis.text.y = element_text(size = 6))
dotplot(go_down, showCategory=15) + ggtitle("GO Biological Process Enrichment: Downregulated Pathways (FDR < 0.1)") + theme(axis.text.y = element_text(size = 6)) + labs(caption = "Note: relaxed threshold (FDR < 0.1), standard threshold (FDR < 0.05) yielded no significant terms")


# Sex-Specifc GO enrichment
# Genes significant in only one sex may reveal sex-specific biological responses.

# Identify male-only and female-only significant genes
male_only <- anti_join(
  filter(hcm_male, FDR < 0.05),
  filter(hcm_female, FDR < 0.05),
  by = "GeneName")

female_only <- anti_join(
  filter(hcm_female, FDR < 0.05),
  filter(hcm_male, FDR < 0.05),
  by = "GeneName")

cat("\nMale-only significant genes:", nrow(male_only), "\n")
cat("Female-only significant genes:", nrow(female_only), "\n")

# Split by direction (results merging upreguated and downregulated genes were less informative)
male_only_up <- male_only %>% filter(logFC > 0) %>% pull(GeneName)
male_only_down <- male_only %>% filter(logFC < 0) %>% pull(GeneName)

female_only_up <- female_only %>% filter(logFC > 0) %>% pull(GeneName)
female_only_down <- female_only %>% filter(logFC < 0) %>% pull(GeneName)

cat("Male-only upregulated genes:", length(male_only_up), "\n")
cat("Male-only downregulated genes:", length(male_only_down), "\n")
cat("Female-only upregulated genes:", length(female_only_up), "\n")
cat("Female-only downregulated genes:", length(female_only_down), "\n")

# Convert directional gene sets to Entrez IDs
# Note: human ortholog annotation (org.Hs.eg.db) used as a proxy due to 
# incomplete feline annotation coverage. This may introduce mapping bias
# and results from clusterProfiler may be incomplete for cat-specific genes.
male_only_up_ids <- bitr(male_only_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
male_only_down_ids <- bitr(male_only_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
female_only_up_ids <- bitr(female_only_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
female_only_down_ids <- bitr(female_only_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Sex-specific background universes
male_universe_ids <- bitr(hcm_male$GeneName,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

female_universe_ids <- bitr(hcm_female$GeneName,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)

# Run GO enrichment for each sex and directional set
go_male_only_up <- enrichGO(gene = male_only_up_ids$ENTREZID,
                            universe = male_universe_ids$ENTREZID,
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", # Biological Process
                            pAdjustMethod = "BH", # Benjamini-Hochberg
                            pvalueCutoff = 0.05,
                            readable = TRUE)

go_male_only_down <- enrichGO(gene = male_only_down_ids$ENTREZID,
                              universe = male_universe_ids$ENTREZID,
                              OrgDb = org.Hs.eg.db, 
                              ont = "BP", # Biological Process
                              pAdjustMethod = "BH", # Benjamini-Hochberg
                              pvalueCutoff = 0.05,
                              readable = TRUE)

go_female_only_up <- enrichGO(gene = female_only_up_ids$ENTREZID,
                              universe = female_universe_ids$ENTREZID,
                              OrgDb = org.Hs.eg.db, 
                              ont = "BP", # Biological Process
                              pAdjustMethod = "BH", # Benjamini-Hochberg
                              pvalueCutoff = 0.05,
                              readable = TRUE)

go_female_only_down <- enrichGO(gene = female_only_down_ids$ENTREZID,
                                universe = female_universe_ids$ENTREZID,
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", # Biological Process
                                pAdjustMethod = "BH", # Benjamini-Hochberg
                                pvalueCutoff = 0.05,
                                readable = TRUE)

cat("Male-only upregulated GO terms:", nrow(as.data.frame(go_male_only_up)), "\n")
cat("Male-only downregulated GO terms:", nrow(as.data.frame(go_male_only_down)), "\n")
cat("Female-only upregulated GO terms:", nrow(as.data.frame(go_female_only_up)), "\n")
cat("Female-only downregulated GO terms:", nrow(as.data.frame(go_female_only_down)), "\n")

# Results:
# Male-only upregulated: 1 GO term (actin filament organization, FDR = 0.025) 
# Too few terms to draw informative conclusions.
#
# Male-only downregulated: significant enrichment for energy metabolism
# (generation of precursor metabolites, cellular respiration, aerobic 
# respiration), suggesting male cats experience more pronounced
# metabolic dysfunction with HCM. Neurological terms also present but
# likely reflect GO annotation overlap.
#
# Female-only upregulated: strong immune activation signal
# (immune effector process, leukocyte differentiation, lymphocyte
# differentiation), suggesting female cats mount a stronger
# adaptive immune response in HCM than males.
#
# Female-only downregulated: predominantly neurological terms
# (synapse organization, axon development), which likely reflect
# GO annotation overlap and are not biologically interpretable in the cardiac context.

# Plot meaningful results only
# Male-only upregulated: only 1 term, therefore not plotted

# Male-only downregulated: strong energy metabolism signal
dotplot(go_male_only_down, showCategory = 15) +
  ggtitle("GO Enrichment: Male-Only Downregulated Genes") +
  theme(axis.text.y = element_text(size = 7))

# Female-only upregulated: strong adaptive immune signal
dotplot(go_female_only_up, showCategory = 15) +
  ggtitle("GO Enrichment: Female-Only Upregulated Genes") +
  theme(axis.text.y = element_text(size = 7))

# Female-only downregulated: not plotted due to significant GO annotation overlap (neurological)

# 8. Quality Control

# Use all shared genes for the QC correlation matrix
all_shared_expr_male <- hcm_male %>%
  filter(GeneName %in% shared$GeneName) %>%
  dplyr::select(GeneName, all_of(male_cols)) %>%
  column_to_rownames("GeneName")

all_shared_expr_female <- hcm_female %>%
  filter(GeneName %in% shared$GeneName) %>%
  dplyr::select(GeneName, all_of(female_cols)) %>%
  column_to_rownames("GeneName")

all_shared_expr <- cbind(all_shared_expr_male, all_shared_expr_female)

# Correlation matrix uses all 730 shared genes for a more stable estimate of sample similarity
# compared to just using the top 30 shared genes. The top 30 genes ranked by normalized FDR were 
# evaluated but produced less clear disease-status clustering.
cor_matrix <- cor(all_shared_expr, method = "pearson")

# Annotation for healthy vs. HCM
disease_annotation <- data.frame(
  Disease = ifelse(colnames(cor_matrix) %in% c("MLVD1", "MLVD17", "MLVD18",
                   "MLVD19", "MLVD2", "FLVD1", "FLVD3", "FLVD4",
                   "FLVD5", "FLVD6", "FLVD7"), "HCM", "Healthy"),
  Sex = ifelse(colnames(cor_matrix) %in% male_cols, "Male", "Female"))

rownames(disease_annotation) <- colnames(cor_matrix)

# Plot
# Default clustering: Euclidean distance, complete linkage
pheatmap(cor_matrix,
         main = "Sample Correlation Matrix",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_col = disease_annotation,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")

# 9. Session Info (for reproducibility)
sink("session_info.txt")
sessionInfo()
sink()
