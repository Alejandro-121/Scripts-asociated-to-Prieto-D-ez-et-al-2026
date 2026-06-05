#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Genotype-based filtering for VCF variant analysis
# Author: alejandro
# Date: 2024

# Load required libraries
library(readr)
library(org.Sc.sgd.db)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tidyr)

# Import data
cat("Loading mutation table...\n")
table <- read_csv("mutation_wide_table_with_alleles.csv")

############################
# Rename problematic cols  #
############################
cat("Renaming columns...\n")

# Rename to change 2-1 and 2-3 to more apt names for R
colnames(table) <- c("CHROM", "POS", "REF", "ALT", "Effect", "Gene", "Transcript",
                     "p2_1_mutation", "p2_3_mutation", "sup.11_mutation", "sup.15_mutation", 
                     "sup.1_mutation", "sup.22_mutation", "sup.23_mutation", "sup.25_mutation", 
                     "sup.27_mutation", "sup.2_mutation", "WT_mutation",
                     "p2_1_allele", "p2_3_allele", "sup.11_allele", "sup.15_allele", 
                     "sup.1_allele", "sup.22_allele", "sup.23_allele", "sup.25_allele", 
                     "sup.27_allele", "sup.2_allele", "WT_allele")

# Reorder table to simplify reading the final output
table_r <- table[, c("CHROM", "POS", "REF", "ALT", "Effect", "Gene", "Transcript",
                     "WT_mutation",
                     "p2_1_mutation",
                     "sup.1_mutation", "sup.2_mutation", "sup.22_mutation", "sup.23_mutation",
                     "p2_3_mutation",
                     "sup.11_mutation", "sup.15_mutation", "sup.25_mutation", "sup.27_mutation",
                     "WT_allele",
                     "p2_1_allele",
                     "sup.1_allele", "sup.2_allele", "sup.22_allele", "sup.23_allele",
                     "p2_3_allele",
                     "sup.11_allele", "sup.15_allele", "sup.25_allele", "sup.27_allele")]

cat("Adding gene annotations from SGD...\n")
annotated <- AnnotationDbi::select(
  org.Sc.sgd.db,
  keys = table_r$Gene,
  columns = c("GENENAME", "DESCRIPTION", "ORF"),
  keytype = "ORF"
)

table_r <- cbind(
  annotated$GENENAME,
  annotated$DESCRIPTION,
  table_r
)

cat("Separating descendants from 2.1 and 2.3...\n")
# Separate descendants from 2.1 and 2.3
# 2-1 -> 1, 2, 22, 23
# 2-3 -> 11, 15, 25, 27

# Remove cols not from 2.1
from_2.1 <- table_r %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.25_mutation,
                -sup.27_mutation,
                -sup.11_allele, -sup.15_allele, -sup.25_allele,
                -sup.27_allele,
                -p2_3_allele,
                -p2_3_mutation)

# Remove cols not from 2.3
from_2.3 <- table_r %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.2_allele, -sup.22_allele, -sup.23_allele,
                -p2_1_allele,
                -p2_1_mutation)

#######
# 001 #
#######
cat("Processing pattern 001 (mutations absent in WT and parent, present in suppressor)...\n")

### 2-1 ###
# sup.1
change_from_2.1_dif_wt_dif_ref_sup1 <- from_2.1 %>%
  filter(WT_mutation == 0, p2_1_mutation == 0, sup.1_mutation == 1) %>%
  dplyr::select(-sup.2_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.2_allele, -sup.22_allele, -sup.23_allele)

# sup.2
change_from_2.1_dif_wt_dif_ref_sup2 <- from_2.1 %>%
  filter(WT_mutation == 0, p2_1_mutation == 0, sup.2_mutation == 1) %>%
  dplyr::select(-sup.1_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.22_allele, -sup.23_allele)

# sup.22
change_from_2.1_dif_wt_dif_ref_sup22 <- from_2.1 %>%
  filter(WT_mutation == 0, p2_1_mutation == 0, sup.22_mutation == 1) %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.2_allele, -sup.23_allele)

# sup.23
change_from_2.1_dif_wt_dif_ref_sup23 <- from_2.1 %>%
  filter(WT_mutation == 0, p2_1_mutation == 0, sup.23_mutation == 1) %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.22_mutation,
                -sup.1_allele, -sup.2_allele, -sup.22_allele)

### 2-3 ###
# sup.11
change_from_2.3_dif_wt_dif_ref_sup11 <- from_2.3 %>%
  filter(WT_mutation == 0, p2_3_mutation == 0, sup.11_mutation == 1) %>%
  dplyr::select(-sup.15_mutation, -sup.25_mutation, -sup.27_mutation,
                -sup.15_allele, -sup.25_allele, -sup.27_allele)

# sup.15
change_from_2.3_dif_wt_dif_ref_sup15 <- from_2.3 %>%
  filter(WT_mutation == 0, p2_3_mutation == 0, sup.15_mutation == 1) %>%
  dplyr::select(-sup.11_mutation, -sup.25_mutation, -sup.27_mutation,
                -sup.11_allele, -sup.25_allele, -sup.27_allele)

# sup.25
change_from_2.3_dif_wt_dif_ref_sup25 <- from_2.3 %>%
  filter(WT_mutation == 0, p2_3_mutation == 0, sup.25_mutation == 1) %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.27_mutation,
                -sup.11_allele, -sup.15_allele, -sup.27_allele)

# sup.27
change_from_2.3_dif_wt_dif_ref_sup27 <- from_2.3 %>%
  filter(WT_mutation == 0, p2_3_mutation == 0, sup.27_mutation == 1) %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.25_mutation,
                -sup.11_allele, -sup.15_allele, -sup.25_allele)

#######
# 101 #
#######
cat("Processing pattern 101 (mutations present in WT and suppressor, absent in parent)...\n")

### 2-1 ###
# sup.1
change_from_2.1_sam_wt_dif_ref_sup1 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 0, sup.1_mutation == 1) %>%
  dplyr::select(-sup.2_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.2_allele, -sup.22_allele, -sup.23_allele)

# sup.2
change_from_2.1_sam_wt_dif_ref_sup2 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 0, sup.2_mutation == 1) %>%
  dplyr::select(-sup.1_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.22_allele, -sup.23_allele)

# sup.22
change_from_2.1_sam_wt_dif_ref_sup22 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 0, sup.22_mutation == 1) %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.2_allele, -sup.23_allele)

# sup.23
change_from_2.1_sam_wt_dif_ref_sup23 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 0, sup.23_mutation == 1) %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.22_mutation,
                -sup.1_allele, -sup.2_allele, -sup.22_allele)

### 2-3 ###
# sup.11
change_from_2.3_sam_wt_dif_ref_sup11 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 0, sup.11_mutation == 1) %>%
  dplyr::select(-sup.15_mutation, -sup.25_mutation, -sup.27_mutation,
                -sup.15_allele, -sup.25_allele, -sup.27_allele)

# sup.15
change_from_2.3_sam_wt_dif_ref_sup15 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 0, sup.15_mutation == 1) %>%
  dplyr::select(-sup.11_mutation, -sup.25_mutation, -sup.27_mutation,
                -sup.11_allele, -sup.25_allele, -sup.27_allele)

# sup.25
change_from_2.3_sam_wt_dif_ref_sup25 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 0, sup.25_mutation == 1) %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.27_mutation,
                -sup.11_allele, -sup.15_allele, -sup.27_allele)

# sup.27
change_from_2.3_sam_wt_dif_ref_sup27 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 0, sup.27_mutation == 1) %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.25_mutation,
                -sup.11_allele, -sup.15_allele, -sup.25_allele)

#######
# 110 #
#######
cat("Processing pattern 110 (mutations present in WT and parent, absent in suppressor - reversions)...\n")

### 2-1 ###
# sup.1
change_from_2.1_dif_wt_sam_ref_sup1 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 1, sup.1_mutation == 0) %>%
  dplyr::select(-sup.2_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.2_allele, -sup.22_allele, -sup.23_allele)

# sup.2
change_from_2.1_dif_wt_sam_ref_sup2 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 1, sup.2_mutation == 0) %>%
  dplyr::select(-sup.1_mutation, -sup.22_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.22_allele, -sup.23_allele)

# sup.22
change_from_2.1_dif_wt_sam_ref_sup22 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 1, sup.22_mutation == 0) %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.23_mutation,
                -sup.1_allele, -sup.2_allele, -sup.23_allele)

# sup.23
change_from_2.1_dif_wt_sam_ref_sup23 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 1, sup.23_mutation == 0) %>%
  dplyr::select(-sup.1_mutation, -sup.2_mutation, -sup.22_mutation,
                -sup.1_allele, -sup.2_allele, -sup.22_allele)

### 2-3 ###
# sup.11
change_from_2.3_dif_wt_sam_ref_sup11 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 1, sup.11_mutation == 0) %>%
  dplyr::select(-sup.15_mutation, -sup.25_mutation, -sup.27_mutation,
                -sup.15_allele, -sup.25_allele, -sup.27_allele)

# sup.15
change_from_2.3_dif_wt_sam_ref_sup15 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 1, sup.15_mutation == 0) %>%
  dplyr::select(-sup.11_mutation, -sup.25_mutation, -sup.27_mutation,
                -sup.11_allele, -sup.25_allele, -sup.27_allele)

# sup.25
change_from_2.3_dif_wt_sam_ref_sup25 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 1, sup.25_mutation == 0) %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.27_mutation,
                -sup.11_allele, -sup.15_allele, -sup.27_allele)

# sup.27
change_from_2.3_dif_wt_sam_ref_sup27 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 1, sup.27_mutation == 0) %>%
  dplyr::select(-sup.11_mutation, -sup.15_mutation, -sup.25_mutation,
                -sup.11_allele, -sup.15_allele, -sup.25_allele)

##########################
# Calculate common mutations
##########################
cat("Calculating common mutations across all suppressors...\n")

# X01 - Common mutations in all suppressors, absent in parent
common_mut_2.1 <- from_2.1 %>%
  filter(p2_1_mutation == 0, sup.1_mutation == 1, sup.2_mutation == 1, 
         sup.22_mutation == 1, sup.23_mutation == 1)

common_mut_2.3 <- from_2.3 %>%
  filter(p2_3_mutation == 0, sup.11_mutation == 1, sup.15_mutation == 1, 
         sup.25_mutation == 1, sup.27_mutation == 1)

# X10 - Common reversions (present in parent, absent in all suppressors)
common_rev_2.1 <- from_2.1 %>%
  filter(p2_1_mutation == 1, sup.1_mutation == 0, sup.2_mutation == 0, 
         sup.22_mutation == 0, sup.23_mutation == 0)

common_rev_2.3 <- from_2.3 %>%
  filter(p2_3_mutation == 1, sup.11_mutation == 0, sup.15_mutation == 0, 
         sup.25_mutation == 0, sup.27_mutation == 0)

cat("\nCommon mutations in 2.1 suppressors:\n")
print(common_mut_2.1)

cat("\nCommon mutations in 2.3 suppressors:\n")
print(common_mut_2.3)

cat("\nCommon reversions in 2.1 suppressors:\n")
print(common_rev_2.1)

cat("\nCommon reversions in 2.3 suppressors:\n")
print(common_rev_2.3)

##########################
# Select common for most interesting patterns
##########################
cat("\nCalculating common patterns 101 and 010...\n")

common_mut_2.1_101 <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 0, sup.1_mutation == 1, 
         sup.2_mutation == 1, sup.22_mutation == 1, sup.23_mutation == 1)

common_mut_2.1_010 <- from_2.1 %>%
  filter(WT_mutation == 0, p2_1_mutation == 1, sup.1_mutation == 0, 
         sup.2_mutation == 0, sup.22_mutation == 0, sup.23_mutation == 0)

common_mut_2.3_101 <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 0, sup.11_mutation == 1, 
         sup.15_mutation == 1, sup.25_mutation == 1, sup.27_mutation == 1)

common_mut_2.3_010 <- from_2.3 %>%
  filter(WT_mutation == 0, p2_3_mutation == 1, sup.11_mutation == 0, 
         sup.15_mutation == 0, sup.25_mutation == 0, sup.27_mutation == 0)

##########################
# Calculate frequency of mutations across suppressors
##########################
cat("\nCalculating mutation frequencies...\n")

common_mut_2.1_101_freq <- from_2.1 %>%
  filter(WT_mutation == 1, p2_1_mutation == 0, 
         sup.1_mutation == 1 | sup.2_mutation == 1 | sup.22_mutation == 1 | sup.23_mutation == 1)

common_mut_2.1_010_freq <- from_2.1 %>%
  filter(WT_mutation == 0, p2_1_mutation == 1, 
         sup.1_mutation == 0 | sup.2_mutation == 0 | sup.22_mutation == 0 | sup.23_mutation == 0)

common_mut_2.3_101_freq <- from_2.3 %>%
  filter(WT_mutation == 1, p2_3_mutation == 0, 
         sup.11_mutation == 1 | sup.15_mutation == 1 | sup.25_mutation == 1 | sup.27_mutation == 1)

common_mut_2.3_010_freq <- from_2.3 %>%
  filter(WT_mutation == 0, p2_3_mutation == 1, 
         sup.11_mutation == 0 | sup.15_mutation == 0 | sup.25_mutation == 0 | sup.27_mutation == 0)

# Calculate how many times each mutation appears
times_010_2.1 <- 4 - (common_mut_2.1_010_freq$sup.1_mutation + 
                      common_mut_2.1_010_freq$sup.2_mutation + 
                      common_mut_2.1_010_freq$sup.22_mutation + 
                      common_mut_2.1_010_freq$sup.23_mutation)

times_101_2.1 <- common_mut_2.1_101_freq$sup.1_mutation + 
                 common_mut_2.1_101_freq$sup.2_mutation + 
                 common_mut_2.1_101_freq$sup.22_mutation + 
                 common_mut_2.1_101_freq$sup.23_mutation

##########################
# Generate plots
##########################
cat("\nGenerating frequency plot for pattern 010 (2.1)...\n")

# Prepare data for plotting
plot_data_010 <- data.frame(
  Gene = common_mut_2.1_010_freq$Gene,
  Freq = times_010_2.1
)

# Save frequency plot
pdf("mutation_frequency_010_2.1.pdf", width = 12, height = 6)
ggplot(plot_data_010, aes(x = Gene, y = Freq)) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Mutation Frequency - Pattern 010 (lineage 2.1)",
       x = "Gene", y = "Number of suppressors")
dev.off()

cat("Generating heatmap for pattern 101 (2.1)...\n")

# Prepare data for heatmap
to_plot_101 <- data.frame(
  Gene_Name = paste(common_mut_2.1_101_freq$Gene, 
                    common_mut_2.1_101_freq$`annotated$GENENAME`),
  sup1 = as.numeric(common_mut_2.1_101_freq$sup.1_mutation),
  sup2 = as.numeric(common_mut_2.1_101_freq$sup.2_mutation),
  sup22 = as.numeric(common_mut_2.1_101_freq$sup.22_mutation),
  sup23 = as.numeric(common_mut_2.1_101_freq$sup.23_mutation)
)

# Reshape to long format
df_long <- to_plot_101 %>%
  pivot_longer(cols = starts_with("sup"), 
               names_to = "SAMPLE", 
               values_to = "VALUE")

# Save heatmap
pdf("mutation_heatmap_101_2.1.pdf", width = 10, height = 8)
ggplot(df_long, aes(x = SAMPLE, y = Gene_Name, fill = VALUE)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(title = "Mutation Heatmap - Pattern 101 (lineage 2.1)", 
       x = "Samples", y = "Genes",
       fill = "Mutation\nPresent")
dev.off()

##########################
# Save filtered results
##########################
cat("\nSaving filtered results to CSV files...\n")

write_csv(common_mut_2.1, "common_mutations_2.1.csv")
write_csv(common_mut_2.3, "common_mutations_2.3.csv")
write_csv(common_rev_2.1, "common_reversions_2.1.csv")
write_csv(common_rev_2.3, "common_reversions_2.3.csv")
write_csv(common_mut_2.1_101, "common_mutations_101_2.1.csv")
write_csv(common_mut_2.1_010, "common_mutations_010_2.1.csv")
write_csv(common_mut_2.3_101, "common_mutations_101_2.3.csv")
write_csv(common_mut_2.3_010, "common_mutations_010_2.3.csv")

cat("\n=== Analysis completed successfully! ===\n")
cat("\nGenerated files:\n")
cat("  - common_mutations_2.1.csv\n")
cat("  - common_mutations_2.3.csv\n")
cat("  - common_reversions_2.1.csv\n")
cat("  - common_reversions_2.3.csv\n")
cat("  - common_mutations_101_2.1.csv\n")
cat("  - common_mutations_010_2.1.csv\n")
cat("  - common_mutations_101_2.3.csv\n")
cat("  - common_mutations_010_2.3.csv\n")
cat("  - mutation_frequency_010_2.1.pdf\n")
cat("  - mutation_heatmap_101_2.1.pdf\n")


