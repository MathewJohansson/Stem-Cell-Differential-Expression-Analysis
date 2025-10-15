

# 1. Analysing Transcriptomic Data. --------------------------------------------


# This is a differential expression analysis comparing 
#   HSPC = Haematopoietic Stem and Progenitor Cells (intermediate stage)
#     with
#   Prog = Progenitors (more differentiated, committed cells)
#     to track blood stem cell differentiation. 

# Research question: Which genes drive differentiation in haematopoietic stem cells? 




# Load packages. 
library(tidyverse)
library(readxl)
library(dplyr)
library(ggrepel)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(GenomicRanges::setdiff)

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("scran")
library(scran)

# BiocManager::install("biomaRt")
library(biomaRt)



# Load the datasets and store them as variables to be used throughout the analysis.
hspc  <- read_csv("/home/johansson/Documents/Bioinformatics/Portfolio Projects/Project 1. Differential Expression Analysis/Stem Cell Analysis/data-raw/secretome_hspc.csv")
prog  <- read_csv("/home/johansson/Documents/Bioinformatics/Portfolio Projects/Project 1. Differential Expression Analysis/Stem Cell Analysis/data-raw/secretome_prog.csv")


# Join the two datasets together by the shared variable "ensemble_gene_id" to make analysis easier.
hspc_prog <- hspc %>%
  left_join(prog,
            by = "ensembl_gene_id")




# 1.1. DISTRIBUTIONS. ----------------------------------------------------------

# Distribution of values across all the data in the file.

# Basic histogram plot to display distribution. 
hspc_prog %>%
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") %>%
  ggplot(aes(x = expr)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Expression", y = "Count") +
  theme_classic()

# 


# Distribution of values across the samples, as a summary table. 
hspc_prog_summary_cell <- hspc_prog %>%
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") %>%
  group_by(cell) %>%
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            sd = sd(expr),
            mean = mean(expr),
            median = median(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            total = sum(expr),
            n_above_zero = sum(expr > 0))
hspc_prog_summary_cell


# Show distribution of expressions in cells, as a basic bar plot. 
hspc_prog_summary_cell %>%
  ggplot(aes(x = cell, y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd),
                  size = 0.1)

# Arrange 'cell' in increasing size of 'mean', as a basic bar plot. 
hspc_prog_summary_cell %>%
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  size = 0.1)


# Distribution of values across the genes, as a summary table. 
hspc_prog_summary_gene <- hspc_prog %>%
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") %>%
  group_by(ensembl_gene_id) %>%
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            sd = sd(expr),
            mean = mean(expr),
            median = median(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            total = sum(expr),
            n_above_zero = sum(expr > 0))
hspc_prog_summary_gene 


# Plotting the logged mean counts for each gene in order of size; genes, unlike cells, are expected 
#   to have large differences in mean expression levels. 
hspc_prog_summary_gene %>%
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  size = 0.1)




# 1.2. QUALITY CONTROL. --------------------------------------------------------

# Find the genes that are 0 in every column of the hspc_prog dataframe (no expression).
hspc_prog %>%
  rowwise() %>%
  filter(sum(c_across(HSPC_001:Prog_852)) == 0)


# Write filtered data to file. 
write_csv(hspc_prog,
          file = "/home/johansson/Documents/Bioinformatics/Portfolio Projects/Project 1. Differential Expression Analysis/Stem Cell Analysis/data-processed/hspc_prog.csv")




# 2. STATISTICAL/DIFFERENTIAL EXPRESSION ANALYSIS. -----------------------------

# 2.1. GENES EXPRESSED IN ONE CELL TYPE. ---------------------------------------


# Genes expressed only in the progenitor cells (those that are 0 in every HSPC cell). 
hspc_prog_only <- hspc_prog %>%
  rowwise() %>%
  filter(sum(c_across(HSPC_001:HSPC_852)) == 0)
hspc_prog_only   # = 0 Prog-only genes.

# Save results.
write_csv(hspc_prog_only, "results/hspc_prog_only.csv")


# Genes expressed only in the HSPC cells. 
hspc_prog_hspc_only <- hspc_prog %>% 
  rowwise() %>% 
  filter(sum(c_across(Prog_001:Prog_852)) == 0)
hspc_prog_hspc_only   # = 0 HSPC-only genes. 

# Save results.
write_csv(hspc_prog_hspc_only, "results/hspc_prog_hspc_only.csv")




# 2.2. PREPARE THE DATA FOR ANALYSIS WITH SCRAN. -------------------------------


# Add the gene ids as the row names, to include more genetic information to the data. 
hspc_prog <- hspc_prog |>
  column_to_rownames("ensembl_gene_id")
hspc_prog


# Create a vector that indicates which column belongs to which cell type. 
n_hspc <- 701
n_prog <- 798

cell_type <- rep(c("hspc","prog"), 
                 times = c(n_hspc, n_prog))
cell_type




# 2.3. DIFFERENTIAL EXPRESSION ANALYSIS. ---------------------------------------


# The scran package is required here. 


# Run the differential expression analysis. 
results_hspc_prog <- findMarkers(hspc_prog, 
                                 cell_type)
results_hspc_prog

# Outputs two dataframes.


results_hspc_prog$prog
# This dataframe is log prog - log hspc (i.e., Prog/HSPC). 
# This means:
#   Pos fold change = prog > hspc. 
#   Neg fold change = hspc > prog.


results_hspc_prog$hspc
# This dataframe is log hspc - log prog.
# This means:
#   Pos fold change = hspc > prog. 
#   Neg fold change = prog > hspc. 

# Both have 423 genes and 5 columns - 1 for gene number, 4 for statistics. 


# Choose one of the dataframes to analyse/display. 

# Extract the results dataframe from the list object and adds gene ids as col. 
hspc_prog_results <- data.frame(results_hspc_prog$prog,
                                ensembl_gene_id = 
                                  row.names(results_hspc_prog$prog))
hspc_prog_results


# Return the ensembl gene ids as a column to the normalised counts. 
hspc_prog <- hspc_prog %>%
  rownames_to_column(var = "ensembl_gene_id")
hspc_prog

hspc_prog_results <- hspc_prog_results %>%
  left_join(hspc_prog, by = "ensembl_gene_id")
hspc_prog_results


# Gene information can then be added. 

# Ensure that the biomaRt package is installed at the top and load the package.

# Connect to the mouse database to access the relevant data on mice genes.
ensembl <- useMart(biomart = "ensembl",
                   dataset = "mmusculus_gene_ensembl")
ensembl


# Display the information on genes that we can retrieve. 
listAttributes(mart = ensembl) %>%
  View()
# = 3,001 possible bits of information showing name and description. 


# Get the gene name, description, and the gene id to enable joining the information with the results. 
gene_info <- getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id",
                                  "external_gene_name",
                                  "description"),
                   values = hspc_prog_results$ensembl_gene_id,
                   mart = ensembl)
gene_info


# Merge gene information with the results. 
hspc_prog_results <- hspc_prog_results %>%
  left_join(gene_info, by = "ensembl_gene_id")


# Save the results to a .csv file in the results folder. 
write_csv(hspc_prog_results, file = "results/hspc_prog_results.csv")




# 3. VISUALISING AND INTERPRETING. ---------------------------------------------


# 3.1. VIEW DATA. --------------------------------------------------------------

glimpse(hspc_prog_results)
#    423 rows - genes.
#  1,507 columns - statistics, gene info, and all 1,499 individual cell expression values.




# 3.2. WRITE THE SIGNIFICANT GENES TO FILE. ------------------------------------

# Filter for significantly differentially expressed genes (FDR ≤ 0.05).
# This is done because only a subset will be truly differentially expressed. 
# Random variation isn't biologically meaningful, shed the noise.
# 
hspc_prog_results_sig0.05 <- hspc_prog_results %>%
  filter(FDR <= 0.05)
hspc_prog_results_sig0.05

write_csv(hspc_prog_results_sig0.05,
          file = "results/hspc_prog_results_sig0.05.csv")




# 3.3. PRINCIPAL COMPONENT ANALYSIS (PCA). -------------------------------------


# Select for results starting with either "HSPC_" or "Prog_".
# t() transforms (required for PCA). 
# Data is then set into a dataframe. 
hspc_prog_trans <- hspc_prog_results %>%
  dplyr::select(starts_with(c("HSPC_", "Prog_"))) %>%
  t() %>%
  data.frame()
hspc_prog_trans


# Set column names as gene ids. 
colnames(hspc_prog_trans) <- hspc_prog_results$ensembl_gene_id
colnames(hspc_prog_trans)


# Perform PCA on the log2 transformed normalised counts. 
pca <- hspc_prog_trans %>%
  prcomp(rank. = 8)
pca

# Summary data. 
summary(pca)


# Create a dataframe of the PCA scores from PC1 and PC2, and add cell_id column of the row names. 
pca_labelled <- data.frame(pca$x,
                           cell_id = row.names(hspc_prog_trans))
pca_labelled


# Extract cell type and number from the cell_id column; useful as we're colouring by cell type. 
pca_labelled <- pca_labelled %>%
  extract(cell_id, 
          remove = FALSE,
          c("cell_type", "cell_number"),
          "([a-zA-Z]{4})_([0-9]{3})")
pca_labelled

# The regular expression "([a-zA-Z..." etc. is explained at the end of this pipeline. 


# Plot the PCA for PC1 and PC2. 
pca_labelled %>%
  ggplot(aes(x= PC1, y = PC2,
             colour = cell_type)) +
  geom_point(alpha = 0.4) +
  theme_classic()
# General cluster visible, but plenty of overlap. 


# Plot the PCA for PC3 and PC4 to display another comparison.
pca_labelled %>%
  ggplot(aes(x = PC3, y = PC4,
             colour = cell_type)) +
  geom_point(alpha = 0.4) +
  theme_classic()
# Greater overlap of PCs in this comparison, less distinction between clusters. 


# Plot PC1 against PC2 again, coloured by cell type. 
# Colours here (viridis) are to aid viewers with colour blindness. 
pca_labelled %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_colour_viridis_d(end = 0.95, begin = 0.15,
                         name = "Cell Type",
                         labels = c("HSPC", "Progenitor")) +
  labs(x = "PC1 (10.7% variance)",
       y = "PC2 (5.5% variance)",
       title = "PCA of HSPC and Progenitor Cells") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11))

# Save the PCA plot in the figures folder. 
ggsave("figures/pca_hspc_prog.png",
       height = 5,
       width = 6,
       units = "in",
       dpi = 300)




# 3.4. VOLCANO PLOTS. -------------------------------------


# Take the hspc_prog_results dataframe and add in a column for log10_FDR, required for volcano plots.
hspc_prog_results <- hspc_prog_results %>%
  mutate(log10_FDR = -log10(FDR))
hspc_prog_results


# Plot a volcano plot for the results. These plots show: 
#   1. Statistical significance   - how confident the change is real. 
#   2. Biological significance    - how large the change actually is. 
# This allows for identification of genes that are significantly different: 
#               High on the y-axis, on the far left or far right.
#               Genes there have both large fold changes and are statistically robust. 
# 
hspc_prog_results %>%
  ggplot(aes(x = summary.logFC,
             y = log10_FDR)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")


# Take the hspc_prog_results and add columns for sig and bigfc.
#     sig = TRUE if a gene is statistically significant, FALSE if not. 
#   bigfc = TRUE if a gene has a large fold change, FALSE if not. 
# 
hspc_prog_results <- hspc_prog_results %>%
  mutate(sig = FDR <= 0.05,
         bigfc = abs(summary.logFC) >= 2)
hspc_prog_results 


# Plot a volcano plot with three colours (only 3 needed, so grey30 is commented).
hspc_prog_results %>%
  ggplot(aes(x = summary.logFC,
             y = log10_FDR,
             colour = interaction(sig, bigfc))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("grey",
                                 "pink",
                                 # "grey30",
                                 "deeppink")) +
  theme_classic() +
  theme(legend.position = "none")


# Plot again, adding labels to the datapoints/values that are: 
#   >0.05 for log10_FDR
#   >2 for summary.logFC
#   <-2 for summary.logFC
#     This highlights results we want to focus on. 
# 
hspc_prog_results %>%
  ggplot(aes(x = summary.logFC,
             y = log10_FDR,
             colour = interaction(sig, bigfc))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("grey",
                                 "pink",
                                 # "grey30",
                                 "deeppink")) +
  geom_text_repel(data = hspc_prog_results %>%
                    filter(bigfc == TRUE, sig == TRUE), 
                  aes(label = external_gene_name),
                  size = 3,
                  max.overlaps = 50) +
  theme_classic() +
  theme(legend.position = "none")
# On this plot, only the significant data with large fold changes have been labelled. 

# Remember: 
#   Pos fold change = up-regulation of Prog. 
#   Neg fold change = down-regulation of Prog (higher in HSPC). 

# If you have forgotten which way around the comparison is: 
#   Examine the gene summary dataframe to see which of the treatments seems to be higher
#     for the positive fold changes. 


# Label one specific gene of interest on the volcano plot. 
hspc_prog_results %>%
  ggplot(aes(x = summary.logFC,
             y = log10_FDR,
             colour = interaction(sig, bigfc))) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey40") +
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed",
             color = "grey40") +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  scale_colour_manual(values = c("grey60", "coral", "firebrick"),
                      name = "Significance",
                      labels = c("Not Significant", "Significant only", "Significant + Large FC")) +
  geom_label_repel(data = hspc_prog_results %>%
                     filter(external_gene_name == "Flt3"),
                   aes(label = external_gene_name), 
                   size = 3.5,
                   box.padding = 0.5,
                   show.legend = FALSE) +
  geom_point(data = hspc_prog_results %>%
               filter(external_gene_name == "Flt3"),
             size = 3,
             show.legend = FALSE) +
  labs(x = "Log2 Fold Change (Prog/HSPC)",
       y = "-log10(FDR)",
       title = "Differential Expression: HSPC vs Progenitor Cells") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 11))

# Save the volcano plot in the figures folder. 
ggsave("figures/volcano_hspc_prog.png",
       height = 5,
       width = 6.5, 
       units = "in",
       dpi = 300)




# 4. REGULAR EXPRESSION EXPLAINED. ---------------------------------------------

# Code from earlier.
pca_labelled <- pca_labelled %>%
  extract(cell_id, 
          remove = FALSE,
          c("cell_type", "cell_number"),
          "([a-zA-Z]{4})_([0-9]{3})")


# First pattern = ([a-zA-Z]{4}) :
# - Brackets - because we want to keep it and put it in cell_type.
# - [a-zA-Z] - any lower or upper case letter. 
# - Square brackets - any chars inside will be matched. 
# - {4} - 4 of them.
# So the first pattern inside the first () will match exactly 4 upper or lower
#    case letters (e.g., Prog, HSPC). 


# Second pattern = _ :
# - Matches the underscore in every cell id that separates cell type from no. 
# - Not in brackets because we don't want to keep it. 


# Third pattern = ([0-9]{3}) : 
# - [0-9] - means any number. 
# - {3} - means 3 of them. 
# - So the second pattern will match exactly 3 numbers (e.g., 001, 851).



# IMPORTANT!!
# For column names LT.HSPC, this has 6 chars with a dot - be careful here. 
# The pattern to match LT.HSPC, Prog, and HSPC = ([a-zA-Z.]{4, 6}) 







