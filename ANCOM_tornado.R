#load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr);
library(ggplot2); library(btools); library(randomcoloR); library(cowplot);
library(pairwiseAdonis); library(picante); library(gridExtra); library(plyr);
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(rstatix); library(ggpubr); library(BiodiversityR); library(ANCOMBC)
library(DT)
library(ggrepel)
library(tidyr)

# setwd
setwd("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/TE Pooling")

#Load some things
source("~/Documents/Grad School/Course Work/Bioinformatics/MergeLowAbund.R")
source("~/Documents/Grad School/Course Work/Bioinformatics/removeNARows.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/g_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/w_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/uw_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/change16STaxaNames.R")

##Data
tepools <- import_biom("SNPconfirmed_AMR_analytic_matrix.biom")
te_individuals <-import_biom("~/Documents/Research/Projects/USDA_NIFA_AMR/Sequencing/TE/SNPconfirmed_AMR_analytic_matrix_phyloseq.biom")
ARG_data <- merge_phyloseq(tepools,te_individuals)

resistome_sampledata <- import_qiime_sample_data("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/TEmetadata_resistome.txt")
sample_names(ARG_data)

pooling_resistome.ps <- merge_phyloseq(ARG_data,resistome_sampledata)

sample_names(pooling_resistome.ps)
rank_names(pooling_resistome.ps)
tax_table(pooling_resistome.ps)
colnames(tax_table(pooling_resistome.ps)) <- c("Type","Class","Mechanism","Group","Gene","SNP")

min(sample_sums(pooling_resistome.ps))
max(sample_sums(pooling_resistome.ps))
sort(sample_sums(pooling_resistome.ps)) #take out under 1000 and HP6-3

pooling_resistome_analysis.ps <- prune_samples(sample_sums(pooling_resistome.ps)>1000, pooling_resistome.ps)
pooling_resistome_analysis.ps <- prune_samples(sample_sums(pooling_resistome_analysis.ps)<32000, pooling_resistome_analysis.ps)
sort(sample_sums(pooling_resistome_analysis.ps))
median(sample_sums(pooling_resistome_analysis.ps))
min(sample_sums(pooling_resistome_analysis.ps))
max(sample_sums(pooling_resistome_analysis.ps))

sum(taxa_sums(pooling_resistome_analysis.ps)==0)
pooling_resistome_analysis.ps <- prune_taxa(taxa_sums(pooling_resistome_analysis.ps)>0,pooling_resistome_analysis.ps)

ARG_resistome_high.ps <- subset_samples(pooling_resistome_analysis.ps, PrevH=="Y")
sample_names(ARG_resistome_high.ps)

ARG_resistome_low.ps <- subset_samples(pooling_resistome_analysis.ps, Prevalence=="Low")
sample_names(ARG_resistome_low.ps)

##ANCOM-BC

ARG_resistome_high_ANCOM.ps <- ARG_resistome_high.ps
colnames(tax_table(ARG_resistome_high_ANCOM.ps)) <- c("Type","Phylum","Class","Family","Gene","SNP")

ARG_resistome_high_ANCOM.ps <- prune_taxa(taxa_sums(ARG_resistome_high_ANCOM.ps)>0,ARG_resistome_high_ANCOM.ps)

sample_data(ARG_resistome_high_ANCOM.ps)$StatsPool <- c("A","D","A","D","A","F","C","E","B","B","E","C","F","E","F","D","C","D","C","C","C","D","D","C","D","C","D","C")

ancom_resistome_high = ancombc2(data = ARG_resistome_high_ANCOM.ps, assay_name = "counts", tax_level = "Family", 
                                  fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                                  p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                  group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                                  global = TRUE, pairwise = TRUE,
                                  dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
dunn_res_resistome_high <- ancom_resistome_high$res_dunn

# Check the structure of the pairwise results
str(dunn_res_resistome_high)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
dunn_res_resistome_high = dunn_res_resistome_high %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

df_fig_poolsize_dunn_resistome_high <- dunn_res_resistome_high %>%
  # Step 1: Apply the condition for the lfc_timepoint_num columns
  dplyr::mutate(across(starts_with("lfc_Number.in.Pool"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x,.x)),
                # Step 2: Round the lfc_timepoint_num columns and create _rounded columns
                across(starts_with("lfc_Number.in.Pool"), ~round(.x, 3), .names = "{.col}_rounded"),
                # Step 3: Assign colors based on diff_timepoint_num and passed_ss columns
                across(starts_with("diff_Number.in.Pool"), 
                       ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "darkred", "darkgrey"), 
                       .names = "{.col}_color")) %>%
  # Step 4: Pivot the _rounded columns into long format
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  # Step 5: Pivot the _color columns into long format
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  # Step 6: Filter to ensure the group and color_group match
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  # Step 7: Select the relevant columns and arrange by taxon
  select(-color_group) %>%
  arrange(taxon)

# Check which levels are available for comparison
as.factor(df_fig_poolsize_dunn_resistome_high$group)

df_fig_poolsize_dunn_resistome_high_filtered <- df_fig_poolsize_dunn_resistome_high %>%
  dplyr::filter(str_ends(group, "_rounded")) %>%
  dplyr::mutate(group = dplyr::case_when(
    group == "Number.in.Pool3_rounded" ~ "Pools of 3 vs Individuals",
    group == "Number.in.Pool6_rounded" ~ "Pools of 6 vs Individuals",
    group == "Number.in.Pool12_rounded" ~ "Pools of 12 vs Individuals",
    TRUE ~ group)) %>%
  dplyr::filter(group %in% c("Pools of 3 vs Individuals", 
                             "Pools of 6 vs Individuals", 
                             "Pools of 12 vs Individuals")) %>%
  dplyr::mutate(group = factor(group, levels = c("Pools of 3 vs Individuals", 
                                                 "Pools of 6 vs Individuals", 
                                                 "Pools of 12 vs Individuals"))) %>%
  droplevels()




levels(as.factor(df_fig_poolsize_dunn_resistome_high_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_resistome_high <- ancom_resistome_high$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_resistome_high[is.na(log_table_resistome_high)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_resistome_high <- exp(log_table_resistome_high)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_resistome_high <- apply(pseudo_counts_resistome_high, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_resistome_high <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_resistome_high <- t(t(pseudo_counts_resistome_high) * (desired_total_count_resistome_high / total_counts_resistome_high))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_resistome_high.df <- as.data.frame(scaled_counts_resistome_high)

# Add a column for taxa names
scaled_counts_resistome_high.df$taxon <- rownames(scaled_counts_resistome_high.df)

# Convert from wide to long format using pivot_longer
scaled_counts_resistome_high.long <- scaled_counts_resistome_high.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
ARG_resistome_high_ANCOM.melt <- psmelt(ARG_resistome_high_ANCOM.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_ARG_high.melt <- ARG_resistome_high_ANCOM.melt %>%
  left_join(scaled_counts_resistome_high.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Group_avg_counts_ARG_high_by_poolsize <- BiasAdj_ARG_high.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunn_resistome_high_filtered_merged <- df_fig_poolsize_dunn_resistome_high_filtered %>%
  left_join(Group_avg_counts_ARG_high_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_ARG_high <- as.data.frame(tax_table(ARG_resistome_high_ANCOM.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_ARG_high <- df_fig_poolsize_dunn_resistome_high_filtered_merged %>%
  left_join(tax_data_ARG_high, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_ARG_high <- df_with_taxonomy_ARG_high %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_ARG_high$group = factor(df_with_taxonomy_ARG_high$group, levels = c("Pools of 3 vs Individuals", 
                                                                                 "Pools of 6 vs Individuals", 
                                                                                 "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_ARG_high <- df_with_taxonomy_ARG_high %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_ARG_high <- ggplot(df_with_taxonomy_ARG_high, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_ARG_high %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 3,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Mechanism", title = " logFC of Gene Group by Pool Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4, scales = "fixed")

main_species_plot_ARG_high

# Create the taxonomy plot data
taxonomy_plot_data_ARG_high <- df_with_taxonomy_ARG_high %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_ARG_high <- taxonomy_plot_data_ARG_high %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_ARG_high <- ggplot(taxonomy_plot_data_ARG_high) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_ARG_high$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_ARG_high  <- plot_grid(
  taxonomy_plot_ARG_high  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_ARG_high  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_ARG_high
# Save the final plot
ggsave("{work on file name}", plot = combined_plot_FECAL_NEG, width = 20, height = 11, dpi = 300)

##Low
##ANCOM-BC
ARG_resistome_low_ANCOM.ps <- ARG_resistome_low.ps
colnames(tax_table(ARG_resistome_low_ANCOM.ps)) <- c("Type","Phylum","Class","Family","Gene","SNP")

sample_data(ARG_resistome_low_ANCOM.ps)$StatsPool <- c("F","D","E","B","D","D","F","A","C","E","B","C","B","C","E","A","F")

ancom_resistome_low = ancombc2(data = ARG_resistome_low_ANCOM.ps, assay_name = "counts", tax_level = "Family", 
                                fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                alpha = 0.05, n_cl = 2, verbose = TRUE,
                                global = TRUE, pairwise = TRUE,
                                dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
pair_res_resistome_low <- ancom_resistome_low$res_pair

# Check the structure of the pairwise results
str(pair_res_resistome_low)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
pair_res_resistome_low = pair_res_resistome_low %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

unique(pair_res_resistome_low$diff_Number.in.Pool6)
unique(pair_res_resistome_low$diff_Number.in.Pool12)
unique(pair_res_resistome_low$diff_Number.in.Pool12_Number.in.Pool6)


df_fig_poolsize_pair_resistome_low <- pair_res_resistome_low %>%
  # Step 1: Apply the condition for the lfc_timepoint_num columns
  dplyr::mutate(across(starts_with("lfc_Number.in.Pool"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x,.x)),
                # Step 2: Round the lfc_timepoint_num columns and create _rounded columns
                across(starts_with("lfc_Number.in.Pool"), ~round(.x, 3), .names = "{.col}_rounded"),
                # Step 3: Assign colors based on diff_timepoint_num and passed_ss columns
                across(starts_with("diff_Number.in.Pool"), 
                       ~ifelse(.x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1, "darkred", "darkgrey"), 
                       .names = "{.col}_color")) %>%
  # Step 4: Pivot the _rounded columns into long format
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  # Step 5: Pivot the _color columns into long format
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  # Step 6: Filter to ensure the group and color_group match
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  # Step 7: Select the relevant columns and arrange by taxon
  select(-color_group) %>%
  arrange(taxon)

# Check which levels are available for comparison
as.factor(df_fig_poolsize_dunn_resistome_low$group)

df_fig_poolsize_dunn_resistome_low_filtered <- df_fig_poolsize_dunn_resistome_low %>%
  dplyr::filter(str_ends(group, "_rounded")) %>%
  dplyr::mutate(group = dplyr::case_when(
    group == "Number.in.Pool3_rounded" ~ "Pools of 3 vs Individuals",
    group == "Number.in.Pool6_rounded" ~ "Pools of 6 vs Individuals",
    group == "Number.in.Pool12_rounded" ~ "Pools of 12 vs Individuals",
    TRUE ~ group)) %>%
  dplyr::filter(group %in% c("Pools of 3 vs Individuals", 
                             "Pools of 6 vs Individuals", 
                             "Pools of 12 vs Individuals")) %>%
  dplyr::mutate(group = factor(group, levels = c("Pools of 3 vs Individuals", 
                                                 "Pools of 6 vs Individuals", 
                                                 "Pools of 12 vs Individuals"))) %>%
  droplevels()




levels(as.factor(df_fig_poolsize_dunn_resistome_low_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_resistome_low <- ancom_resistome_low$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_resistome_low[is.na(log_table_resistome_low)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_resistome_low <- exp(log_table_resistome_low)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_resistome_low <- apply(pseudo_counts_resistome_low, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_resistome_low <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_resistome_low <- t(t(pseudo_counts_resistome_low) * (desired_total_count_resistome_low / total_counts_resistome_low))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_resistome_low.df <- as.data.frame(scaled_counts_resistome_low)

# Add a column for taxa names
scaled_counts_resistome_low.df$taxon <- rownames(scaled_counts_resistome_low.df)

# Convert from wide to long format using pivot_longer
scaled_counts_resistome_low.long <- scaled_counts_resistome_low.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
ARG_resistome_low_ANCOM.melt <- psmelt(ARG_resistome_low_ANCOM.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_ARG_low.melt <- ARG_resistome_low_ANCOM.melt %>%
  left_join(scaled_counts_resistome_low.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Group_avg_counts_ARG_low_by_poolsize <- BiasAdj_ARG_low.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunn_resistome_low_filtered_merged <- df_fig_poolsize_dunn_resistome_low_filtered %>%
  left_join(Group_avg_counts_ARG_low_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_ARG_low <- as.data.frame(tax_table(ARG_resistome_low_ANCOM.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_ARG_low <- df_fig_poolsize_dunn_resistome_low_filtered_merged %>%
  left_join(tax_data_ARG_low, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_ARG_low <- df_with_taxonomy_ARG_low %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_ARG_low$group = factor(df_with_taxonomy_ARG_low$group, levels = c("Pools of 3 vs Individuals", 
                                                                                     "Pools of 6 vs Individuals", 
                                                                                     "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_ARG_low <- df_with_taxonomy_ARG_low %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_ARG_low <- ggplot(df_with_taxonomy_ARG_low, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_ARG_low %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 3,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Mechanism", title = " logFC of Gene Group by Pool Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4, scales = "fixed")

main_species_plot_ARG_low

# Create the taxonomy plot data
taxonomy_plot_data_ARG_low <- df_with_taxonomy_ARG_low %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_ARG_low <- taxonomy_plot_data_ARG_low %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_ARG_low <- ggplot(taxonomy_plot_data_ARG_low) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_ARG_low$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_ARG_low  <- plot_grid(
  taxonomy_plot_ARG_low  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_ARG_low  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_ARG_low
# Save the final plot
ggsave("{work on file name}", plot = combined_plot_FECAL_NEG, width = 20, height = 11, dpi = 300)

