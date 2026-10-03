##ANCOMBC-2_tornado dunnet

##HS_16S
df_dunnet_16S_HS = dunnet_16S_HS %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_fig_poolsize_dunnet_16S_HS <- df_dunnet_16S_HS %>%
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
as.factor(df_fig_poolsize_dunnet_16S_HS$group)

df_fig_poolsize_dunnet_16S_HS_filtered <- df_fig_poolsize_dunnet_16S_HS %>%
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




levels(as.factor(df_fig_poolsize_dunnet_16S_HS_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_HS <- ancom_Family_16S_HS$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_HS[is.na(log_table_16S_HS)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_HS <- exp(log_table_16S_HS)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_HS <- apply(pseudo_counts_16S_HS, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_HS <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_HS <- t(t(pseudo_counts_16S_HS) * (desired_total_count_16S_HS / total_counts_16S_HS))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_HS.df <- as.data.frame(scaled_counts_16S_HS)

# Add a column for taxa names
scaled_counts_16S_HS.df$taxon <- rownames(scaled_counts_16S_HS.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_HS.long <- scaled_counts_16S_HS.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
#HS.micro_data.melt <- psmelt(HSpool_16S.ps)
## try glom to genus
HSpool_16S_family.ps <- tax_glom(HSpool_16S.ps, taxrank = "Family")
HS.micro_data.melt <- psmelt(HSpool_16S_family.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_HS.micro_data.melt <- HS.micro_data.melt %>%
  left_join(scaled_counts_16S_HS.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_HS_by_poolsize <- BiasAdj_HS.micro_data.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = F)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_HS_filtered_merged_data <- df_fig_poolsize_dunnet_16S_HS_filtered %>%
  left_join(Family_avg_counts_16S_HS_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_HS <- as.data.frame(tax_table(HSpool_16S.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_HS <- df_fig_poolsize_dunnet_16S_HS_filtered_merged_data %>%
  left_join(tax_data_16S_HS, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_HS <- df_with_taxonomy_16S_HS %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_HS$group = factor(df_with_taxonomy_16S_HS$group, levels = c("Pools of 3 vs Individuals", 
                                                                                 "Pools of 6 vs Individuals", 
                                                                                 "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_HS <- df_with_taxonomy_16S_HS %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_HS <- ggplot(na.omit(df_with_taxonomy_16S_HS), aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_HS %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 8,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 22, colour = "black"),
    axis.text.x = element_text(size = 16, angle = 45, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4)

main_species_plot_16S_HS
ggsave("Results/16S/Figures/DA_micro_family_HS_simple.tiff", plot = main_species_plot_16S_HS, width = 20, height = 11, dpi = 300)



# Create the taxonomy plot data
taxonomy_plot_data_16S_HS <- df_with_taxonomy_16S_HS %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_HS <- taxonomy_plot_data_16S_HS %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_HS <- ggplot(taxonomy_plot_data_16S_HS) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_HS$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )

ggsave(filename = "Results/16S/Figures/DA_taxonomy_HS.tiff", taxonomy_plot_16S_HS,dpi = 300)


combined_plot_16S_HS  <- plot_grid(
  taxonomy_plot_16S_HS  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_HS  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_HS

ggsave("Results/16S/Figures/DA_micro_family_HS.jpg", plot = combined_plot_16S_HS, width = 20, height = 11, dpi = 300)



##HR_16S
df_dunnet_16S_HR = dunnet_16S_HR %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_fig_poolsize_dunnet_16S_HR <- df_dunnet_16S_HR %>%
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
as.factor(df_fig_poolsize_dunnet_16S_HR$group)

df_fig_poolsize_dunnet_16S_HR_filtered <- df_fig_poolsize_dunnet_16S_HR %>%
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




levels(as.factor(df_fig_poolsize_dunnet_16S_HR_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_HR <- ancom_Family_16S_HR$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_HR[is.na(log_table_16S_HR)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_HR <- exp(log_table_16S_HR)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_HR <- apply(pseudo_counts_16S_HR, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_HR <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_HR <- t(t(pseudo_counts_16S_HR) * (desired_total_count_16S_HR / total_counts_16S_HR))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_HR.df <- as.data.frame(scaled_counts_16S_HR)

# Add a column for taxa names
scaled_counts_16S_HR.df$taxon <- rownames(scaled_counts_16S_HR.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_HR.long <- scaled_counts_16S_HR.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
HRpool_16S_Family.ps <- tax_glom(HRpool_16S.ps, taxrank = "Family")
HR.micro_data.melt <- psmelt(HRpool_16S_Family.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_HR.micro_data.melt <- HR.micro_data.melt %>%
  left_join(scaled_counts_16S_HR.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_HR_by_poolsize <- BiasAdj_HR.micro_data.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_HR_filtered_merged_data <- df_fig_poolsize_dunnet_16S_HR_filtered %>%
  left_join(Family_avg_counts_16S_HR_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_HR <- as.data.frame(tax_table(HSpool_16S.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_HR <- df_fig_poolsize_dunnet_16S_HR_filtered_merged_data %>%
  left_join(tax_data_16S_HR, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_HR <- df_with_taxonomy_16S_HR %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_HR$group = factor(df_with_taxonomy_16S_HR$group, levels = c("Pools of 3 vs Individuals", 
                                                                                 "Pools of 6 vs Individuals", 
                                                                                 "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_HR <- df_with_taxonomy_16S_HR %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_HR <- ggplot(na.omit(df_with_taxonomy_16S_HR), aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_HR %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 8,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 22, colour = "black"),
    axis.text.x = element_text(size = 20, angle = 45, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4, scales = "fixed")

main_species_plot_16S_HR
ggsave(filename = "Results/16S/Figures/DA_micro_family_HR_simple.tiff",plot = main_species_plot_16S_HR, dpi = 300, width = 20, height = 11)

# Create the taxonomy plot data
taxonomy_plot_data_16S_HR <- df_with_taxonomy_16S_HR %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_HR <- taxonomy_plot_data_16S_HR %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_HR <- ggplot(taxonomy_plot_data_16S_HR) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 4.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_HR$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )

ggsave(filename = "Results/16S/Figures/DA_taxonomy_HR.tiff",plot = taxonomy_plot_16S_HR,dpi = 300)

combined_plot_16S_HR  <- plot_grid(
  taxonomy_plot_16S_HR  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_HR  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_HR
ggsave("Results/16S/Figures/DA_micro_family_HR.jpg", plot = combined_plot_16S_HR, width = 20, height = 11, dpi = 300)



##LS_16S
df_dunnet_16S_LS = dunnet_16S_LS %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_fig_poolsize_dunnet_16S_LS <- df_dunnet_16S_LS %>%
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
as.factor(df_fig_poolsize_dunnet_16S_LS$group)

df_fig_poolsize_dunnet_16S_LS_filtered <- df_fig_poolsize_dunnet_16S_LS %>%
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




levels(as.factor(df_fig_poolsize_dunnet_16S_LS_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_LS <- ancom_Family_16S_LS$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_LS[is.na(log_table_16S_LS)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_LS <- exp(log_table_16S_LS)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_LS <- apply(pseudo_counts_16S_LS, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_LS <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_LS <- t(t(pseudo_counts_16S_LS) * (desired_total_count_16S_LS / total_counts_16S_LS))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_LS.df <- as.data.frame(scaled_counts_16S_LS)

# Add a column for taxa names
scaled_counts_16S_LS.df$taxon <- rownames(scaled_counts_16S_LS.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_LS.long <- scaled_counts_16S_LS.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
#LS.micro_data.melt <- psmelt(LSpool_16S.ps)
## try glom to Family
LSpool_16S_family.ps <- tax_glom(LSpool_16S.ps, taxrank = "Family")
LS.micro_data.melt <- psmelt(LSpool_16S_family.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_LS.micro_data.melt <- LS.micro_data.melt %>%
  left_join(scaled_counts_16S_LS.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_LS_by_poolsize <- BiasAdj_LS.micro_data.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = F)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_LS_filtered_merged_data <- df_fig_poolsize_dunnet_16S_LS_filtered %>%
  left_join(Family_avg_counts_16S_LS_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_LS <- as.data.frame(tax_table(LSpool_16S.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_LS <- df_fig_poolsize_dunnet_16S_LS_filtered_merged_data %>%
  left_join(tax_data_16S_LS, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_LS <- df_with_taxonomy_16S_LS %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_LS$group = factor(df_with_taxonomy_16S_LS$group, levels = c("Pools of 3 vs Individuals", 
                                                                                 "Pools of 6 vs Individuals", 
                                                                                 "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_LS <- df_with_taxonomy_16S_LS %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_LS <- ggplot(na.omit(df_with_taxonomy_16S_LS), aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_LS %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 8,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 22, colour = "black"),
    axis.text.x = element_text(size = 18, angle = 45, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4)

main_species_plot_16S_LS
ggsave(filename = "Results/16S/Figures/DA_micro_family_LS_simple.tiff",plot = main_species_plot_16S_LS, dpi = 300, height = 11, width = 20)


# Create the taxonomy plot data
taxonomy_plot_data_16S_LS <- df_with_taxonomy_16S_LS %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_LS <- taxonomy_plot_data_16S_LS %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_LS <- ggplot(taxonomy_plot_data_16S_LS) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 4.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_LS$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )

ggsave(filename = "Results/16S/Figures/DA_taxonomy_LS.tiff", plot=taxonomy_plot_16S_LS,dpi = 300)

combined_plot_16S_LS  <- plot_grid(
  taxonomy_plot_16S_LS  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_LS  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_LS
ggsave("Results/16S/Figures/DA_micro_family_LS.jpg", plot = combined_plot_16S_LS, width = 20, height = 11, dpi = 300)



##LR_16S
df_dunnet_16S_LR = dunnet_16S_LR %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_fig_poolsize_dunnet_16S_LR <- df_dunnet_16S_LR %>%
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
as.factor(df_fig_poolsize_dunnet_16S_LR$group)

df_fig_poolsize_dunnet_16S_LR_filtered <- df_fig_poolsize_dunnet_16S_LR %>%
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




levels(as.factor(df_fig_poolsize_dunnet_16S_LR_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_LR <- ancom_Family_16S_LR$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_LR[is.na(log_table_16S_LR)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_LR <- exp(log_table_16S_LR)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_LR <- apply(pseudo_counts_16S_LR, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_LR <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_LR <- t(t(pseudo_counts_16S_LR) * (desired_total_count_16S_LR / total_counts_16S_LR))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_LR.df <- as.data.frame(scaled_counts_16S_LR)

# Add a column for taxa names
scaled_counts_16S_LR.df$taxon <- rownames(scaled_counts_16S_LR.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_LR.long <- scaled_counts_16S_LR.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
LRpool_16S_Family.ps <- tax_glom(LRpool_16S.ps, taxrank = "Family")
LR.micro_data.melt <- psmelt(LRpool_16S_Family.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_LR.micro_data.melt <- LR.micro_data.melt %>%
  left_join(scaled_counts_16S_LR.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_LR_by_poolsize <- BiasAdj_LR.micro_data.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_LR_filtered_merged_data <- df_fig_poolsize_dunnet_16S_LR_filtered %>%
  left_join(Family_avg_counts_16S_LR_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_LR <- as.data.frame(tax_table(HSpool_16S.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_LR <- df_fig_poolsize_dunnet_16S_LR_filtered_merged_data %>%
  left_join(tax_data_16S_LR, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_LR <- df_with_taxonomy_16S_LR %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_LR$group = factor(df_with_taxonomy_16S_LR$group, levels = c("Pools of 3 vs Individuals", 
                                                                                 "Pools of 6 vs Individuals", 
                                                                                 "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_LR <- df_with_taxonomy_16S_LR %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_LR <- ggplot(na.omit(df_with_taxonomy_16S_LR), aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_LR %>% dplyr::filter(color == "darkred"),
                  aes(label = taxon),
                  color = "darkred",
                  size = 6,                  # Adjust label size as needed
                  position = position_dodge(width = 0.5),
                  box.padding = 0.4,         # Moderate padding to allow some vertical movement
                  point.padding = 0.4,       # Padding around points for some horizontal spacing
                  max.overlaps = Inf,        # Allow ggrepel to try placing all labels
                  segment.color = "grey70") +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 4, scales = "fixed")

main_species_plot_16S_LR
ggsave("Results/16S/Figures/DA_micro_family_LR_simple.tiff", plot = main_species_plot_16S_LR, width = 20, height = 11, dpi = 300)



# Create the taxonomy plot data
taxonomy_plot_data_16S_LR <- df_with_taxonomy_16S_LR %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_LR <- taxonomy_plot_data_16S_LR %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_LR <- ggplot(taxonomy_plot_data_16S_LR) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_LR$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )

ggsave(filename = "Results/16S/Figures/DA_taxonomy_LR.tiff", taxonomy_plot_16S_LR,dpi = 300)

combined_plot_16S_LR  <- plot_grid(
  taxonomy_plot_16S_LR  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_LR  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,3),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_LR
ggsave("Results/16S/Figures/DA_micro_family_LR.jpg", plot = combined_plot_16S_LR, width = 20, height = 11, dpi = 300)

