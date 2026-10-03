####Core microbiome analysis

##packages
library(microbiome)


##data
HSpool_16S.ps
HSpool_16S_individual.ps <- subset_samples(HSpool_16S.ps,Number.in.Pool=="1")
sum(taxa_sums(HSpool_16S_individual.ps)==0)
HSpool_16S_individual.ps <- prune_taxa(taxa_sums(HSpool_16S_individual.ps)>0,HSpool_16S_individual.ps)
sum(sample_sums(HSpool_16S_individual.ps)==0)

HSpool_16S_individual.css <- phyloseq_transform_css(HSpool_16S_individual.ps,log = F)
HSpool_16S_individual.ra <- transform_sample_counts(HSpool_16S_individual.css, function(x) {x/sum(x)}*100)

##now get core-done with ra data
core_members(HSpool_16S_individual.ra, detection = 0, prevalence = 50/100) #gives otus, but there are 46
HSpool_16S_individual.core <- core(HSpool_16S_individual.ps, detection = 0, prevalence = 0.5)

taxa(HSpool_16S_individual.core)

tax_table(HSpool_16S_individual.core)

##what if we do on "raw"/css data?
core_members(HSpool_16S_individual.ps, detection = 0, prevalence = 0.5) #75
core_members(HSpool_16S_individual.css, detection = 0, prevalence = 0.5) #75

###Going on with ps, then will normalize...that ok?
HSpool_16S_3.ps <- subset_samples(HSpool_16S.ps,Number.in.Pool=="3")
sum(taxa_sums(HSpool_16S_3.ps)==0)
HSpool_16S_3.ps <- prune_taxa(taxa_sums(HSpool_16S_3.ps)>0,HSpool_16S_3.ps)
sum(sample_sums(HSpool_16S_3.ps)==0)

core_members(HSpool_16S_3.ps, detection=0, prevalence = 0.5)
HSpool_16S_3.core <- core(HSpool_16S_3.ps, detection = 0, prevalence = 0.5)

HSpool_16S_6.ps <- subset_samples(HSpool_16S.ps, Number.in.Pool=="6")
HSpool_16S_6.ps <- prune_taxa(taxa_sums(HSpool_16S_6.ps)>0,HSpool_16S_6.ps)
core_members(HSpool_16S_6.ps,detection = 0, prevalence = 0.5)
HSpool_16S_6.core <- core(HSpool_16S_6.ps, detection = 0, prevalence = 0.5)

HSpool_16S_12.ps <- subset_samples(HSpool_16S.ps, Number.in.Pool=="12")
HSpool_16S_12.ps <- prune_taxa(taxa_sums(HSpool_16S_12.ps)>0,HSpool_16S_12.ps)
core_members(HSpool_16S_12.ps,detection = 0, prevalence = 0.5)
HSpool_16S_12.core <- core(HSpool_16S_12.ps, detection = 0, prevalence = 0.5)


otu_table_core_1 <- otu_table(HSpool_16S_individual.core,taxa_are_rows = T)
otu_table_core_3 <- otu_table(HSpool_16S_3.core, taxa_are_rows = T)
otu_table_core_6 <- otu_table(HSpool_16S_6.core, taxa_are_rows = T)
otu_table_core_12 <- otu_table(HSpool_16S_12.core, taxa_are_rows = T)

tax_table_core_1 <- tax_table(HSpool_16S_individual.core)
tax_table_core_3 <- tax_table(HSpool_16S_3.core)
tax_table_core_6 <- tax_table(HSpool_16S_6.core)
tax_table_core_12 <- tax_table(HSpool_16S_12.core)


otu_table_core <- merge_phyloseq(otu_table_core_1,otu_table_core_3,otu_table_core_6,otu_table_core_12)
tax_table_core <- merge_phyloseq(tax_table_core_1,tax_table_core_3,tax_table_core_6,tax_table_core_12)

core_microbiome_high <- merge_phyloseq(otu_table_core,tax_table_core)

HSpool_16S_sampledata <- sample_data(HSpool_16S.ps)
sample_names(HSpool_16S_sampledata)

HSpool_16S_core.ps <- merge_phyloseq(core_microbiome_high,HSpool_16S_sampledata)
rank_names(HSpool_16S_core.ps)


##let's do unweighted unifrac
#double check no taxa not found in any sample
sum(taxa_sums(HSpool_16S_core.ps)==0)
sum(sample_sums(HSpool_16S_core.ps)==0)

tree_16S <- read_tree("../DATA/16S/tree.nwk")

HSpool_16S_core.ps <- merge_phyloseq(HSpool_16S_core.ps,tree_16S)

HSpool_16S_core.css <- phyloseq_transform_css(HSpool_16S_core.ps,log = F)
HSpool_16S_core.css.df <- as(sample_data(HSpool_16S_core.css),"data.frame")
HSpool_16S_core_uwunifrac.dist <- uwunifrac(HSpool_16S_core.css)
HSpool_16S_core_uwunifrac.ord <- ordinate(HSpool_16S_core.css,method = "NMDS",distance = HS_16S_core_uwunifrac.dist)
#HSpool_16S_core_uwunifrac.pcoa <- ordinate(HSpool_16S_core.css, method = "PCoA", distance = HS_16S_core_uwunifrac.dist)

HSpool_16S_core_uwunifrac_nmdsplot <- ordiplot(HS_16S_core_uwunifrac.ord$points)
HSpool_16S_core_uwunifrac_pcoaplot <- ordiplot(HS_16S_core_uwunifrac.pcoa$vectors)
HSpool_16S_core_uwunifrac_siteslong_nmds <-sites.long(HS_16S_core_uwunifrac_nmdsplot,HSpool_16S_core.css.df)
HSpool_16S_core_uwunifrac_siteslong_pcoa <- sites.long(HS_16S_core_uwunifrac_pcoaplot, HSpool_16S_core.css.df)


HSpool_16S_core_uwunifrac_centroids_nmds <- envfit(HS_16S_core_uwunifrac.ord~HSpool_16S_core.css.df$Number.in.Pool)
HSpool_16S_core_uwunifrac_centroids_nmds
#HSpool_16S_core_uwunifrac_centroids_pcoa <- envfit(HS_16S_core_uwunifrac.pcoa$vectors~HSpool_16S_core.css.df$Number.in.Pool)
#HSpool_16S_core_uwunifrac_centroids_pcoa


HSpool_16S_core_uwunifrac_NMDS_col1 <- c(-0.0765,0.3823,0.1933,0.3613)
HSpool_16S_core_uwunifrac_NMDS_col2 <- c(0.0046,0.0552,-0.0591,-0.0602)
HSpool_16S_core_uwunifrac_centroids.df <-data.frame(pool_col,HSpool_16S_core_uwunifrac_NMDS_col1,HSpool_16S_core_uwunifrac_NMDS_col2)

HSpool_16S_core_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Unweighted Unifrac") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = HSpool_16S_core_uwunifrac_siteslong_nmds, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2)) +
  stat_ellipse(data = HSpool_16S_core_uwunifrac_siteslong_nmds, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = HSpool_16S_core_uwunifrac_centroids.df, aes(x=HSpool_16S_core_uwunifrac_NMDS_col1, y=HSpool_16S_core_uwunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = HSpool_16S_core_uwunifrac_centroids.df, aes(x=HSpool_16S_core_uwunifrac_NMDS_col1, y=HSpool_16S_core_uwunifrac_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        plot.title=element_blank())

plot_ordination(HSpool_16S_core.css, subseth_wunifrac.pcoa,type = "samples", color = "Number.in.Pool", shape = "Number.in.Pool") +
  theme_bw() +
  #labs(x= "NMDS1", y= "NMDS2") +
  scale_shape_manual(values=c(18,19,19,19))+
  stat_ellipse(geom = "polygon", aes(fill= Number.in.Pool), alpha = 0.5, lty = 2, size = 1, level =0.90) +
  geom_point(data = subseth_wunifrac_PCoA_centroids.df, aes(x=subseth_wunifrac_PCoA_col1, y=subseth_wunifrac_PCoA_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subseth_wunifrac_PCoA_centroids.df, aes(x=subseth_wunifrac_PCoA_col1, y=subseth_wunifrac_PCoA_col2, label = pool_col), colour = "white", size = c(3,3,3,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  #scale_y_continuous(limits = c(-0.17,0.2)) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour= "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))

ggsave("hs_wunifrac.tiff", plot = hs_wunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")

HSpool_16S_core.css.df$Number.in.Pool <- as.factor(HSpool_16S_core.css.df$Number.in.Pool)

adonis2(HSpool_16S_core_uwunifrac.dist ~ Number.in.Pool, HSpool_16S_core.css.df, nperm = 9999)

HSpool_16S_core_uwunifrac_pair <- pairwise.adonis2(HSpool_16S_core_uwunifrac.dist ~ Number.in.Pool, HSpool_16S_core.css.df, nperm = 9999, p.adjust.methods="BH") # NS

write.csv(HSpool_16S_core_uwunifrac_pair,"Results/16S/Stats/core_bd_uwunifrac_pairadonis.csv")


HSpool_16S_core_uwunifrac.disper <- betadisper(HSpool_16S_core_uwunifrac.dist, HSpool_16S_core.css.df$Number.in.Pool)
plot(HSpool_16S_core_uwunifrac.disper)
boxplot(HSpool_16S_core_uwunifrac.disper)
TukeyHSD(HSpool_16S_core_uwunifrac.disper)
anova(HSpool_16S_core_uwunifrac.disper)
permutest(HSpool_16S_core_uwunifrac.disper, permutations = 9999, pairwise=F)
HSpool_16S_core_uwunifrac.permdisp <- permutest(HSpool_16S_core_uwunifrac.disper, permutations = 9999, pairwise = T)
HSpool_16S_core_uwunifrac.permdisp #
write.csv(HSpool_16S_core_uwunifrac.permdisp[["pairwise"]][["permuted"]],"Results/16S/Stats/core_bd_uwunifrac_permdisp.csv")

#Doing weighted unifrac, because unweighted was interesting
HSpool_16S_core_wunifrac.dist <- wunifrac(HSpool_16S_core.css)
HSpool_16S_core_wunifrac.ord <- ordinate(HSpool_16S_core.css,method = "NMDS",distance = HSpool_16S_core_wunifrac.dist)
#HSpool_16S_core_wunifrac.pcoa <- ordinate(HSpool_16S_core.css, method = "PCoA", distance = HS_16S_core_wunifrac.dist)

HSpool_16S_core_wunifrac_nmdsplot <- ordiplot(HSpool_16S_core_wunifrac.ord$points)
#HSpool_16S_core_wunifrac_pcoaplot <- ordiplot(HS_16S_core_wunifrac.pcoa$vectors)
HSpool_16S_core_wunifrac_siteslong_nmds <-sites.long(HSpool_16S_core_wunifrac_nmdsplot,HSpool_16S_core.css.df)
#HSpool_16S_core_wunifrac_siteslong_pcoa <- sites.long(HS_16S_core_wunifrac_pcoaplot, HSpool_16S_core.css.df)


HSpool_16S_core_wunifrac_centroids_nmds <- envfit(HSpool_16S_core_wunifrac.ord~HSpool_16S_core.css.df$Number.in.Pool)
HSpool_16S_core_wunifrac_centroids_nmds
#HSpool_16S_core_wunifrac_centroids_pcoa <- envfit(HS_16S_core_wunifrac.pcoa$vectors~HSpool_16S_core.css.df$Number.in.Pool)
#HSpool_16S_core_wunifrac_centroids_pcoa


HSpool_16S_core_wunifrac_NMDS_col1 <- c(-0.0026,0.0003,0.0331,0.0027)
HSpool_16S_core_wunifrac_NMDS_col2 <- c(-0.0024,0.0097,0.0229,-0.0002)
HSpool_16S_core_wunifrac_centroids.df <-data.frame(pool_col,HSpool_16S_core_wunifrac_NMDS_col1,HSpool_16S_core_wunifrac_NMDS_col2)

HSpool_16S_core_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Weighted Unifrac") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = HSpool_16S_core_wunifrac_siteslong_nmds, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2)) +
  stat_ellipse(data = HSpool_16S_core_wunifrac_siteslong_nmds, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = HSpool_16S_core_wunifrac_centroids.df, aes(x=HSpool_16S_core_wunifrac_NMDS_col1, y=HSpool_16S_core_wunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = HSpool_16S_core_wunifrac_centroids.df, aes(x=HSpool_16S_core_wunifrac_NMDS_col1, y=HSpool_16S_core_wunifrac_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title.x = element_blank(),
        plot.title=element_blank())

plot_ordination(HSpool_16S_core.css, subseth_wunifrac.pcoa,type = "samples", color = "Number.in.Pool", shape = "Number.in.Pool") +
  theme_bw() +
  #labs(x= "NMDS1", y= "NMDS2") +
  scale_shape_manual(values=c(18,19,19,19))+
  stat_ellipse(geom = "polygon", aes(fill= Number.in.Pool), alpha = 0.5, lty = 2, size = 1, level =0.90) +
  geom_point(data = subseth_wunifrac_PCoA_centroids.df, aes(x=subseth_wunifrac_PCoA_col1, y=subseth_wunifrac_PCoA_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subseth_wunifrac_PCoA_centroids.df, aes(x=subseth_wunifrac_PCoA_col1, y=subseth_wunifrac_PCoA_col2, label = pool_col), colour = "white", size = c(3,3,3,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  #scale_y_continuous(limits = c(-0.17,0.2)) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour= "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))

ggsave("hs_wunifrac.tiff", plot = hs_wunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")

adonis2(HSpool_16S_core_wunifrac.dist ~ Number.in.Pool, HSpool_16S_core.css.df, nperm = 9999)

HSpool_16S_core_wunifrac_pair <- pairwise.adonis2(HSpool_16S_core_wunifrac.dist ~ Number.in.Pool, HSpool_16S_core.css.df, nperm = 9999, p.adjust.methods="BH") # NS

write.csv(HSpool_16S_core_wunifrac_pair,"Results/16S/Stats/core_bd_wunifrac_pairadonis.csv")


HSpool_16S_core_wunifrac.disper <- betadisper(HSpool_16S_core_wunifrac.dist, HSpool_16S_core.css.df$Number.in.Pool)
plot(HSpool_16S_core_wunifrac.disper)
boxplot(HSpool_16S_core_wunifrac.disper)
TukeyHSD(HSpool_16S_core_wunifrac.disper)
anova(HSpool_16S_core_wunifrac.disper)
permutest(HSpool_16S_core_wunifrac.disper, permutations = 9999, pairwise=F)
HSpool_16S_core_wunifrac.permdisp <- permutest(HSpool_16S_core_wunifrac.disper, permutations = 9999, pairwise = T)
HSpool_16S_core_wunifrac.permdisp #
write.csv(HSpool_16S_core_wunifrac.permdisp[["pairwise"]][["permuted"]],"Results/16S/Stats/core_bd_wunifrac_permdisp.csv")


##Needs phylogenetic tree
sample_data(HSpool_16S_core.ps)
###ANCOMBC
HS_16S_core_ancom_Family <- ancombc2(data = HSpool_16S_core.ps, assay_name = "counts", tax_level = "Family", 
         fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
         p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
         group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
         alpha = 0.05, n_cl = 2, verbose = TRUE,
         global = TRUE, pairwise = TRUE,
         dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
pairwise_ancom_HS_16S_core <- HS_16S_core_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(pairwise_ancom_HS_16S_core)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_pairwise_ancom_HS_16S_core = pairwise_ancom_HS_16S_core %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_pairwise_ancom_HS_16S_core <- df_pairwise_ancom_HS_16S_core %>%
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
as.factor(df_pairwise_ancom_HS_16S_core$group)

df_pairwise_ancom_HS_16S_core_filtered <- df_pairwise_ancom_HS_16S_core %>%
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




levels(as.factor(df_pairwise_ancom_HS_16S_core_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_HS_core <- HS_16S_core_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_HS_core[is.na(log_table_16S_HS_core)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_HS_core <- exp(log_table_16S_HS_core)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_HS_core <- apply(pseudo_counts_16S_HS_core, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_HS_core <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_HS_core <- t(t(pseudo_counts_16S_HS_core) * (desired_total_count_16S_HS_core / total_counts_16S_HS_core))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_HS_core.df <- as.data.frame(scaled_counts_16S_HS_core)

# Add a column for taxa names
scaled_counts_16S_HS_core.df$taxon <- rownames(scaled_counts_16S_HS_core.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_HS_core.long <- scaled_counts_16S_HS_core.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
HS.micro_data_core.melt <- psmelt(HSpool_16S_core.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_HS.micro_data_core.melt <- HS.micro_data_core.melt %>%
  left_join(scaled_counts_16S_HS_core.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_HS_core_by_poolsize <- BiasAdj_HS.micro_data_core.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_HS_filtered_merged_data_core <- df_pairwise_ancom_HS_16S_core_filtered %>%
  left_join(Family_avg_counts_16S_HS_core_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_HS_core <- as.data.frame(tax_table(HSpool_16S_core.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_HS_core <- df_fig_poolsize_dunnet_16S_HS_filtered_merged_data_core %>%
  left_join(tax_data_16S_HS_core, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_HS_core <- df_with_taxonomy_16S_HS_core %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_HS_core$group = factor(df_with_taxonomy_16S_HS_core$group, levels = c("Pools of 3 vs Individuals", 
                                                                                 "Pools of 6 vs Individuals", 
                                                                                 "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_HS_core <- df_with_taxonomy_16S_HS_core %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_HS_core <- ggplot(df_with_taxonomy_16S_HS_core, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_HS_core %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by timepoint_num") +
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

main_species_plot_16S_HS_core

# Create the taxonomy plot data
taxonomy_plot_data_16S_HS_core <- df_with_taxonomy_16S_HS_core %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_HS_core <- taxonomy_plot_data_16S_HS_core %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_HS_core <- ggplot(taxonomy_plot_data_16S_HS_core) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_HS_core$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_HS_core  <- plot_grid(
  taxonomy_plot_16S_HS_core  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_HS_core  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_HS_core
ggsave("Results/16S/Figures/DA_micro_family_HS_core.jpg", plot = combined_plot_16S_HS_core, width = 20, height = 11, dpi = 300)



##Lowsubset
##data
LSpool_16S.ps
LSpool_16S_individual.ps <- subset_samples(LSpool_16S.ps,Number.in.Pool=="1")
sum(taxa_sums(LSpool_16S_individual.ps)==0)
LSpool_16S_individual.ps <- prune_taxa(taxa_sums(LSpool_16S_individual.ps)>0,LSpool_16S_individual.ps)
sum(sample_sums(LSpool_16S_individual.ps)==0)

LSpool_16S_individual.css <- phyloseq_transform_css(LSpool_16S_individual.ps,log = F)
LSpool_16S_individual.ra <- transform_sample_counts(LSpool_16S_individual.css, function(x) {x/sum(x)}*100)

##now get core-done with ra data
core_members(LSpool_16S_individual.ra, detection = 0, prevalence = 50/100) #gives otus, but there are 46
LSpool_16S_individual.core <- core(LSpool_16S_individual.ps, detection = 0, prevalence = 0.5)

taxa(LSpool_16S_individual.core)

tax_table(LSpool_16S_individual.core)

##what if we do on "raw"/css data?
core_members(LSpool_16S_individual.ps, detection = 0, prevalence = 0.5) #75
core_members(LSpool_16S_individual.css, detection = 0, prevalence = 0.5) #75

###Going on with ps, then will normalize...that ok?
LSpool_16S_3.ps <- subset_samples(LSpool_16S.ps,Number.in.Pool=="3")
sum(taxa_sums(LSpool_16S_3.ps)==0)
LSpool_16S_3.ps <- prune_taxa(taxa_sums(LSpool_16S_3.ps)>0,LSpool_16S_3.ps)
sum(sample_sums(LSpool_16S_3.ps)==0)

core_members(LSpool_16S_3.ps, detection=0, prevalence = 0.5)
LSpool_16S_3.core <- core(LSpool_16S_3.ps, detection = 0, prevalence = 0.5)

LSpool_16S_6.ps <- subset_samples(LSpool_16S.ps, Number.in.Pool=="6")
LSpool_16S_6.ps <- prune_taxa(taxa_sums(LSpool_16S_6.ps)>0,LSpool_16S_6.ps)
core_members(LSpool_16S_6.ps,detection = 0, prevalence = 0.5)
LSpool_16S_6.core <- core(LSpool_16S_6.ps, detection = 0, prevalence = 0.5)

LSpool_16S_12.ps <- subset_samples(LSpool_16S.ps, Number.in.Pool=="12")
LSpool_16S_12.ps <- prune_taxa(taxa_sums(LSpool_16S_12.ps)>0,LSpool_16S_12.ps)
core_members(LSpool_16S_12.ps,detection = 0, prevalence = 0.5)
LSpool_16S_12.core <- core(LSpool_16S_12.ps, detection = 0, prevalence = 0.5)


otu_table_core_1 <- otu_table(LSpool_16S_individual.core,taxa_are_rows = T)
otu_table_core_3 <- otu_table(LSpool_16S_3.core, taxa_are_rows = T)
otu_table_core_6 <- otu_table(LSpool_16S_6.core, taxa_are_rows = T)
otu_table_core_12 <- otu_table(LSpool_16S_12.core, taxa_are_rows = T)

tax_table_core_1 <- tax_table(LSpool_16S_individual.core)
tax_table_core_3 <- tax_table(LSpool_16S_3.core)
tax_table_core_6 <- tax_table(LSpool_16S_6.core)
tax_table_core_12 <- tax_table(LSpool_16S_12.core)


otu_table_core <- merge_phyloseq(otu_table_core_1,otu_table_core_3,otu_table_core_6,otu_table_core_12)
tax_table_core <- merge_phyloseq(tax_table_core_1,tax_table_core_3,tax_table_core_6,tax_table_core_12)

core_microbiome_LS <- merge_phyloseq(otu_table_core,tax_table_core)

LSpool_16S_sampledata <- sample_data(LSpool_16S.ps)
sample_names(LSpool_16S_sampledata)

LSpool_16S_core.ps <- merge_phyloseq(core_microbiome_LS,LSpool_16S_sampledata)
rank_names(LSpool_16S_core.ps)

##let's do unweighted unifrac
#double check no taxa not found in any sample
sum(taxa_sums(LSpool_16S_core.ps)==0)
sum(sample_sums(LSpool_16S_core.ps)==0)

LSpool_16S_core.css <- phyloseq_transform_css(LSpool_16S_core.ps,log = F)
LS_16S_core_uwunifrac.dist <- uwunifrac(LSpool_16S_core.css)

##Needs phylogenetic tree

###ANCOMBC
LS_16S_core_ancom_Family <- ancombc2(data = LSpool_16S_core.ps, assay_name = "counts", tax_level = "Family", 
                                     fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                     group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = TRUE, pairwise = TRUE,
                                     dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
pairwise_ancom_LS_16S_core <- LS_16S_core_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(pairwise_ancom_LS_16S_core)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_pairwise_ancom_LS_16S_core = pairwise_ancom_LS_16S_core %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_pairwise_ancom_LS_16S_core <- df_pairwise_ancom_LS_16S_core %>%
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
as.factor(df_pairwise_ancom_LS_16S_core$group)

df_pairwise_ancom_LS_16S_core_filtered <- df_pairwise_ancom_LS_16S_core %>%
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




levels(as.factor(df_pairwise_ancom_LS_16S_core_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_LS_core <- LS_16S_core_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_LS_core[is.na(log_table_16S_LS_core)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_LS_core <- exp(log_table_16S_LS_core)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_LS_core <- apply(pseudo_counts_16S_LS_core, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_LS_core <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_LS_core <- t(t(pseudo_counts_16S_LS_core) * (desired_total_count_16S_LS_core / total_counts_16S_LS_core))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_LS_core.df <- as.data.frame(scaled_counts_16S_LS_core)

# Add a column for taxa names
scaled_counts_16S_LS_core.df$taxon <- rownames(scaled_counts_16S_LS_core.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_LS_core.long <- scaled_counts_16S_LS_core.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
LS.micro_data_core.melt <- psmelt(LSpool_16S_core.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_LS.micro_data_core.melt <- LS.micro_data_core.melt %>%
  left_join(scaled_counts_16S_LS_core.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_LS_core_by_poolsize <- BiasAdj_LS.micro_data_core.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_LS_filtered_merged_data_core <- df_pairwise_ancom_LS_16S_core_filtered %>%
  left_join(Family_avg_counts_16S_LS_core_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_LS_core <- as.data.frame(tax_table(LSpool_16S_core.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_LS_core <- df_fig_poolsize_dunnet_16S_LS_filtered_merged_data_core %>%
  left_join(tax_data_16S_LS_core, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_LS_core <- df_with_taxonomy_16S_LS_core %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_LS_core$group = factor(df_with_taxonomy_16S_LS_core$group, levels = c("Pools of 3 vs Individuals", 
                                                                                           "Pools of 6 vs Individuals", 
                                                                                           "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_LS_core <- df_with_taxonomy_16S_LS_core %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_LS_core <- ggplot(df_with_taxonomy_16S_LS_core, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_LS_core %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by timepoint_num") +
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

main_species_plot_16S_LS_core

# Create the taxonomy plot data
taxonomy_plot_data_16S_LS_core <- df_with_taxonomy_16S_LS_core %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_LS_core <- taxonomy_plot_data_16S_LS_core %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_LS_core <- ggplot(taxonomy_plot_data_16S_LS_core) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_LS_core$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_LS_core  <- plot_grid(
  taxonomy_plot_16S_LS_core  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_LS_core  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_LS_core

##Random
#High
##data
HRpool_16S.ps
HRpool_16S_individual.ps <- subset_samples(HRpool_16S.ps,Number.in.Pool=="1")
sum(taxa_sums(HRpool_16S_individual.ps)==0)
HRpool_16S_individual.ps <- prune_taxa(taxa_sums(HRpool_16S_individual.ps)>0,HRpool_16S_individual.ps)
sum(sample_sums(HRpool_16S_individual.ps)==0)

HRpool_16S_individual.css <- phyloseq_transform_css(HRpool_16S_individual.ps,log = F)
HRpool_16S_individual.ra <- transform_sample_counts(HRpool_16S_individual.css, function(x) {x/sum(x)}*100)

##now get core-done with ra data
core_members(HRpool_16S_individual.ra, detection = 0, prevalence = 50/100) #gives otus, but there are 46
HRpool_16S_individual.core <- core(HRpool_16S_individual.ps, detection = 0, prevalence = 0.5)

taxa(HRpool_16S_individual.core)

tax_table(HRpool_16S_individual.core)

##what if we do on "raw"/css data?
core_members(HRpool_16S_individual.ps, detection = 0, prevalence = 0.5) #75
core_members(HRpool_16S_individual.css, detection = 0, prevalence = 0.5) #75

###Going on with ps, then will normalize...that ok?
HRpool_16S_3.ps <- subset_samples(HRpool_16S.ps,Number.in.Pool=="3")
sum(taxa_sums(HRpool_16S_3.ps)==0)
HRpool_16S_3.ps <- prune_taxa(taxa_sums(HRpool_16S_3.ps)>0,HRpool_16S_3.ps)
sum(sample_sums(HRpool_16S_3.ps)==0)

core_members(HRpool_16S_3.ps, detection=0, prevalence = 0.5)
HRpool_16S_3.core <- core(HRpool_16S_3.ps, detection = 0, prevalence = 0.5)

HRpool_16S_6.ps <- subset_samples(HRpool_16S.ps, Number.in.Pool=="6")
HRpool_16S_6.ps <- prune_taxa(taxa_sums(HRpool_16S_6.ps)>0,HRpool_16S_6.ps)
core_members(HRpool_16S_6.ps,detection = 0, prevalence = 0.5)
HRpool_16S_6.core <- core(HRpool_16S_6.ps, detection = 0, prevalence = 0.5)

HRpool_16S_12.ps <- subset_samples(HRpool_16S.ps, Number.in.Pool=="12")
HRpool_16S_12.ps <- prune_taxa(taxa_sums(HRpool_16S_12.ps)>0,HRpool_16S_12.ps)
core_members(HRpool_16S_12.ps,detection = 0, prevalence = 0.5)
HRpool_16S_12.core <- core(HRpool_16S_12.ps, detection = 0, prevalence = 0.5)


otu_table_core_1_HR <- otu_table(HRpool_16S_individual.core,taxa_are_rows = T)
otu_table_core_3_HR <- otu_table(HRpool_16S_3.core, taxa_are_rows = T)
otu_table_core_6_HR <- otu_table(HRpool_16S_6.core, taxa_are_rows = T)
otu_table_core_12_HR <- otu_table(HRpool_16S_12.core, taxa_are_rows = T)

tax_table_core_1_HR <- tax_table(HRpool_16S_individual.core)
tax_table_core_3_HR <- tax_table(HRpool_16S_3.core)
tax_table_core_6_HR <- tax_table(HRpool_16S_6.core)
tax_table_core_12_HR <- tax_table(HRpool_16S_12.core)


otu_table_core_HR <- merge_phyloseq(otu_table_core_1_HR,otu_table_core_3_HR,otu_table_core_6_HR,otu_table_core_12_HR)
tax_table_core_HR <- merge_phyloseq(tax_table_core_1_HR,tax_table_core_3_HR,tax_table_core_6_HR,tax_table_core_12_HR)

core_microbiome_HR <- merge_phyloseq(otu_table_core_HR,tax_table_core_HR)

HRpool_16S_sampledata <- sample_data(HRpool_16S.ps)
sample_names(HRpool_16S_sampledata)

HRpool_16S_core.ps <- merge_phyloseq(core_microbiome_HR,HRpool_16S_sampledata)
rank_names(HRpool_16S_core.ps)

##let's do unweighted unifrac
#double check no taxa not found in any sample
sum(taxa_sums(HRpool_16S_core.ps)==0)
sum(sample_sums(HRpool_16S_core.ps)==0)

HRpool_16S_core.css <- phyloseq_transform_css(HRpool_16S_core.ps,log = F)
HR_16S_core_uwunifrac.dist <- uwunifrac(HRpool_16S_core.css)

##Needs phylogenetic tree

###ANCOMBC
HR_16S_core_ancom_Family <- ancombc2(data = HRpool_16S_core.ps, assay_name = "counts", tax_level = "Family", 
                                     fix_formula = "Number.in.Pool", rand_formula = NULL,
                                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                     group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = TRUE, pairwise = TRUE,
                                     dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
pairwise_ancom_HR_16S_core <- HR_16S_core_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(pairwise_ancom_HR_16S_core)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_pairwise_ancom_HR_16S_core = pairwise_ancom_HR_16S_core %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_pairwise_ancom_HR_16S_core <- df_pairwise_ancom_HR_16S_core %>%
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
as.factor(df_pairwise_ancom_HR_16S_core$group)

df_pairwise_ancom_HR_16S_core_filtered <- df_pairwise_ancom_HR_16S_core %>%
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




levels(as.factor(df_pairwise_ancom_HR_16S_core_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_HR_core <- HR_16S_core_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_HR_core[is.na(log_table_16S_HR_core)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_HR_core <- exp(log_table_16S_HR_core)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_HR_core <- apply(pseudo_counts_16S_HR_core, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_HR_core <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_HR_core <- t(t(pseudo_counts_16S_HR_core) * (desired_total_count_16S_HR_core / total_counts_16S_HR_core))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_HR_core.df <- as.data.frame(scaled_counts_16S_HR_core)

# Add a column for taxa names
scaled_counts_16S_HR_core.df$taxon <- rownames(scaled_counts_16S_HR_core.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_HR_core.long <- scaled_counts_16S_HR_core.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
HR.micro_data_core.melt <- psmelt(HRpool_16S_core.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_HR.micro_data_core.melt <- HR.micro_data_core.melt %>%
  left_join(scaled_counts_16S_HR_core.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_HR_core_by_poolsize <- BiasAdj_HR.micro_data_core.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_HR_filtered_merged_data_core <- df_pairwise_ancom_HR_16S_core_filtered %>%
  left_join(Family_avg_counts_16S_HR_core_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_HR_core <- as.data.frame(tax_table(HRpool_16S_core.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_HR_core <- df_fig_poolsize_dunnet_16S_HR_filtered_merged_data_core %>%
  left_join(tax_data_16S_HR_core, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_HR_core <- df_with_taxonomy_16S_HR_core %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_HR_core$group = factor(df_with_taxonomy_16S_HR_core$group, levels = c("Pools of 3 vs Individuals", 
                                                                                           "Pools of 6 vs Individuals", 
                                                                                           "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_HR_core <- df_with_taxonomy_16S_HR_core %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_HR_core <- ggplot(df_with_taxonomy_16S_HR_core, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_HR_core %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by timepoint_num") +
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

main_species_plot_16S_HR_core

# Create the taxonomy plot data
taxonomy_plot_data_16S_HR_core <- df_with_taxonomy_16S_HR_core %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_HR_core <- taxonomy_plot_data_16S_HR_core %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_HR_core <- ggplot(taxonomy_plot_data_16S_HR_core) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_HR_core$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_HR_core  <- plot_grid(
  taxonomy_plot_16S_HR_core  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_HR_core  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_HR_core

##Lowsubset
##data
LRpool_16S.ps
LRpool_16S_individual.ps <- subset_samples(LRpool_16S.ps,Number.in.Pool=="1")
sum(taxa_sums(LRpool_16S_individual.ps)==0)
LRpool_16S_individual.ps <- prune_taxa(taxa_sums(LRpool_16S_individual.ps)>0,LRpool_16S_individual.ps)
sum(sample_sums(LRpool_16S_individual.ps)==0)

LRpool_16S_individual.css <- phyloseq_transform_css(LRpool_16S_individual.ps,log = F)
LRpool_16S_individual.ra <- transform_sample_counts(LRpool_16S_individual.css, function(x) {x/sum(x)}*100)

##now get core-done with ra data
core_members(LRpool_16S_individual.ra, detection = 0, prevalence = 50/100) #gives otus, but there are 46
LRpool_16S_individual.core <- core(LRpool_16S_individual.ps, detection = 0, prevalence = 0.5)

taxa(LRpool_16S_individual.core)

tax_table(LRpool_16S_individual.core)

##what if we do on "raw"/css data?
core_members(LRpool_16S_individual.ps, detection = 0, prevalence = 0.5) #75
core_members(LRpool_16S_individual.css, detection = 0, prevalence = 0.5) #75

###Going on with ps, then will normalize...that ok?
LRpool_16S_3.ps <- subset_samples(LRpool_16S.ps,Number.in.Pool=="3")
sum(taxa_sums(LRpool_16S_3.ps)==0)
LRpool_16S_3.ps <- prune_taxa(taxa_sums(LRpool_16S_3.ps)>0,LRpool_16S_3.ps)
sum(sample_sums(LRpool_16S_3.ps)==0)

core_members(LRpool_16S_3.ps, detection=0, prevalence = 0.5)
LRpool_16S_3.core <- core(LRpool_16S_3.ps, detection = 0, prevalence = 0.5)

LRpool_16S_6.ps <- subset_samples(LRpool_16S.ps, Number.in.Pool=="6")
sum(taxa_sums(LRpool_16S_6.ps)==0)
LRpool_16S_6.ps <- prune_taxa(taxa_sums(LRpool_16S_6.ps)>0,LRpool_16S_6.ps)
core_members(LRpool_16S_6.ps,detection = 0, prevalence = 0.5)
LRpool_16S_6.core <- core(LRpool_16S_6.ps, detection = 0, prevalence = 0.5)

LRpool_16S_12.ps <- subset_samples(LRpool_16S.ps, Number.in.Pool=="12")
sum(taxa_sums(LRpool_16S_12.ps)==0)
LRpool_16S_12.ps <- prune_taxa(taxa_sums(LRpool_16S_12.ps)>0,LRpool_16S_12.ps)
core_members(LRpool_16S_12.ps,detection = 0, prevalence = 0.5)
LRpool_16S_12.core <- core(LRpool_16S_12.ps, detection = 0, prevalence = 0.5)


otu_table_core_1_LR <- otu_table(LRpool_16S_individual.core,taxa_are_rows = T)
otu_table_core_3_LR <- otu_table(LRpool_16S_3.core, taxa_are_rows = T)
otu_table_core_6_LR <- otu_table(LRpool_16S_6.core, taxa_are_rows = T)
otu_table_core_12_LR <- otu_table(LRpool_16S_12.core, taxa_are_rows = T)

tax_table_core_1_LR <- tax_table(LRpool_16S_individual.core)
tax_table_core_3_LR <- tax_table(LRpool_16S_3.core)
tax_table_core_6_LR <- tax_table(LRpool_16S_6.core)
tax_table_core_12_LR <- tax_table(LRpool_16S_12.core)


otu_table_core_LR <- merge_phyloseq(otu_table_core_1_LR,otu_table_core_3_LR,otu_table_core_6_LR,otu_table_core_12_LR)
tax_table_core_LR <- merge_phyloseq(tax_table_core_1_LR,tax_table_core_3_LR,tax_table_core_6_LR,tax_table_core_12_LR)

core_microbiome_LR <- merge_phyloseq(otu_table_core_LR,tax_table_core_LR)

LRpool_16S_sampledata <- sample_data(LRpool_16S.ps)
sample_names(LRpool_16S_sampledata)
sample_names(core_microbiome_LR)
sample_names(LRpool_16S_sampledata)


LRpool_16S_core.ps <- merge_phyloseq(core_microbiome_LR,LRpool_16S_sampledata)
rank_names(LRpool_16S_core.ps)
sample_names(LRpool_16S_core.ps)

##let's do unweighted unifrac
#double check no taxa not found in any sample
sum(taxa_sums(LRpool_16S_core.ps)==0)
sum(sample_sums(LRpool_16S_core.ps)==0)

LRpool_16S_core.css <- phyloseq_transform_css(LRpool_16S_core.ps,log = F)
LR_16S_core_uwunifrac.dist <- uwunifrac(LRpool_16S_core.css)

##Needs phylogenetic tree

str(sample_data(LRpool_16S_core.ps))

###ANCOMBC
LR_16S_core_ancom_Family <- ancombc2(data = LRpool_16S_core.ps, assay_name = "counts", tax_level = "Family", 
                                     fix_formula = "Number.in.Pool", rand_formula = NULL,
                                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                     group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = TRUE, pairwise = TRUE,
                                     dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
pairwise_ancom_LR_16S_core <- LR_16S_core_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(pairwise_ancom_LR_16S_core)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_pairwise_ancom_LR_16S_core = pairwise_ancom_LR_16S_core %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_pairwise_ancom_LR_16S_core <- df_pairwise_ancom_LR_16S_core %>%
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
as.factor(df_pairwise_ancom_LR_16S_core$group)

df_pairwise_ancom_LR_16S_core_filtered <- df_pairwise_ancom_LR_16S_core %>%
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




levels(as.factor(df_pairwise_ancom_LR_16S_core_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_LR_core <- LR_16S_core_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_LR_core[is.na(log_table_16S_LR_core)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_LR_core <- exp(log_table_16S_LR_core)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_LR_core <- apply(pseudo_counts_16S_LR_core, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_LR_core <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_LR_core <- t(t(pseudo_counts_16S_LR_core) * (desired_total_count_16S_LR_core / total_counts_16S_LR_core))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_LR_core.df <- as.data.frame(scaled_counts_16S_LR_core)

# Add a column for taxa names
scaled_counts_16S_LR_core.df$taxon <- rownames(scaled_counts_16S_LR_core.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_LR_core.long <- scaled_counts_16S_LR_core.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
LR.micro_data_core.melt <- psmelt(LRpool_16S_core.ps)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_LR.micro_data_core.melt <- LR.micro_data_core.melt %>%
  left_join(scaled_counts_16S_LR_core.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_LR_core_by_poolsize <- BiasAdj_LR.micro_data_core.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_LR_filtered_merged_data_core <- df_pairwise_ancom_LR_16S_core_filtered %>%
  left_join(Family_avg_counts_16S_LR_core_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_LR_core <- as.data.frame(tax_table(LRpool_16S_core.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_LR_core <- df_fig_poolsize_dunnet_16S_LR_filtered_merged_data_core %>%
  left_join(tax_data_16S_LR_core, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_LR_core <- df_with_taxonomy_16S_LR_core %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_LR_core$group = factor(df_with_taxonomy_16S_LR_core$group, levels = c("Pools of 3 vs Individuals", 
                                                                                           "Pools of 6 vs Individuals", 
                                                                                           "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_LR_core <- df_with_taxonomy_16S_LR_core %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_LR_core <- ggplot(df_with_taxonomy_16S_LR_core, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_LR_core %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by timepoint_num") +
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

main_species_plot_16S_LR_core

# Create the taxonomy plot data
taxonomy_plot_data_16S_LR_core <- df_with_taxonomy_16S_LR_core %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_LR_core <- taxonomy_plot_data_16S_LR_core %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_LR_core <- ggplot(taxonomy_plot_data_16S_LR_core) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_LR_core$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_LR_core  <- plot_grid(
  taxonomy_plot_16S_LR_core  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_LR_core  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_LR_core

##Relative abundance plots
#trying core microbiome on overall phyloseq object instead of splitting and recombining
HSpool_16S_all.core <- core(HSpool_16S.ps, detection = 0, prevalence = 0.5)

HS_16S_core_all_ancom_Family <- ancombc2(data = HSpool_16S_all.core, assay_name = "counts", tax_level = "Family", 
                                     fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                                     p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                     group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = TRUE, pairwise = TRUE,
                                     dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
dunn_ancom_HS_16S_core_all <- HS_16S_core_all_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(dunn_ancom_HS_16S_core_all)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_dunn_ancom_HS_16S_core_all =dunn_ancom_HS_16S_core_all %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_dunn_ancom_HS_16S_core_all <- df_dunn_ancom_HS_16S_core_all %>%
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
as.factor(df_dunn_ancom_HS_16S_core_all$group)

df_dunn_ancom_HS_16S_core_all_filtered <- df_dunn_ancom_HS_16S_core_all %>%
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




levels(as.factor(df_dunn_ancom_HS_16S_core_all_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_HS_core_all <- HS_16S_core_all_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_HS_core_all[is.na(log_table_16S_HS_core_all)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_HS_core_all <- exp(log_table_16S_HS_core_all)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_HS_core_all <- apply(pseudo_counts_16S_HS_core_all, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_HS_core_all <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_HS_core_all <- t(t(pseudo_counts_16S_HS_core_all) * (desired_total_count_16S_HS_core_all / total_counts_16S_HS_core_all))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_HS_core_all.df <- as.data.frame(scaled_counts_16S_HS_core_all)

# Add a column for taxa names
scaled_counts_16S_HS_core_all.df$taxon <- rownames(scaled_counts_16S_HS_core_all.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_HS_core_all.long <- scaled_counts_16S_HS_core_all.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
HS.micro_data_core_all.melt <- psmelt(HSpool_16S_all.core)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_HS.micro_data_core_all.melt <- HS.micro_data_core_all.melt %>%
  left_join(scaled_counts_16S_HS_core_all.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_HS_core_all_by_poolsize <- BiasAdj_HS.micro_data_core_all.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_HS_filtered_merged_data_core_all <- df_dunn_ancom_HS_16S_core_all_filtered %>%
  left_join(Family_avg_counts_16S_HS_core_all_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_HS_core_all <- as.data.frame(tax_table(HSpool_16S_all.core))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_HS_core_all <- df_fig_poolsize_dunnet_16S_HS_filtered_merged_data_core_all %>%
  left_join(tax_data_16S_HS_core_all, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_HS_core_all <- df_with_taxonomy_16S_HS_core_all %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_HS_core_all$group = factor(df_with_taxonomy_16S_HS_core_all$group, levels = c("Pools of 3 vs Individuals", 
                                                                                           "Pools of 6 vs Individuals", 
                                                                                           "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_HS_core_all <- df_with_taxonomy_16S_HS_core_all %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_HS_core_all <- ggplot(df_with_taxonomy_16S_HS_core_all, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_HS_core_all %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size, HIGH-SUB") +
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

main_species_plot_16S_HS_core_all

# Create the taxonomy plot data
taxonomy_plot_data_16S_HS_core_all <- df_with_taxonomy_16S_HS_core_all %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_HS_core_all <- taxonomy_plot_data_16S_HS_core_all %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_HS_core_all <- ggplot(taxonomy_plot_data_16S_HS_core_all) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_HS_core_all$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_HS_core_all  <- plot_grid(
  taxonomy_plot_16S_HS_core_all  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_HS_core_all  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widths = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_HS_core_all
ggsave("Results/16S/Figures/DA_micro_family_HS_core_all.jpg", plot = combined_plot_16S_HS_core_all, width = 20, height = 11, dpi = 300)

#trying core microbiome on overall phyloseq object instead of splitting and recombining
HRpool_16S_all.core <- core(HRpool_16S.ps, detection = 0, prevalence = 0.5)

HR_16S_core_all_ancom_Family <- ancombc2(data = HRpool_16S_all.core, assay_name = "counts", tax_level = "Family", 
                                         fix_formula = "Number.in.Pool", rand_formula = NULL,
                                         p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                         group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                         alpha = 0.05, n_cl = 2, verbose = TRUE,
                                         global = TRUE, pairwise = TRUE,
                                         dunnet = TRUE, trend = FALSE)

   # Extract the pairwise test results
dunn_ancom_HR_16S_core_all <- HR_16S_core_all_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(dunn_ancom_HR_16S_core_all)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_dunn_ancom_HR_16S_core_all =dunn_ancom_HR_16S_core_all %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_dunn_ancom_HR_16S_core_all <- df_dunn_ancom_HR_16S_core_all %>%
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
as.factor(df_dunn_ancom_HR_16S_core_all$group)

df_dunn_ancom_HR_16S_core_all_filtered <- df_dunn_ancom_HR_16S_core_all %>%
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




levels(as.factor(df_dunn_ancom_HR_16S_core_all_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_HR_core_all <- HR_16S_core_all_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_HR_core_all[is.na(log_table_16S_HR_core_all)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_HR_core_all <- exp(log_table_16S_HR_core_all)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_HR_core_all <- apply(pseudo_counts_16S_HR_core_all, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_HR_core_all <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_HR_core_all <- t(t(pseudo_counts_16S_HR_core_all) * (desired_total_count_16S_HR_core_all / total_counts_16S_HR_core_all))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_HR_core_all.df <- as.data.frame(scaled_counts_16S_HR_core_all)

# Add a column for taxa names
scaled_counts_16S_HR_core_all.df$taxon <- rownames(scaled_counts_16S_HR_core_all.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_HR_core_all.long <- scaled_counts_16S_HR_core_all.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
HR.micro_data_core_all.melt <- psmelt(HRpool_16S_all.core)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_HR.micro_data_core_all.melt <- HR.micro_data_core_all.melt %>%
  left_join(scaled_counts_16S_HR_core_all.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_HR_core_all_by_poolsize <- BiasAdj_HR.micro_data_core_all.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_HR_filtered_merged_data_core_all <- df_dunn_ancom_HR_16S_core_all_filtered %>%
  left_join(Family_avg_counts_16S_HR_core_all_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_HR_core_all <- as.data.frame(tax_table(HRpool_16S_all.core))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_HR_core_all <- df_fig_poolsize_dunnet_16S_HR_filtered_merged_data_core_all %>%
  left_join(tax_data_16S_HR_core_all, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_HR_core_all <- df_with_taxonomy_16S_HR_core_all %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_HR_core_all$group = factor(df_with_taxonomy_16S_HR_core_all$group, levels = c("Pools of 3 vs Individuals", 
                                                                                                   "Pools of 6 vs Individuals", 
                                                                                                   "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_HR_core_all <- df_with_taxonomy_16S_HR_core_all %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_HR_core_all <- ggplot(df_with_taxonomy_16S_HR_core_all, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_HR_core_all %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size, HIGH-RAND") +
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

main_species_plot_16S_HR_core_all

# Create the taxonomy plot data
taxonomy_plot_data_16S_HR_core_all <- df_with_taxonomy_16S_HR_core_all %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_HR_core_all <- taxonomy_plot_data_16S_HR_core_all %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_HR_core_all <- ggplot(taxonomy_plot_data_16S_HR_core_all) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_HR_core_all$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_HR_core_all  <- plot_grid(
  taxonomy_plot_16S_HR_core_all  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_HR_core_all  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widtHR = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_HR_core_all
ggsave("Results/16S/Figures/DA_micro_family_HR_core_all.jpg", plot = combined_plot_16S_HR_core_all, width = 20, height = 11, dpi = 300)


#Getting ancom done
LSpool_16S_all.core <- core(LSpool_16S.ps, detection = 0, prevalence = 0.5)
LRpool_16S_all.core <- core(LRpool_16S.ps,detection = 0, prevalence = 0.5)

LS_16S_core_all_ancom_Family <- ancombc2(data = LSpool_16S_all.core, assay_name = "counts", tax_level = "Family", 
                                         fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                                         p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                         group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                         alpha = 0.05, n_cl = 2, verbose = TRUE,
                                         global = TRUE, pairwise = TRUE,
                                         dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
dunn_ancom_LS_16S_core_all <- LS_16S_core_all_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(dunn_ancom_LS_16S_core_all)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_dunn_ancom_LS_16S_core_all =dunn_ancom_LS_16S_core_all %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_dunn_ancom_LS_16S_core_all <- df_dunn_ancom_LS_16S_core_all %>%
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
as.factor(df_dunn_ancom_LS_16S_core_all$group)

df_dunn_ancom_LS_16S_core_all_filtered <- df_dunn_ancom_LS_16S_core_all %>%
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




levels(as.factor(df_dunn_ancom_LS_16S_core_all_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_LS_core_all <- LS_16S_core_all_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_LS_core_all[is.na(log_table_16S_LS_core_all)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_LS_core_all <- exp(log_table_16S_LS_core_all)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_LS_core_all <- apply(pseudo_counts_16S_LS_core_all, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_LS_core_all <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_LS_core_all <- t(t(pseudo_counts_16S_LS_core_all) * (desired_total_count_16S_LS_core_all / total_counts_16S_LS_core_all))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_LS_core_all.df <- as.data.frame(scaled_counts_16S_LS_core_all)

# Add a column for taxa names
scaled_counts_16S_LS_core_all.df$taxon <- rownames(scaled_counts_16S_LS_core_all.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_LS_core_all.long <- scaled_counts_16S_LS_core_all.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
LS.micro_data_core_all.melt <- psmelt(LSpool_16S_all.core)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_LS.micro_data_core_all.melt <- LS.micro_data_core_all.melt %>%
  left_join(scaled_counts_16S_LS_core_all.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_LS_core_all_by_poolsize <- BiasAdj_LS.micro_data_core_all.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_LS_filtered_merged_data_core_all <- df_dunn_ancom_LS_16S_core_all_filtered %>%
  left_join(Family_avg_counts_16S_LS_core_all_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_LS_core_all <- as.data.frame(tax_table(LSpool_16S_all.core))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_LS_core_all <- df_fig_poolsize_dunnet_16S_LS_filtered_merged_data_core_all %>%
  left_join(tax_data_16S_LS_core_all, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_LS_core_all <- df_with_taxonomy_16S_LS_core_all %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_LS_core_all$group = factor(df_with_taxonomy_16S_LS_core_all$group, levels = c("Pools of 3 vs Individuals", 
                                                                                                   "Pools of 6 vs Individuals", 
                                                                                                   "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_LS_core_all <- df_with_taxonomy_16S_LS_core_all %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_LS_core_all <- ggplot(df_with_taxonomy_16S_LS_core_all, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_LS_core_all %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size, LOW-SUB") +
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

main_species_plot_16S_LS_core_all

# Create the taxonomy plot data
taxonomy_plot_data_16S_LS_core_all <- df_with_taxonomy_16S_LS_core_all %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_LS_core_all <- taxonomy_plot_data_16S_LS_core_all %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_LS_core_all <- ggplot(taxonomy_plot_data_16S_LS_core_all) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_LS_core_all$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_LS_core_all  <- plot_grid(
  taxonomy_plot_16S_LS_core_all  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_LS_core_all  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widtLS = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_LS_core_all
ggsave("Results/16S/Figures/DA_micro_family_LS_core_all.jpg", plot = combined_plot_16S_LS_core_all, width = 20, height = 11, dpi = 300)



LR_16S_core_all_ancom_Family <- ancombc2(data = LRpool_16S_all.core, assay_name = "counts", tax_level = "Family", 
                                         fix_formula = "Number.in.Pool", rand_formula = NULL,
                                         p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                         group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                                         alpha = 0.05, n_cl = 2, verbose = TRUE,
                                         global = TRUE, pairwise = TRUE,
                                         dunnet = TRUE, trend = FALSE)

# Extract the pairwise test results
dunn_ancom_LR_16S_core_all <- LR_16S_core_all_ancom_Family$res_dunn

# Check the structure of the pairwise results
str(dunn_ancom_LR_16S_core_all)

# Filter the results to focus on comparisons between timepoint_nums
# For example, if you want to look at a specific species or taxon, you can filter accordingly
df_dunn_ancom_LR_16S_core_all =dunn_ancom_LR_16S_core_all %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 


df_dunn_ancom_LR_16S_core_all <- df_dunn_ancom_LR_16S_core_all %>%
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
as.factor(df_dunn_ancom_LR_16S_core_all$group)

df_dunn_ancom_LR_16S_core_all_filtered <- df_dunn_ancom_LR_16S_core_all %>%
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




levels(as.factor(df_dunn_ancom_LR_16S_core_all_filtered$group))

#### FECAL_NEG Dotplot - but with bias corrected abundances #######

# Assuming 'ancom_output_FECAL_NEG$bias_correct_log_table' is your bias-corrected log-transformed data
log_table_16S_LR_core_all <- LR_16S_core_all_ancom_Family$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_16S_LR_core_all[is.na(log_table_16S_LR_core_all)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_16S_LR_core_all <- exp(log_table_16S_LR_core_all)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_16S_LR_core_all <- apply(pseudo_counts_16S_LR_core_all, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_16S_LR_core_all <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_16S_LR_core_all <- t(t(pseudo_counts_16S_LR_core_all) * (desired_total_count_16S_LR_core_all / total_counts_16S_LR_core_all))

# Convert scaled_counts_FECAL_NEG to a data frame
scaled_counts_16S_LR_core_all.df <- as.data.frame(scaled_counts_16S_LR_core_all)

# Add a column for taxa names
scaled_counts_16S_LR_core_all.df$taxon <- rownames(scaled_counts_16S_LR_core_all.df)

# Convert from wide to long format using pivot_longer
scaled_counts_16S_LR_core_all.long <- scaled_counts_16S_LR_core_all.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
LR.micro_data_core_all.melt <- psmelt(LRpool_16S_all.core)

# Merge data
# Merge the melted data with the scaled counts on 'taxon' and 'Sample'
BiasAdj_LR.micro_data_core_all.melt <- LR.micro_data_core_all.melt %>%
  left_join(scaled_counts_16S_LR_core_all.long, by = c("Family" = "taxon", "Sample" = "Sample"))

Family_avg_counts_16S_LR_core_all_by_poolsize <- BiasAdj_LR.micro_data_core_all.melt %>%
  dplyr::group_by(Family, Number.in.Pool) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE)) %>%
  ungroup() %>%
  dplyr::mutate(Number.in.Pool_label = dplyr::case_when(
    Number.in.Pool == 3 ~ "Pools of 3 vs Individuals",
    Number.in.Pool == 6 ~ "Pools of 6 vs Individuals",
    Number.in.Pool == 12 ~ "Pools of 12 vs Individuals",
    TRUE ~ as.character(Number.in.Pool)  # Keep other timepoint_num values unchanged
  ))



df_fig_poolsize_dunnet_16S_LR_filtered_merged_data_core_all <- df_dunn_ancom_LR_16S_core_all_filtered %>%
  left_join(Family_avg_counts_16S_LR_core_all_by_poolsize, by = c("group" = "Number.in.Pool_label", "taxon" = "Family"))

# Extract taxonomic data from the phyloseq object
tax_data_16S_LR_core_all <- as.data.frame(tax_table(LRpool_16S_all.core))

# Join the taxonomic data with your data frame
df_with_taxonomy_16S_LR_core_all <- df_fig_poolsize_dunnet_16S_LR_filtered_merged_data_core_all %>%
  left_join(tax_data_16S_LR_core_all, by = c("taxon" = "Family"),relationship = "many-to-many")  # Assuming 'Family' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_16S_LR_core_all <- df_with_taxonomy_16S_LR_core_all %>%
  arrange(Phylum) %>%  # Sort the data by Class
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order

df_with_taxonomy_16S_LR_core_all$group = factor(df_with_taxonomy_16S_LR_core_all$group, levels = c("Pools of 3 vs Individuals", 
                                                                                                   "Pools of 6 vs Individuals", 
                                                                                                   "Pools of 12 vs Individuals"))

# Ensure unique taxon for labeling purposes
# Can try to fix it later, has to do with how I'm merging the taxonomy to the counts which creates duplicates for the lower taxa levels.
# This is not necassry for the AMR data.
df_with_taxonomy_16S_LR_core_all <- df_with_taxonomy_16S_LR_core_all %>%
  group_by(group, taxon) %>%
  filter(row_number() == 1) %>%
  ungroup()

## Trying to add taxonomy
main_species_plot_16S_LR_core_all <- ggplot(df_with_taxonomy_16S_LR_core_all, aes(x = value, y = Class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Vertical line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  # Add labels for points where color is "darkred" with balanced horizontal and vertical repelling
  geom_text_repel(data = df_with_taxonomy_16S_LR_core_all %>% dplyr::filter(color == "darkred"),
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
  labs(x = "Log-Fold Change", y = "Class", title = " logFC of Family by Pool Size, LOW-RAND") +
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

main_species_plot_16S_LR_core_all

# Create the taxonomy plot data
taxonomy_plot_data_16S_LR_core_all <- df_with_taxonomy_16S_LR_core_all %>%
  distinct(Phylum, Class, taxon)

# Modify the data to create new columns with the "label_" prefix
taxonomy_plot_data_16S_LR_core_all <- taxonomy_plot_data_16S_LR_core_all %>%
  dplyr::group_by(Phylum) %>%
  dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
  ungroup() 

# Create the updated taxonomy plot
taxonomy_plot_16S_LR_core_all <- ggplot(taxonomy_plot_data_16S_LR_core_all) +
  geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
  #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
  scale_y_discrete(limits = levels(taxonomy_plot_data_16S_LR_core_all$Class)) +  # Ensure y-axis matches the taxon order
  theme_void() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
  )


combined_plot_16S_LR_core_all  <- plot_grid(
  taxonomy_plot_16S_LR_core_all  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
  main_species_plot_16S_LR_core_all  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
  ncol = 2, 
  rel_widtLR = c(0.5,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb"
)

combined_plot_16S_LR_core_all
ggsave("Results/16S/Figures/DA_micro_family_LR_core_all.jpg", plot = combined_plot_16S_LR_core_all, width = 20, height = 11, dpi = 300)



