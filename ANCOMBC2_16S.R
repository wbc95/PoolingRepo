###DataSetup
data_16S_50k #158 samples, 5 dropped

sum(taxa_sums(data_16S_50k)==0) #42 taxa
any(sample_sums(data_16S_50k)==0)


HSpool_16S.ps <- subset_samples(data_16S_50k, PrevH=="Y"&PoolingType!="Random")
View(sample_data(HSpool_16S.ps))

sample_data(HSpool_16S.ps)$StatsPool<-c("A","B","C","D","E","F","A","B","C","D","E","F","A","B","C","D","F","B","E","A","E","B","E","E","F","C","F","A","F",
                                        "D","D","E","C","F","F","F","D","C","A","B","A","F","C","B","F","B","F","F","D","D","B","D","E","B","C","E","A","D",
                                        "D","C","E","A","A","B","C","F","D","E","C","A","D","B","C","D","F","A","E","E","A","B","A","C","C","C","E","A","D","B")
sample_data(HSpool_16S.ps)$StatsPool <- factor(sample_data(HSpool_16S.ps)$StatsPool, levels = c("A","B","C","D","E","F"))

sample_data(HSpool_16S.ps)$Number.in.Pool <- factor(sample_data(HSpool_16S.ps)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(HSpool_16S.ps)

#Running model
ancom_Family_16S_HS = ancombc2(data = HSpool_16S.ps, assay_name = "counts", tax_level = "Family",
                          fix_formula = "StatsPool + Number.in.Pool", rand_formula = NULL,
                          p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                          group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          global = TRUE, pairwise = TRUE,
                          dunnet = TRUE, trend = FALSE,
                          iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                          trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                      nrow = 2,
                                                                      byrow = TRUE),matrix(c(-1, 0, 1, -1),
                                                                                           nrow = 2,
                                                                                           byrow = TRUE)),
                                               node = list(2, 2),
                                               solver = "ECOS",
                                               B = 100))


ancom_Phylum_16S_HS = ancombc2(data = HSpool_16S.ps, assay_name = "counts", tax_level = "Phylum",
                               fix_formula = "StatsPool + Number.in.Pool", rand_formula = NULL,
                               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               dunnet = TRUE, trend = FALSE,
                               iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                           nrow = 2,
                                                                           byrow = TRUE),matrix(c(-1, 0, 1, -1),
                                                                                                nrow = 2,
                                                                                                byrow = TRUE)),
                                                    node = list(2, 2),
                                                    solver = "ECOS",
                                                    B = 100))

res_pair_16S_HS = ancom_Family_16S_HS$res_pair
#glob_HS_16S = ancom_Family_16S_HS$res_global
dunnet_16S_HS = ancom_Family_16S_HS$res_dunn
BCabundancetable_16S_HS = ancom_Family_16S_HS$bias_correct_log_table
dunnet_16S_HS_phyl =ancom_Phylum_16S_HS$res_dunn


write.csv(dunnet_16S_HS,"dunnet_16S_HS.csv")
write.csv(BCabundancetable_16S_HS, "BCLogAbundance_16S_HS.csv")


df_fig_dunnet_16S_HS1 = dunnet_16S_HS %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc3,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

diff_16S_HS_3 = dunnet_16S_HS %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_16S_HS2 = dunnet_16S_HS %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12==1) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_16S_HS = df_fig_dunnet_16S_HS1 %>%
  dplyr::left_join(df_fig_dunnet_16S_HS2, by=c("taxon","group"))

df_fig_dunnet_16S_HS$group = recode(df_fig_dunnet_16S_HS$group,
                              `lfc1`= "Pools of 3",
                              `lfc2`= "Pools of 6",
                              `lfc3`= "Pools of 12")

df_fig_dunnet_16S_HS$group =factor(df_fig_dunnet_16S_HS$group,
                                   levels = c("Pools of 3",
                                              "Pools of 6",
                                              "Pools of 12"))

lo_16S_HS_dunn = floor(min(df_fig_dunnet_16S_HS$value))
up_16S_HS_dunn = ceiling(max(df_fig_dunnet_16S_HS$value))
mid_16S_HS_dunn = (lo_16S_HS_dunn+up_16S_HS_dunn)/2

fig_dunn_16S_HS = df_fig_dunnet_16S_HS %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_HS_dunn,up_16S_HS_dunn),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

df_fig_res_pair_16S_HS1 = res_pair_16S_HS %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1|
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0),
                lfc4=ifelse(diff_Number.in.Pool6_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool6_Number.in.Pool3,2),0),
                lfc5=ifelse(diff_Number.in.Pool12_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool3, 2),0),
                lfc6=ifelse(diff_Number.in.Pool12_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool6, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc6,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_HS2 = res_pair_16S_HS %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12 ==1 |
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black"),
                lfc4= ifelse(passed_ss_Number.in.Pool6_Number.in.Pool3 ==1 & diff_Number.in.Pool6_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc5= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool3 ==1 & diff_Number.in.Pool12_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc6= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool6 ==1 & diff_Number.in.Pool12_Number.in.Pool6 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_HS = df_fig_res_pair_16S_HS1 %>%
  dplyr::left_join(df_fig_res_pair_16S_HS2, by=c("taxon","group"))

df_fig_res_pair_16S_HS$group = recode(df_fig_res_pair_16S_HS$group,
                                    `lfc1`= "3 vs 1",
                                    `lfc2`= "6 vs 1",
                                    `lfc3`= "12 vs 1",
                                    `lfc4`= "6 vs 3",
                                    `lfc5`= "12 vs 3",
                                    `lfc6`= "12 vs 6")

df_fig_res_pair_16S_HS$group =factor(df_fig_res_pair_16S_HS$group,
                                   levels = c("3 vs 1",
                                              "6 vs 1",
                                              "12 vs 1",
                                              "6 vs 3",
                                              "12 vs 3",
                                              "12 vs 6"))

lo_16S_HS_pair = floor(min(df_fig_res_pair_16S_HS$value))
up_16S_HS_pair = ceiling(max(df_fig_res_pair_16S_HS$value))
mid_16S_HS_pair = (lo_16S_HS_pair+up_16S_HS_pair)/2

fig_pair_16S_HS = df_fig_res_pair_16S_HS %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_HS_pair,up_16S_HS_pair),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

##TRying at the phylum level
df_fig_dunnet_16S_HS_phyl1 = dunnet_16S_HS_phyl %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc3,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

diff_16S_HS_3 = dunnet_16S_HS_phyl %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_16S_HS_phyl2 = dunnet_16S_HS_phyl %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_16S_HS_phyl = df_fig_dunnet_16S_HS_phyl1 %>%
  dplyr::left_join(df_fig_dunnet_16S_HS_phyl2, by=c("taxon","group"))

df_fig_dunnet_16S_HS_phyl$group = recode(df_fig_dunnet_16S_HS_phyl$group,
                                    `lfc1`= "Pools of 3",
                                    `lfc2`= "Pools of 6",
                                    `lfc3`= "Pools of 12")

df_fig_dunnet_16S_HS_phyl$group =factor(df_fig_dunnet_16S_HS_phyl$group,
                                   levels = c("Pools of 3",
                                              "Pools of 6",
                                              "Pools of 12"))

lo_16S_HS_dunn_phyl = floor(min(df_fig_dunnet_16S_HS_phyl$value))
up_16S_HS_dunn_phyl = ceiling(max(df_fig_dunnet_16S_HS_phyl$value))
mid_16S_HS_dunn_phyl = (lo_16S_HS_dunn+up_16S_HS_dunn)/2

fig_dunn_16S_HS = df_fig_dunnet_16S_HS_phyl %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_HS_dunn_phyl,up_16S_HS_dunn_phyl),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_HS_16S = res_pair_HS_16S %>%
  dplyr::select(taxon, contains("Number.in.Pool"))


df_HS_fig_glob = glob_HS_16S%>%
  filter(passed_ss ==1 & diff_abn==1)

# Extract the pairwise comparison results
df_fig_HS_16S = df_HS_16S %>%
  filter((diff_Number.in.Pool3 == 1 & passed_ss_Number.in.Pool3==1)| (diff_Number.in.Pool6 == 1 & passed_ss_Number.in.Pool6==1)|
           (diff_Number.in.Pool12 == 1 & passed_ss_Number.in.Pool12==1)| (diff_Number.in.Pool6_Number.in.Pool3  == 1 & passed_ss_Number.in.Pool6_Number.in.Pool3==1)|
           (diff_Number.in.Pool12_Number.in.Pool3  == 1 & passed_ss_Number.in.Pool12_Number.in.Pool3==1)| (diff_Number.in.Pool12_Number.in.Pool6 == 1 & passed_ss_Number.in.Pool12_Number.in.Pool6==1 )) %>%
  mutate(lfc_Number.in.Pool3 = ifelse(diff_Number.in.Pool3 == 1,  # 3v1
                                    lfc_Number.in.Pool3, 0),
         lfc_Number.in.Pool6 = ifelse(diff_Number.in.Pool6 == 1,  # 6v1
                                      lfc_Number.in.Pool6, 0),
         lfc_Number.in.Pool12 = ifelse(diff_Number.in.Pool12 == 1,  # 12v1
                                      lfc_Number.in.Pool12, 0),
         lfc_Number.in.Pool6_Number.in.Pool3 = ifelse(diff_Number.in.Pool6_Number.in.Pool3 == 1, #6v3
                                                      lfc_Number.in.Pool6_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool3 = ifelse(diff_Number.in.Pool12_Number.in.Pool3 == 1,  # 12v3
                                      lfc_Number.in.Pool12_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool6 = ifelse(diff_Number.in.Pool12_Number.in.Pool6 == 1,  # 12v6
                                       lfc_Number.in.Pool12_Number.in.Pool6, 0),
  ) %>%
  transmute(taxon, 
            `3v1` = round(lfc_Number.in.Pool3, 2),
            `6v1` = round(lfc_Number.in.Pool6, 2),
            `12v1` = round(lfc_Number.in.Pool12, 2),
            `6v3` = round(lfc_Number.in.Pool6_Number.in.Pool3, 2),
            `12v3` = round(lfc_Number.in.Pool12_Number.in.Pool3, 2),
            `12v6` = round(lfc_Number.in.Pool12_Number.in.Pool6,2)) %>%
  pivot_longer(cols = `3v1`:`6v1`:`12v1`:`6v3`:`12v3`:`12v6`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Change order of X-value
df_fig_HS_16S
df_fig_HS_16S$group <- factor(df_fig_HS_16S$group, levels = c("3v1","6v1","12v1","6v3","12v3","12v6"))
unique(df_fig_HS_16S$taxon)

unique(df_HS_fig_glob$taxon)

df_fig_HS_16S_filtered <- subset(df_fig_HS_16S, taxon=="unclassified Unassigned" | taxon== "Acholeplasmataceae"| taxon=="RF39"| taxon=="Mycoplasmataceae"|
                                   taxon=="Lactobacillaceae"| taxon=="Streptococcaceae"| taxon=="Selenomonadaceae"| taxon=="Campylobacteraceae"|
                                   taxon== "Microbacteriaceae" | taxon=="Neisseriaceae"| taxon=="Comamonadaceae"| taxon=="Moraxellaceae"|
                                   taxon=="Succinivibrionaceae"| taxon=="Methanobacteriaceae"| taxon=="Spirochaetaceae"| taxon=="Chitinophagaceae")

# Plot heatmap
lo = floor(min(df_fig_HS_16S$value))
up = ceiling(max(df_fig_HS_16S$value))
mid = (lo + up)/2



fig_HS_16S = df_fig_HS_16S %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by combined order") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_HS_16S

tab_pa_HS_16S = ancom_Family_1_HS$zero_ind



#filtering
lo_HS_16S_filtered = floor(min(df_fig_HS_16S_filtered$value))
up_HS_16S_filtered = ceiling(max(df_fig_HS_16S_filtered$value))
mid_HS_16S_filtered = (lo_HS_16S_filtered + up_HS_16S_filtered)/2

  
fig_HS_16S_filtered = df_fig_HS_16S_filtered %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo_HS_16S_filtered, up_HS_16S_filtered),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by combined order") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_HS_16S_filtered

##LOW SUBSET
LSpool_16S.ps <- subset_samples(data_16S_50k, PrevL=="Y"&PoolingType!="Random")
View(sample_data(LSpool_16S.ps))

sample_data(LSpool_16S.ps)$LSPool

sample_data(LSpool_16S.ps)$StatsPool<-c("A","B","C","D","E","F","A","B","C","D","E","F","A","B","C","D","E","F","D","F","E","E","F","A","B","D","F","D","E",
                                        "C","C","F","D","F","C","B","A","E","A","A","E","F","F","F","A","D","D","A","E","F","B","A","B","E","A","E","C","A",
                                        "F","B","C","A","C","C","B","C","B","D","C","D","E","E","A","B","D","B","B","F","F","A","B","D","E","C","B","D","C","C")
sample_data(LSpool_16S.ps)$StatsPool <- factor(sample_data(LSpool_16S.ps)$StatsPool, levels = c("A","B","C","D","E","F"))

sample_data(LSpool_16S.ps)$Number.in.Pool <- factor(sample_data(LSpool_16S.ps)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(LSpool_16S.ps)

#Running model

#Already run, but changed name
#ancom_Family_16S_LS <- ancom_Family_16_LS
ancom_Family_16S_LS = ancombc2(data = LSpool_16S.ps, assay_name = "counts", tax_level = "Family",
                             fix_formula = "StatsPool + Number.in.Pool", rand_formula = NULL,
                             p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                             group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                             global = TRUE, pairwise = TRUE,
                             dunnet = TRUE, trend = FALSE,
                             iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                             em_control = list(tol = 1e-5, max_iter = 100),
                             lme_control = lme4::lmerControl(),
                             mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                             trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                         nrow = 2,
                                                                         byrow = TRUE),matrix(c(-1, 0, 1, -1),
                                                                                              nrow = 2,
                                                                                              byrow = TRUE)),
                                                  node = list(2, 2),
                                                  solver = "ECOS",
                                                  B = 100))


res_pair_16S_LS = ancom_Family_16S_LS$res_pair
glob_16S_LS = ancom_Family_16S_LS$res_global

dunnet_16S_LS = ancom_Family_16S_LS$res_dunn
BCabundancetable_16S_LS


df_fig_dunnet_16S_LS1 = dunnet_16S_LS %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc3,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_16S_LS2 = dunnet_16S_LS %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_16S_LS = df_fig_dunnet_16S_LS1 %>%
  dplyr::left_join(df_fig_dunnet_16S_LS2, by=c("taxon","group"))

df_fig_dunnet_16S_LS$group = recode(df_fig_dunnet_16S_LS$group,
                                    `lfc1`= "Pools of 3",
                                    `lfc2`= "Pools of 6",
                                    `lfc3`= "Pools of 12")

df_fig_dunnet_16S_LS$group =factor(df_fig_dunnet_16S_LS$group,
                                   levels = c("Pools of 3",
                                              "Pools of 6",
                                              "Pools of 12"))

lo_16S_LS_dunn = floor(min(df_fig_dunnet_16S_LS$value))
up_16S_LS_dunn = ceiling(max(df_fig_dunnet_16S_LS$value))
mid_16S_LS_dunn = (lo_16S_LS_dunn+up_16S_LS_dunn)/2

fig_dunn_16S_LS = df_fig_dunnet_16S_LS %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_LS_dunn,up_16S_LS_dunn),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = "none")+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

###Paired
df_fig_res_pair_16S_LS1 = res_pair_16S_LS %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1|
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0),
                lfc4=ifelse(diff_Number.in.Pool6_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool6_Number.in.Pool3,2),0),
                lfc5=ifelse(diff_Number.in.Pool12_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool3, 2),0),
                lfc6=ifelse(diff_Number.in.Pool12_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool6, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc6,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_LS2 = res_pair_16S_LS %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12 ==1 |
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black"),
                lfc4= ifelse(passed_ss_Number.in.Pool6_Number.in.Pool3 ==1 & diff_Number.in.Pool6_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc5= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool3 ==1 & diff_Number.in.Pool12_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc6= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool6 ==1 & diff_Number.in.Pool12_Number.in.Pool6 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_LS = df_fig_res_pair_16S_LS1 %>%
  dplyr::left_join(df_fig_res_pair_16S_LS2, by=c("taxon","group"))

df_fig_res_pair_16S_LS$group = recode(df_fig_res_pair_16S_LS$group,
                                      `lfc1`= "3 vs 1",
                                      `lfc2`= "6 vs 1",
                                      `lfc3`= "12 vs 1",
                                      `lfc4`= "6 vs 3",
                                      `lfc5`= "12 vs 3",
                                      `lfc6`= "12 vs 6")

df_fig_res_pair_16S_LS$group =factor(df_fig_res_pair_16S_LS$group,
                                     levels = c("3 vs 1",
                                                "6 vs 1",
                                                "12 vs 1",
                                                "6 vs 3",
                                                "12 vs 3",
                                                "12 vs 6"))

lo_16S_LS_pair = floor(min(df_fig_res_pair_16S_LS$value))
up_16S_LS_pair = ceiling(max(df_fig_res_pair_16S_LS$value))
mid_16S_LS_pair = (lo_16S_LS_pair+up_16S_LS_pair)/2

fig_pair_16S_LS = df_fig_res_pair_16S_LS %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_LS_pair,up_16S_LS_pair),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_HS_16S = res_pair_HS_16S %>%
  dplyr::select(taxon, contains("Number.in.Pool"))


df_LS_fig_glob = glob_LS_16S%>%
  filter(passed_ss ==1 & diff_abn==1)

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_LS_16S = res_pair_LS_16S %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

# Extract the pairwise comparison results
df_fig_LS_16S = df_LS_16S %>%
  filter((diff_Number.in.Pool3 == 1 & passed_ss_Number.in.Pool3==1)| (diff_Number.in.Pool6 == 1 & passed_ss_Number.in.Pool6==1)|
           (diff_Number.in.Pool12 == 1 & passed_ss_Number.in.Pool12==1)| (diff_Number.in.Pool6_Number.in.Pool3  == 1 & passed_ss_Number.in.Pool6_Number.in.Pool3==1)|
           (diff_Number.in.Pool12_Number.in.Pool3  == 1 & passed_ss_Number.in.Pool12_Number.in.Pool3==1)| (diff_Number.in.Pool12_Number.in.Pool6 == 1 & passed_ss_Number.in.Pool12_Number.in.Pool6==1 )) %>%
  mutate(lfc_Number.in.Pool3 = ifelse(diff_Number.in.Pool3 == 1,  # 3v1
                                      lfc_Number.in.Pool3, 0),
         lfc_Number.in.Pool6 = ifelse(diff_Number.in.Pool6 == 1,  # 6v1
                                      lfc_Number.in.Pool6, 0),
         lfc_Number.in.Pool12 = ifelse(diff_Number.in.Pool12 == 1,  # 12v1
                                       lfc_Number.in.Pool12, 0),
         lfc_Number.in.Pool6_Number.in.Pool3 = ifelse(diff_Number.in.Pool6_Number.in.Pool3 == 1, #6v3
                                                      lfc_Number.in.Pool6_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool3 = ifelse(diff_Number.in.Pool12_Number.in.Pool3 == 1,  # 12v3
                                                       lfc_Number.in.Pool12_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool6 = ifelse(diff_Number.in.Pool12_Number.in.Pool6 == 1,  # 12v6
                                                       lfc_Number.in.Pool12_Number.in.Pool6, 0),
  ) %>%
  transmute(taxon, 
            `3v1` = round(lfc_Number.in.Pool3, 2),
            `6v1` = round(lfc_Number.in.Pool6, 2),
            `12v1` = round(lfc_Number.in.Pool12, 2),
            `6v3` = round(lfc_Number.in.Pool6_Number.in.Pool3, 2),
            `12v3` = round(lfc_Number.in.Pool12_Number.in.Pool3, 2),
            `12v6` = round(lfc_Number.in.Pool12_Number.in.Pool6,2)) %>%
  pivot_longer(cols = `3v1`:`6v1`:`12v1`:`6v3`:`12v3`:`12v6`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Change order of X-value
df_fig_LS_16S
df_fig_LS_16S$group <- factor(df_fig_LS_16S$group, levels = c("3v1","6v1","12v1","6v3","12v3","12v6"))




# Plot heatmap
lo = floor(min(df_fig_LS_16S$value))
up = ceiling(max(df_fig_LS_16S$value))
mid = (lo + up)/2


fig_LS_16S = df_fig_LS_16S %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by combined order") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_LS_16S

unique(df_LS_fig_glob$taxon)

df_fig_LS_16S_filtered <- subset(df_fig_LS_16S, taxon=="unclassified Bacteria" | taxon==  "unclassified Unassigned"| taxon==    "Family_XI"| taxon== "Eubacteriaceae" | 
                                   taxon== "Erysipelatoclostridiaceae"| taxon==  "Streptococcaceae" | taxon== "Acidaminococcaceae"| taxon==  "Propionibacteriaceae" | 
                                   taxon== "Microbacteriaceae"| taxon== "Neisseriaceae"| taxon== "Pasteurellaceae"| taxon== "Succinivibrionaceae"|
                                   taxon== "Pseudomonadaceae"| taxon==  "Prevotellaceae" | taxon==  "Lachnospiraceae")

#filtering
lo_LS_16S_filtered = floor(min(df_fig_LS_16S_filtered$value))
up_LS_16S_filtered = ceiling(max(df_fig_LS_16S_filtered$value))
mid_LS_16S_filtered = (lo_LS_16S_filtered + up_LS_16S_filtered)/2


fig_LS_16S_filtered = df_fig_LS_16S_filtered %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo_LS_16S_filtered, up_LS_16S_filtered),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by combined order") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_LS_16S_filtered


########
HRpool_16S.ps <- subset_samples(data_16S_50k, PrevH=="Y"&PoolingType!="Subset")
sample_data(HRpool_16S.ps)$HRPool


sample_data(HRpool_16S.ps)$Number.in.Pool <- factor(sample_data(HRpool_16S.ps)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(HRpool_16S.ps)

#Running model
ancom_Family_16S_HR = ancombc2(data = HRpool_16S.ps, assay_name = "counts", tax_level = "Family",
                             fix_formula = "Number.in.Pool", rand_formula = NULL,
                             p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                             group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                             global = TRUE, pairwise = TRUE,
                             dunnet = TRUE, trend = FALSE,
                             iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                             em_control = list(tol = 1e-5, max_iter = 100),
                             lme_control = lme4::lmerControl(),
                             mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                             trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                         nrow = 2,
                                                                         byrow = TRUE),matrix(c(-1, 0, 1, -1),
                                                                                              nrow = 2,
                                                                                              byrow = TRUE)),
                                                  node = list(2, 2),
                                                  solver = "ECOS",
                                                  B = 100))


res_pair_16S_HR = ancom_Family_16S_HR$res_pair
glob_HR_16S = ancom_Family_16S_HR$res_global
dunn_16S_HR = ancom_Family_16S_HR$res_dunn

df_HR_16S_fig_glob = glob_HR_16S%>%
  filter(passed_ss ==1 & diff_abn==1)

unique(df_HR_16S_fig_glob$taxon)


#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons

                
df_fig_dunnet_16S_HR1 = dunn_16S_HR %>%
                  dplyr::filter(diff_Number.in.Pool3==1 |
                                  diff_Number.in.Pool6==1|
                                  diff_Number.in.Pool12==1)%>%
                  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                                            round(lfc_Number.in.Pool3,2),0),
                                lfc2=ifelse(diff_Number.in.Pool6==1,
                                            round(lfc_Number.in.Pool6, 2),0),
                                lfc3=ifelse(diff_Number.in.Pool12==1,
                                            round(lfc_Number.in.Pool12, 2),0)) %>%
                  tidyr::pivot_longer(cols=lfc1:lfc3,
                                      names_to = "group",values_to = "value") %>%
                  dplyr::arrange(taxon)

# Extract the pairwise comparison results
df_fig_dunnet_16S_HR2 = dunnet_16S_HR %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_16S_HR = df_fig_dunnet_16S_HR1 %>%
  dplyr::left_join(df_fig_dunnet_16S_HR2, by=c("taxon","group"))

df_fig_dunnet_16S_HR$group = recode(df_fig_dunnet_16S_HR$group,
                                    `lfc1`= "Pools of 3",
                                    `lfc2`= "Pools of 6",
                                    `lfc3`= "Pools of 12")

df_fig_dunnet_16S_HR$group =factor(df_fig_dunnet_16S_HR$group,
                                   levels = c("Pools of 3",
                                              "Pools of 6",
                                              "Pools of 12"))

lo_16S_HR_dunn = floor(min(df_fig_dunnet_16S_HR$value))
up_16S_HR_dunn = ceiling(max(df_fig_dunnet_16S_HR$value))
mid_16S_HR_dunn = (lo_16S_HR_dunn+up_16S_HR_dunn)/2

fig_dunn_16S_HR = df_fig_dunnet_16S_HR %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_HR_dunn,up_16S_HR_dunn),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

# Change order of X-value
df_fig_HR_16S
df_fig_HR_16S$group <- factor(df_fig_HR_16S$group, levels = c("3v1","6v1","12v1","6v3","12v3","12v6"))

df_fig_HR_16S_filtered <- subset(df_fig_HR_16S, taxon== "unclassified Unassigned" |taxon== "Acholeplasmataceae"|taxon=="RF39"|taxon=="Mycoplasmataceae"|
                                   taxon=="Lactobacillaceae"|taxon=="Streptococcaceae"|taxon== "Selenomonadaceae"|taxon== "Campylobacteraceae"|
                                   taxon=="Microbacteriaceae"|taxon=="Neisseriaceae"|taxon=="Comamonadaceae"|taxon=="Moraxellaceae"|
                                   taxon=="Succinivibrionaceae"|taxon=="Methanobacteriaceae"|taxon=="Spirochaetaceae"|taxon=="Chitinophagaceae")




# Plot heatmap
lo = floor(min(df_fig_HR_16S$value))
up = ceiling(max(df_fig_HR_16S$value))
mid = (lo + up)/2


fig_HR_16S = df_fig_HR_16S %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 

  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 3) +
  labs(x = NULL, y = NULL, title = "LogFC by Pool Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_HR_16S

lo_HR_16S_filtered = floor (min(df_fig_HR_16S_filtered$value))
up_HR_16S_filtered =  ceiling(max(df_fig_HR_16S_filtered$value))
mid_HS_16S_filtered=(lo_HR_16S_filtered + up_HR_16S_filtered)/2

fig_HR_16S_filtered = df_fig_HR_16S_filtered %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo_HR_16S_filtered, up_HR_16S_filtered),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 3) +
  labs(x = NULL, y = NULL, title = "LogFC by Pool Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_HR_16S_filtered




##LOW RANDOM
LRpool_16S.ps <- subset_samples(data_16S_50k, PrevL=="Y"&PoolingType!="Subset")
View(sample_data(LRpool_16S.ps))

sample_data(LRpool_16S.ps)$LSPool

sample_data(LRpool_16S.ps)$StatsPool<-c("A","B","C","D","E","F","A","B","C","D","E","F","A","B","C","D","E","F","D","F","E","E","F","A","B","D","F","D","E",
                                        "C","C","F","D","F","C","B","A","E","A","A","E","F","F","F","A","D","D","A","E","F","B","A","B","E","A","E","C","A",
                                        "F","B","C","A","C","C","B","C","B","D","C","D","E","E","A","B","D","B","B","F","F","A","B","D","E","C","B","D","C","C")
sample_data(LRpool_16S.ps)$StatsPool <- factor(sample_data(LRpool_16S.ps)$StatsPool, levels = c("A","B","C","D","E","F"))

sample_data(LRpool_16S.ps)$Number.in.Pool <- factor(sample_data(LRpool_16S.ps)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(LRpool_16S.ps)

#Running model
ancom_Family_16S_LR = ancombc2(data = LRpool_16S.ps, assay_name = "counts", tax_level = "Family",
                              fix_formula = "Number.in.Pool", rand_formula = NULL,
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = TRUE,
                              dunnet = TRUE, trend = FALSE,
                              iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                          nrow = 2,
                                                                          byrow = TRUE),matrix(c(-1, 0, 1, -1),
                                                                                               nrow = 2,
                                                                                               byrow = TRUE)),
                                                   node = list(2, 2),
                                                   solver = "ECOS",
                                                   B = 100))


res_pair_16S_LR = ancom_Family_16_LS$res_pair
res_glob_LR_16S = ancom_Family_16_LS$res_global
dunnet_16S_LR = ancom_Family_16S_LR$res_dunn
BCabundancetable_16S_LR = ancom_Family_16S_LR$bias_correct_log_table

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_LR_16S = res_pair_LR_16S %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

# Extract the pairwise comparison results
df_fig_dunnet_16S_LR1 = dunnet_16S_LR %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc3,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

# Extract the pairwise comparison results
df_fig_dunnet_16S_LR2 = dunnet_16S_LR %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_16S_LR = df_fig_dunnet_16S_LR1 %>%
  dplyr::left_join(df_fig_dunnet_16S_LR2, by=c("taxon","group"))

df_fig_dunnet_16S_LR$group = recode(df_fig_dunnet_16S_LR$group,
                                    `lfc1`= "Pools of 3",
                                    `lfc2`= "Pools of 6",
                                    `lfc3`= "Pools of 12")

df_fig_dunnet_16S_LR$group =factor(df_fig_dunnet_16S_LR$group,
                                   levels = c("Pools of 3",
                                              "Pools of 6",
                                              "Pools of 12"))

lo_16S_LR_dunn = floor(min(df_fig_dunnet_16S_LR$value))
up_16S_LR_dunn = ceiling(max(df_fig_dunnet_16S_LR$value))
mid_16S_LR_dunn = (lo_16S_LR_dunn+up_16S_LR_dunn)/2

fig_dunn_16S_LR = df_fig_dunnet_16S_LR %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_LR_dunn,up_16S_LR_dunn),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

fig_dunn_16S_LR

###RAND
df_fig_res_pair_16S_HR1 = res_pair_16S_HR %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1|
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0),
                lfc4=ifelse(diff_Number.in.Pool6_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool6_Number.in.Pool3,2),0),
                lfc5=ifelse(diff_Number.in.Pool12_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool3, 2),0),
                lfc6=ifelse(diff_Number.in.Pool12_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool6, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc6,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_HR2 = res_pair_16S_HR %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12 ==1 |
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black"),
                lfc4= ifelse(passed_ss_Number.in.Pool6_Number.in.Pool3 ==1 & diff_Number.in.Pool6_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc5= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool3 ==1 & diff_Number.in.Pool12_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc6= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool6 ==1 & diff_Number.in.Pool12_Number.in.Pool6 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_HR = df_fig_res_pair_16S_HR1 %>%
  dplyr::left_join(df_fig_res_pair_16S_HR2, by=c("taxon","group"))

df_fig_res_pair_16S_HR$group = recode(df_fig_res_pair_16S_HR$group,
                                      `lfc1`= "3 vs 1",
                                      `lfc2`= "6 vs 1",
                                      `lfc3`= "12 vs 1",
                                      `lfc4`= "6 vs 3",
                                      `lfc5`= "12 vs 3",
                                      `lfc6`= "12 vs 6")

df_fig_res_pair_16S_HR$group =factor(df_fig_res_pair_16S_HR$group,
                                     levels = c("3 vs 1",
                                                "6 vs 1",
                                                "12 vs 1",
                                                "6 vs 3",
                                                "12 vs 3",
                                                "12 vs 6"))

lo_16S_HR_pair = floor(min(df_fig_res_pair_16S_HR$value))
up_16S_HR_pair = ceiling(max(df_fig_res_pair_16S_HR$value))
mid_16S_HR_pair = (lo_16S_HR_pair+up_16S_HR_pair)/2

fig_pair_16S_HR = df_fig_res_pair_16S_HR %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_HR_pair,up_16S_HR_pair),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

df_fig_res_pair_16S_LR1 = res_pair_16S_LR %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1|
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0),
                lfc4=ifelse(diff_Number.in.Pool6_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool6_Number.in.Pool3,2),0),
                lfc5=ifelse(diff_Number.in.Pool12_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool3, 2),0),
                lfc6=ifelse(diff_Number.in.Pool12_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool12_Number.in.Pool6, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc6,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_LR2 = res_pair_16S_LR %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12 ==1 |
                  diff_Number.in.Pool6_Number.in.Pool3==1 |
                  diff_Number.in.Pool12_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool12_Number.in.Pool6 ==1) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black"),
                lfc4= ifelse(passed_ss_Number.in.Pool6_Number.in.Pool3 ==1 & diff_Number.in.Pool6_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc5= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool3 ==1 & diff_Number.in.Pool12_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc6= ifelse(passed_ss_Number.in.Pool12_Number.in.Pool6 ==1 & diff_Number.in.Pool12_Number.in.Pool6 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_res_pair_16S_LR = df_fig_res_pair_16S_LR1 %>%
  dplyr::left_join(df_fig_res_pair_16S_LR2, by=c("taxon","group"))

df_fig_res_pair_16S_LR$group = recode(df_fig_res_pair_16S_LR$group,
                                      `lfc1`= "3 vs 1",
                                      `lfc2`= "6 vs 1",
                                      `lfc3`= "12 vs 1",
                                      `lfc4`= "6 vs 3",
                                      `lfc5`= "12 vs 3",
                                      `lfc6`= "12 vs 6")

df_fig_res_pair_16S_LR$group =factor(df_fig_res_pair_16S_LR$group,
                                     levels = c("3 vs 1",
                                                "6 vs 1",
                                                "12 vs 1",
                                                "6 vs 3",
                                                "12 vs 3",
                                                "12 vs 6"))

lo_16S_LR_pair = floor(min(df_fig_res_pair_16S_LR$value))
up_16S_LR_pair = ceiling(max(df_fig_res_pair_16S_LR$value))
mid_16S_LR_pair = (lo_16S_LR_pair+up_16S_LR_pair)/2

fig_pair_16S_LR = df_fig_res_pair_16S_LR %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_16S_LR_pair,up_16S_LR_pair),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(fig_pair_16S_HS, fig_pair_16S_HR, fig_pair_16S_LS, fig_pair_16S_LR, nrow = 2, ncol = 2, labels = c("HIGH-SUB","HIGH-RAND","LOW-SUB","LOW-RAND"),heights = c(1.5,1), common.legend = T)


####AMR####
##High
high_metaSNV.ps <- subset_samples(ARG_noSNP_ANCOMBC, PrevH=="Y")
any(sample_sums(high_metaSNV.ps)==0)
any(taxa_sums(high_metaSNV.ps)==0)
high_metaSNV.ps <- prune_taxa(taxa_sums(high_metaSNV.ps)>0, high_metaSNV.ps)

sample_data(high_metaSNV.ps)

View(sample_data(high_metaSNV.ps))

sample_data(high_metaSNV.ps)$StatsPool<-c("A","B",
                                        "C","C","C","C","C","C","C","C","C","C","C","C","C",
                                        "D","D","D","D","D","D","D","D","D","D","D","D","D",
                                        "E","F","A","B","C","D","E","F","A","B","C","D","E","F")
sample_data(high_metaSNV.ps)$StatsPool <- factor(sample_data(high_metaSNV.ps)$StatsPool, levels = c("A","B","C","D","E","F"))
tax_table(high_metaSNV.ps)

##Low
low_metaSNV.ps <- subset_samples(ARG_noSNP_ANCOMBC, Prevalence =="Low")
sample_data(low_metaSNV.ps)
low_metaSNV.ps <- subset_samples(low_metaSNV.ps, SampleType=="Pool")
any(sample_sums(low_metaSNV.ps)==0)
any(taxa_sums(low_metaSNV.ps)==0)
low_metaSNV.ps <- prune_taxa(taxa_sums(low_metaSNV.ps)>0, low_metaSNV.ps)
View(sample_data(low_metaSNV.ps))
sample_data(low_metaSNV.ps)$StatsPool <- c("A","B","C","D","E","F","A","B","C","D","E","F","A","B","C","D","E","F")
tax_table(low_metaSNV.ps)

##Let's go!
ancom_output_high_amr_1 = ancombc2(data = high_metaSNV.ps, assay_name = "counts", tax_level = "Species",
                               fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = TRUE, pairwise = TRUE,
                               dunnet = TRUE, trend = FALSE,
                               iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),
                                                    node = list(2, 2),solver = "ECOS",B = 100))

##Extracting pairwise comparisons
res_pair_high_amr = ancom_output_high_amr_1$res_pair
BCabundance_amr_hs =ancom_output_high_amr_1$bias_correct_log_table

dunnet_amr_HS = ancom_output_high_amr_1$res_dunn
unique(dunnet_amr_HS$diff_Number.in.Pool3)
unique(dunnet_amr_HS$diff_Number.in.Pool6)
unique(dunnet_amr_HS$diff_Number.in.Pool12)

df_fig_dunnet_amr_HS1 = dunnet_amr_HS %>%
  dplyr::filter(diff_Number.in.Pool3==1 |
                  diff_Number.in.Pool6==1|
                  diff_Number.in.Pool12==1)%>%
  dplyr::mutate(lfc1=ifelse(diff_Number.in.Pool3==1,
                            round(lfc_Number.in.Pool3,2),0),
                lfc2=ifelse(diff_Number.in.Pool6==1,
                            round(lfc_Number.in.Pool6, 2),0),
                lfc3=ifelse(diff_Number.in.Pool12==1,
                            round(lfc_Number.in.Pool12, 2),0)) %>%
  tidyr::pivot_longer(cols=lfc1:lfc3,
                      names_to = "group",values_to = "value") %>%
  dplyr::arrange(taxon)

diff_amr_HS_3 = dunnet_amr_HS %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_amr_HS2 = dunnet_amr_HS %>%
  dplyr::filter(diff_Number.in.Pool3 ==1 |
                  diff_Number.in.Pool6 ==1 |
                  diff_Number.in.Pool12) %>%
  dplyr::mutate(lfc1= ifelse(passed_ss_Number.in.Pool3 ==1 & diff_Number.in.Pool3 ==1,
                             "aquamarine3","black"),
                lfc2= ifelse(passed_ss_Number.in.Pool6 ==1 & diff_Number.in.Pool6 ==1,
                             "aquamarine3","black"),
                lfc3= ifelse(passed_ss_Number.in.Pool12 ==1 & diff_Number.in.Pool12 ==1,
                             "aquamarine3","black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3,
                      names_to = "group",values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_dunnet_amr_HS = df_fig_dunnet_amr_HS1 %>%
  dplyr::left_join(df_fig_dunnet_amr_HS2, by=c("taxon","group"))

df_fig_dunnet_amr_HS$group = recode(df_fig_dunnet_amr_HS$group,
                                    `lfc1`= "Pools of 3",
                                    `lfc2`= "Pools of 6",
                                    `lfc3`= "Pools of 12")

df_fig_dunnet_amr_HS$group =factor(df_fig_dunnet_amr_HS$group,
                                   levels = c("Pools of 3",
                                              "Pools of 6",
                                              "Pools of 12"))

lo_amr_HS_dunn = floor(min(df_fig_dunnet_amr_HS$value))
up_amr_HS_dunn = ceiling(max(df_fig_dunnet_amr_HS$value))
mid_amr_HS_dunn = (lo_amr_HS_dunn+up_amr_HS_dunn)/2

fig_dunn_amr_HS = df_fig_dunnet_amr_HS %>%
  ggplot(aes(x=group, y=taxon, fill=value)) +
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue",high = "red",mid = "white",
                       na.value = "white", midpoint = 0, limit=c(lo_amr_HS_dunn,up_amr_HS_dunn),
                       name = NULL)+
  geom_text(aes(group,taxon,label=value, color=color),size=4)+
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

##Things to define
###lfc_Number.in.Pool3 lfc_Number.in.Pool6 lfc_Number.in.Pool12 lfc_Number.in.Pool6_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool6

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_pool_high_amr = res_pair_high_amr %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

df_pool_high_amr$lfc_Number.in.Pool3

# Extract the pairwise comparison results
df_fig_pool_high_amr = df_pool_high_amr %>%
  filter(diff_Number.in.Pool3 == 1 | diff_Number.in.Pool6 == 1 | diff_Number.in.Pool12 ==1
         | diff_Number.in.Pool6_Number.in.Pool3  == 1 | diff_Number.in.Pool12_Number.in.Pool3  == 1 
         | diff_Number.in.Pool12_Number.in.Pool6 == 1 ) %>%
  mutate(lfc_Number.in.Pool3 = ifelse(diff_Number.in.Pool3 == 1,  # This compares the 2nd time point to Start_CONV
                                      lfc_Number.in.Pool3, 0),
         lfc_Number.in.Pool6 = ifelse(diff_Number.in.Pool6== 1, 
                                      lfc_Number.in.Pool6, 0),
         lfc_Number.in.Pool12 = ifelse(diff_Number.in.Pool12 == 1, 
                                       lfc_Number.in.Pool12, 0),
         lfc_Number.in.Pool6_Number.in.Pool3 = ifelse(diff_Number.in.Pool6_Number.in.Pool3== 1, 
                                                      lfc_Number.in.Pool6_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool3 = ifelse(diff_Number.in.Pool12_Number.in.Pool3== 1, 
                                                       lfc_Number.in.Pool12_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool6 = ifelse(diff_Number.in.Pool12_Number.in.Pool6 == 1,
                                                       lfc_Number.in.Pool12_Number.in.Pool6, 0)
  ) %>%
  transmute(taxon, 
            `3 v 1` = round(lfc_Number.in.Pool3, 2),
            `6 v 1` = round(lfc_Number.in.Pool6, 2),
            `12 v 1` = round(lfc_Number.in.Pool12, 2),
            `6 v 3` = round(lfc_Number.in.Pool6_Number.in.Pool3, 2),
            `12 v 3` = round(lfc_Number.in.Pool12_Number.in.Pool3, 2),
            `12 v 6` = round(lfc_Number.in.Pool12_Number.in.Pool6, 2)) %>%
  pivot_longer(cols = `3 v 1`:`6 v 1`:`12 v 1`:`6 v 3`:`12 v 3`:`12 v 6`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Change order of X-value
df_fig_pool_high_amr
df_fig_pool_high$group <- factor(df_fig_pool_high$group, levels = c("3 v 1","6 v 1","12 v 1","6 v 3","12 v 3","12 v 6"))




# Plot heatmap
lo = floor(min(df_fig_pool_high$value))
up = ceiling(max(df_fig_pool_high$value))
mid = (lo + up)/2


fig_ARG_high = df_fig_pool_high %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by Pool Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_ARG_high


##Let's go!
ancom_output_low_1 = ancombc2(data = low_metaSNV.ps, assay_name = "counts", tax_level = "Species",
                               fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = TRUE, pairwise = TRUE,
                               dunnet = TRUE, trend = FALSE,
                               iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),
                                                    node = list(2, 2),solver = "ECOS",B = 100))

##Extracting pairwise comparisons
res_pair_low = ancom_output_low_1$res_pair

##Things to define
###lfc_Number.in.Pool3 lfc_Number.in.Pool6 lfc_Number.in.Pool12 lfc_Number.in.Pool6_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool6

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_pool_low = res_pair_low %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

df_pool_low$lfc_Number.in.Pool6
unique(df_pool_low$diff_Number.in.Pool12_Number.in.Pool6)

# Extract the pairwise comparison results
df_fig_pool_low = df_pool_low %>%
  filter(diff_Number.in.Pool6 == 1 | diff_Number.in.Pool12 ==1
         | diff_Number.in.Pool12_Number.in.Pool6 == 1 ) %>%
  mutate(lfc_Number.in.Pool6 = ifelse(diff_Number.in.Pool6== 1, 
                                      lfc_Number.in.Pool6, 0),
         lfc_Number.in.Pool12 = ifelse(diff_Number.in.Pool12 == 1, 
                                       lfc_Number.in.Pool12, 0),
         lfc_Number.in.Pool12_Number.in.Pool6 = ifelse(diff_Number.in.Pool12_Number.in.Pool6 == 1,
                                                       lfc_Number.in.Pool12_Number.in.Pool6, 0)
  ) %>%
  transmute(taxon, 
            `6 v 3` = round(lfc_Number.in.Pool6, 2),
            `12 v 3` = round(lfc_Number.in.Pool12, 2),
            `12 v 6` = round(lfc_Number.in.Pool12_Number.in.Pool6, 2)) %>%
  pivot_longer(cols = `6 v 3`:`12 v 3`:`12 v 6`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Change order of X-value
df_fig_pool_low
df_fig_pool_low$group <- factor(df_fig_pool_low$group, levels = c("6 v 3","12 v 3","12 v 6"))




# Plot heatmap
lo = floor(min(df_fig_pool_low$value))
up = ceiling(max(df_fig_pool_low$value))
mid = (lo + up)/2


fig_ARG_low = df_fig_pool_low %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", low = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by Pool Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_ARG_low

###MHDATA####
PSVdata


high_MH.ps <- subset_samples(PSVdata, PrevH=="Y")
any(sample_sums(high_MH.ps)==0)
any(taxa_sums(high_MH.ps)==0)

sample_data(high_MH.ps)
sample_names(high_MH.ps)

sample_data(high_MH.ps)$StatsPool<-c("F","D","B","C","D","D","B",
                                          "C","E","A","C","D","C","D",
                                          "B","D","D","C","F","A","F",
                                          "C","C","C","D","D","D","C",
                                          "C","C","C","C","D","D","C",
                                          "C","D","E","E","A","D","D")
sample_data(high_MH.ps)$StatsPool <- factor(sample_data(high_MH.ps)$StatsPool, levels = c("A","B","C","D","E","F"))
sample_data(high_MH.ps)$Number.in.Pool <-factor(sample_data(high_MH.ps)$Number.in.Pool, levels = c(1,3,6,12))


##Low
low_MH.ps <- subset_samples(PSVdata, Prevalence =="Low")
sample_data(low_MH.ps)

any(sample_sums(low_MH.ps)==0)
any(taxa_sums(low_MH.ps)==0)
sample_names(low_MH.ps)

sample_data(low_MH.ps)$StatsPool <- c("E","A","D","C","D","D","F","C","B","E","B",
                                           "F","A","C","E","F","A","B")
sample_data(low_MH.ps)$StatsPool <- factor(sample_data(low_MH.ps)$StatsPool, levels = c("A","B","C","D","E","F"))
sample_data(low_MH.ps)$Number.in.Pool <-factor(sample_data(low_MH.ps)$Number.in.Pool, levels = c(1,3,6,12))

high_MH_ancom <- high_MH.ps
colnames(tax_table(high_MH_ancom)) <- c("Family","Genus","Species","All")

##Let's go!
ancom_output_high_MH = ancombc2(data = high_MH_ancom, assay_name = "counts", tax_level = "Species",
                               fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                               group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                               global = TRUE, pairwise = TRUE,
                               dunnet = TRUE, trend = FALSE,
                               iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                               trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),
                                                                    matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),
                                                    node = list(2, 2),solver = "ECOS",B = 100))

##Extracting pairwise comparisons
res_pair_high_MH = ancom_output_high_MH$res_pair

##Things to define
###lfc_Number.in.Pool3 lfc_Number.in.Pool6 lfc_Number.in.Pool12 lfc_Number.in.Pool6_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool6

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_pool_high_MH = res_pair_high_MH %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 



# Extract the pairwise comparison results
df_fig_pool_high_MH = df_pool_high_MH %>%
  filter((diff_Number.in.Pool3 == 1 & passed_ss_Number.in.Pool3==1)| (diff_Number.in.Pool6 == 1 & passed_ss_Number.in.Pool6==1)|
           (diff_Number.in.Pool12 == 1 & passed_ss_Number.in.Pool12==1)| (diff_Number.in.Pool6_Number.in.Pool3  == 1 & passed_ss_Number.in.Pool6_Number.in.Pool3==1)|
           (diff_Number.in.Pool12_Number.in.Pool3  == 1 & passed_ss_Number.in.Pool12_Number.in.Pool3==1)| (diff_Number.in.Pool12_Number.in.Pool6 == 1 & passed_ss_Number.in.Pool12_Number.in.Pool6==1 )) %>%
  mutate(lfc_Number.in.Pool3 = ifelse(diff_Number.in.Pool3 == 1,  # This compares the 2nd time point to Start_CONV
                                      lfc_Number.in.Pool3, 0),
         lfc_Number.in.Pool6 = ifelse(diff_Number.in.Pool6== 1, 
                                      lfc_Number.in.Pool6, 0),
         lfc_Number.in.Pool12 = ifelse(diff_Number.in.Pool12 == 1, 
                                       lfc_Number.in.Pool12, 0),
         lfc_Number.in.Pool6_Number.in.Pool3 = ifelse(diff_Number.in.Pool6_Number.in.Pool3== 1, 
                                                      lfc_Number.in.Pool6_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool3 = ifelse(diff_Number.in.Pool12_Number.in.Pool3== 1, 
                                                       lfc_Number.in.Pool12_Number.in.Pool3, 0),
         lfc_Number.in.Pool12_Number.in.Pool6 = ifelse(diff_Number.in.Pool12_Number.in.Pool6 == 1,
                                                       lfc_Number.in.Pool12_Number.in.Pool6, 0)
  ) %>%
  transmute(taxon, 
            `3 v 1` = round(lfc_Number.in.Pool3, 2),
            `6 v 1` = round(lfc_Number.in.Pool6, 2),
            `12 v 1` = round(lfc_Number.in.Pool12, 2),
            `6 v 3` = round(lfc_Number.in.Pool6_Number.in.Pool3, 2),
            `12 v 3` = round(lfc_Number.in.Pool12_Number.in.Pool3, 2),
            `12 v 6` = round(lfc_Number.in.Pool12_Number.in.Pool6, 2)) %>%
  pivot_longer(cols = `3 v 1`:`6 v 1`:`12 v 1`:`6 v 3`:`12 v 3`:`12 v 6`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Change order of X-value
df_fig_pool_high_MH
df_fig_pool_high$group <- factor(df_fig_pool_high$group, levels = c("3 v 1","6 v 1","12 v 1","6 v 3","12 v 3","12 v 6"))


##Let's go!
ancom_output_low_MH = ancombc2(data = low_MH.ps, assay_name = "counts", tax_level = "Species",
                              fix_formula = "Number.in.Pool+StatsPool", rand_formula = NULL,
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                              group = "Number.in.Pool", struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = TRUE, pairwise = TRUE,
                              dunnet = TRUE, trend = FALSE,
                              iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                              trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),
                                                                   matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),
                                                   node = list(2, 2),solver = "ECOS",B = 100))

##Extracting pairwise comparisons
res_pair_low = ancom_output_low_1$res_pair

##Things to define
###lfc_Number.in.Pool3 lfc_Number.in.Pool6 lfc_Number.in.Pool12 lfc_Number.in.Pool6_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool3 lfc_Number.in.Pool12_Number.in.Pool6

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_pool_low = res_pair_low %>%
  dplyr::select(taxon, contains("Number.in.Pool")) 

df_pool_low$lfc_Number.in.Pool6
unique(df_pool_low$diff_Number.in.Pool12_Number.in.Pool6)

# Extract the pairwise comparison results
df_fig_pool_low = df_pool_low %>%
  filter(diff_Number.in.Pool6 == 1 | diff_Number.in.Pool12 ==1
         | diff_Number.in.Pool12_Number.in.Pool6 == 1 ) %>%
  mutate(lfc_Number.in.Pool6 = ifelse(diff_Number.in.Pool6== 1, 
                                      lfc_Number.in.Pool6, 0),
         lfc_Number.in.Pool12 = ifelse(diff_Number.in.Pool12 == 1, 
                                       lfc_Number.in.Pool12, 0),
         lfc_Number.in.Pool12_Number.in.Pool6 = ifelse(diff_Number.in.Pool12_Number.in.Pool6 == 1,
                                                       lfc_Number.in.Pool12_Number.in.Pool6, 0)
  ) %>%
  transmute(taxon, 
            `6 v 3` = round(lfc_Number.in.Pool6, 2),
            `12 v 3` = round(lfc_Number.in.Pool12, 2),
            `12 v 6` = round(lfc_Number.in.Pool12_Number.in.Pool6, 2)) %>%
  pivot_longer(cols = `6 v 3`:`12 v 3`:`12 v 6`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Change order of X-value
df_fig_pool_low
df_fig_pool_low$group <- factor(df_fig_pool_low$group, levels = c("6 v 3","12 v 3","12 v 6"))




# Plot heatmap
lo = floor(min(df_fig_pool_low$value))
up = ceiling(max(df_fig_pool_low$value))
mid = (lo + up)/2


fig_ARG_low = df_fig_pool_low %>%
  ggplot(aes(x = group, y = reorder(taxon, value), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", low = "red", mid = "white",
                       midpoint = 0, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "LogFC by Pool Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_ARG_low

