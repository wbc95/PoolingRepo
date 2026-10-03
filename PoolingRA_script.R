#load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr);
library(ggplot2); library(btools); library(randomcoloR); library(cowplot);
library(pairwiseAdonis); library(picante); library(gridExtra); library(plyr);
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(rstatix); library(ggpubr); library(BiodiversityR)

# setwd
setwd("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/TE Pooling")

#Load some things
source("~/Documents/Grad School/Course Work/Bioinformatics/MergeLowAbund.R")
source("~/Documents/Grad School/Course Work/Bioinformatics/removeNARows.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/g_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/w_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/uw_unifrac.R")
source("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/change16STaxaNames.R")

# import data
tepools <- import_biom("SNPconfirmed_AMR_analytic_matrix.biom")
te_individuals <-import_biom("~/Documents/Research/Projects/USDA_NIFA_AMR/Sequencing/TE/SNPconfirmed_AMR_analytic_matrix_phyloseq.biom")
ARG_data <- merge_phyloseq(tepools,te_individuals)
te_mh <- import_biom("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/TE Pooling/MHEnrichment/MH_count_matrix_normalized.biom")




sample_names(te_mh)

sample_names(te_mh) <- c("H_12-01","H_12-02","H_12-03","H-12-3-SC-142","H-12-3-SC-151","H-12-3-SC-163","H-12-3-SC-171","H-12-3-SC-187","H-12-3-SC-197",
                         "H-12-3-SC-202","H-12-3-SC-209","H-12-3-SC-218","H-12-3-SC-230","H-12-3-SC-233","H-12-3-SC-234","H_12-04","H-12-4-SC-147","H-12-4-SC-148",
                         "H-12-4-SC-162","H-12-4-SC-180","H-12-4-SC-181","H-12-4-SC-183","H-12-4-SC-193","H-12-4-SC-194","H-12-4-SC-207","H-12-4-SC-215","H-12-4-SC-220",
                         "H-12-4-SC-237","H_12-05","H_12-06","HP_03-01","HP_03-02","HP_03-03","HP_03-04","HP_03-05","HP_03-06",
                         "HP_06-01","HP_06-02","HP_06-03","HP_06-04","HP_06-05","HP_06-06","L_12-01","L_12-02","L_12-03",
                         "L_12-04","L_12-05","L_12-06","LP_03-01","LP_03-02","LP_03-03","LP_03-04","LP_03-05","LP_03-06",
                         "LP_06-01","LP_06-02","LP_06-03","LP_06-04","LP_06-05","LP_06-06")

sample_names(map_file)
sample_names(tephyloseq)
write.csv(sample_names(tephyloseq),"samplenamesformetadata_allTE.csv")

resistome_sampledata <- import_qiime_sample_data("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/TEmetadata_resistome.txt")
sample_names(ARG_data)

pooling_resistome.ps <- merge_phyloseq(ARG_data,resistome_sampledata)

tedata # 319 taxa AMR genes (1 is sequins) AMR+Mh genes, 142 samples
sample_names(tedata)
data <- merge_phyloseq(tedata, map_file)

data_mh <-merge_phyloseq(te_mh,map_file)
sample_names(data_mh)

sample_names(data)
sample_names(data_mh)
sample_names(data) <- sample_data(data)$MEG_ID
sample_names(data_mh) <-sample_data(data_mh)$MEG_ID

sample_data(data)$Number.in.Pool <- as.factor(sample_data(data)$Number.in.Pool)
sample_data(data_mh)$Number.in.Pool <- as.factor(sample_data(data_mh)$Number.in.Pool)

# check the names of our ranks
rank_names(data) # "Rank1" - "Rank6" not ideal, lets change em
colnames(tax_table(data)) <- c("Type","Class","Mechanism","Group","Gene","SNP")
rank_names(data) # beauty, now they are named properly

rank_names(data_mh)
colnames(tax_table(data_mh)) <- c("Species","Strain","All")
rank_names(data_mh)


### need to split up the genes needing SNP confirmation:
#data_SNPconfirm <- subset_taxa(data, SNP=="yes", T) # 36 of the 319 genes
#data_SNPconfirm
#data_noSNP <- subset_taxa(data, SNP=="no") # 283 of the 319 genes
#data_noSNP

# some QC checks
min(sample_sums(data)) # 1
max(sample_sums(data)) # 32170
mean(sample_sums(data)) #4984.083
median(sample_sums(data)) #2839.5
sort(sample_sums(data)) # Trim at 1000
write.csv(sort(sample_sums(data)),"datacounts.csv")

min(sample_sums(data_mh)) # 10000 (Normalized by multiplying x 10000, so this is 1)
data_mh <- subset_samples(data_mh, MEG_ID!="SC_193")
any(taxa_sums(data_mh)==0)
sample_data(data_mh)$MEG_ID
write.csv(sort(sample_sums(data_mh)),"data_mhcounts.csv") 

data_250 <- prune_samples(sample_sums(data)>250, data)
data_63 <- subset_samples(data, HP_06.03=="Y")
sum(taxa_sums(data_63)==0)
data_63 <- prune_taxa(taxa_sums(data_63)>0,data_63)
data_63.css <- phyloseq_transform_css(data_63, log = F)
data_63.css.df <- as(sample_data(data_63),"data.frame")
rel_abund_63 <- transform_sample_counts(data_63.css, function(x) {x/sum(x)}*100)

rel_abund_63_group <- tax_glom(rel_abund_63, taxrank = "Group") %>%
  psmelt()

median(sample_data(data)$InputReads)
min(sample_data(data)$InputReads)
median(sample_data(data)$HostRemoved)


data_trimmed_1k <- prune_samples(sample_sums(data)>1000,data)
any(taxa_sums(data_trimmed_1k)==0)
median(sample_sums(data_trimmed_1k))
mean(sample_sums(data_trimmed_1k))
sort(sample_sums(data_trimmed_1k))

data_trimmed_middle <- prune_samples(sample_sums(data_trimmed_1k)<32000, data_trimmed_1k)
any(taxa_sums(data_trimmed_middle)==0)
data_trimmed_middle <- prune_taxa(taxa_sums(data_trimmed_middle)>0, data_trimmed_middle)
median(sample_sums(data_trimmed_middle)) #3757
mean(sample_sums(data_trimmed_middle)) #5844
sort(sample_sums(data_trimmed_middle))
data_trimmed_middle #293 taxa in 45 samples
sample_data(data_trimmed_middle)

data_trimmed_indiv <- subset_samples(data_trimmed_middle, Number.in.Pool==1)
sum(taxa_sums(data_trimmed_indiv)==0)
data_trimmed_indiv <- prune_taxa(taxa_sums(data_trimmed_indiv)>0,data_trimmed_indiv)
ntaxa(data_trimmed_indiv)
data_trimmed_pools <- subset_samples(data_trimmed_middle, Number.in.Pool!=1)
sum(taxa_sums(data_trimmed_pools)==0)
data_trimmed_pools <- prune_taxa(taxa_sums(data_trimmed_pools)>0, data_trimmed_pools)
ntaxa(data_trimmed_pools)

sample_data(data_trimmed_1k)$Number.in.Pool <- as.factor(sample_data(data_trimmed_1k)$Number.in.Pool)

write.csv(sample_sums(data_trimmed_1k),"datacounts_trimmed.csv")
write.csv(sample_sums(data_trimmed_middle), "datacounts_middle.csv")

# any taxa with no counted reads
sum(taxa_sums(data_trimmed_1k)==0) # 64
data_trimmed_1k <- prune_taxa(taxa_sums(data_trimmed_1k) > 0, data_trimmed_1k)
any(sample_sums(data_trimmed_1k)==0) # nope


##### Sequencing depth (based on reads passing QC that were put into ARM++)
data.df <- as(sample_data(data_trimmed_1k), "data.frame")
data_middle.df <-as(sample_data(data_trimmed_middle),"data.frame")

View(sample_data(data_trimmed_1k))
##Number in pool palette
#when individual swabs add this back: "#B3DD79", 

numberinpoolpalette <- c("#30123BFF","#1AE4B6FF" ,"#FABA39FF" ,"#7A0403FF")
poolonly_palette<- c("#1AE4B6FF" ,"#FABA39FF" ,"#7A0403FF")

poolcomparisons <- list(c("1","3"),c("1","6"),c("1","12"),c("3","6"),c("3","12"),c("6","12"))
poolonly_comparisons <- list(c("3","6"),c("3","12"),c("6","12"))
# just compare the "pools"
ggplot(data.df, aes(x= Number.in.Pool, y= InputReads, fill = Number.in.Pool, alpha=Number.in.Pool,color=Number.in.Pool)) + theme_bw() +
  geom_errorbar(stat = "summary", color="black")+
  geom_bar(stat="summary") +
  geom_point()+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_color_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))
  #stat_compare_means(label.y = 80000000)#+ #Differences, KW: p=0.02
 # stat_compare_means(comparisons = poolcomparisons, method = "wilcox", p.adjust.methods="bonferoni")

sample_sum_df_ARG <- data.frame(ARGs = sample_sums(data_trimmed_middle))
sample_sum_meta_ARG <- as(sample_data(data_trimmed_middle), "data.frame") # making into DF for metadata
sample_sum_stats_ARG <- cbind(sample_sum_df_ARG, sample_sum_meta_ARG)
kruskal_test(sample_sum_stats_ARG, ARGs~Number.in.Pool)

sample_sum_stats_ARG$Number.in.Pool <- as.factor(sample_sum_stats_ARG$Number.in.Pool)

sample_sum_stats_ARG %>%
  group_by(Number.in.Pool)%>%
  summarise(min=min(ARGs), med=median(ARGs), max=max(ARGs))



ggplot(data_middle.df, aes(x= Number.in.Pool, y= InputReads, fill = Number.in.Pool, alpha=Number.in.Pool,color=Number.in.Pool)) + theme_bw() +
  geom_errorbar(stat = "summary", color="black")+
  geom_bar(stat="summary") +
  scale_fill_manual(values=numberinpoolpalette)+
  scale_color_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))

kruskal_test(data.df,InputReads~Number.in.Pool) #P=0.0199
dunn_test(data.df, InputReads~Number.in.Pool, p.adjust.method = "BH")
# 1 v 12 different (p=0.01)

kruskal_test(data_middle.df, InputReads~Number.in.Pool) #0.0234
dunn_test(data_middle.df, InputReads~Number.in.Pool, p.adjust.method = "BH")
#1 v 12 (P=0.0172)

#Contains Mh culture positive or notor not
kruskal_test(data.df, InputReads~Contains.CP)




# differences

#Individuals
indi_mh <- subset_samples(data_mh, Number.in.Pool==1)
tax_table(indi_mh)

indi_mh <- transform_sample_counts(indi_mh, function(x) round({x}*10, digits = 0))
high_mh_alpha <- transform_sample_counts(high_mh, function(x) round({x}*10, digits = 0))
otu_table(high_mh_alpha)
low_mh_alpha <- transform_sample_counts(low_mh, function(x) round({x}*10, digits = 0))
otu_table(low_mh_alpha)

data_alpha_mh <- transform_sample_counts(data_mh, function(x) round({x}*10, digits = 0))
otu_table(data_alpha_mh)

##HP
high <- subset_samples(data_trimmed_1k, PrevH=="Y") #29 samples, 251 taxa
sum(taxa_sums(high)==0)
high <- prune_taxa(taxa_sums(high)>0,high) 
any(sample_sums(high)==0)
View(sample_data(high))


high_250 <- subset_samples(data_250, PrevH=="Y")
sum(taxa_sums(high_250)==0)
high_250 <- prune_taxa(taxa_sums(high_250)>0,high_250)

high_trimmed <- subset_samples(data_trimmed_middle, PrevH=="Y") #28 samples, 293 taxa
sum(taxa_sums(high_trimmed)==0)
high_trimmed <- prune_taxa(taxa_sums(high_trimmed)>0,high_trimmed) 
any(sample_sums(high_trimmed)==0)
View(sample_data(high_trimmed))

hpools <-subset_samples(high, Number.in.Pool!="1") #16 samples, 234 taxa
sum(taxa_sums(hpools)==0)
hpools <- prune_taxa(taxa_sums(hpools)>0,hpools)
any(sample_sums(hpools)==0)

hpools_trimmed <-subset_samples(high_trimmed, Number.in.Pool!="1") #16 samples, 234 taxa
sum(taxa_sums(hpools_trimmed)==0)
hpools_trimmed <- prune_taxa(taxa_sums(hpools_trimmed)>0,hpools)
any(sample_sums(hpools_trimmed)==0)

##LP
low <- subset_samples(data_trimmed_1k,PrevL=="Y") #24 samples, 270 taxa
sum(taxa_sums(low)==0)
low <- prune_taxa(taxa_sums(low)>0,low)
any(sample_sums(low)==0)

low_trimmed <- subset_samples(data_trimmed_middle, PrevL=="Y")
low_trimmed <- subset_samples(low_trimmed, Number.in.Pool!=1)
sum(taxa_sums(low_trimmed)==0)
low_trimmed <- prune_taxa(taxa_sums(low_trimmed)>0,low_trimmed) #24 samples, 270 taxa
any(sample_sums(low_trimmed)==0)
View(sample_data(low_trimmed))

lpools <-subset_samples(low, Number.in.Pool!="1") #17 samples, 265 taxa
sum(taxa_sums(lpools)==0)
lpools <- prune_taxa(taxa_sums(lpools)>0,lpools)
any(sample_sums(lpools)==0)

#pools
pools_amr <- subset_samples(data_trimmed_middle, SampleType=="Pool")
sum(taxa_sums(pools_amr)==0) #7 args found only inidividuals
pools_amr <- prune_taxa(taxa_sums(pools_amr)>0, pools_amr)

#Checking input reads
high.df <- as(sample_data(high), "data.frame")
kruskal_test(high.df,InputReads~Number.in.Pool) #NS

high_trimmed.df <- as(sample_data(high_trimmed), "data.frame")
kruskal_test(high_trimmed.df,InputReads~Number.in.Pool) #NS, P=0.104
high_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(InputReads), median=median(InputReads), n=n(),sd=sd(InputReads))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
high_trimmed.df %>%
  reframe(mean=mean(InputReads), median=median(InputReads), n=n(),sd=sd(InputReads))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)


kruskal_test(high_trimmed.df,HostRemoved~Number.in.Pool)
dunn_test(high_trimmed.df, HostRemoved~Number.in.Pool, p.adjust.method = "BH")
high_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(HostRemoved), median=median(HostRemoved), n=n(),sd=sd(HostRemoved))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
high_trimmed.df %>%
  reframe(mean=mean(HostRemoved), median=median(HostRemoved), n=n(),sd=sd(HostRemoved))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(high_trimmed.df,MH_align~Number.in.Pool)
dunn_test(high_trimmed.df, MH_align~Number.in.Pool, p.adjust.method = "BH")
high_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(MH_align), median=median(MH_align), n=n(),sd=sd(MH_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

high_trimmed.df %>%
  reframe(mean=mean(MH_align), median=median(MH_align), n=n(),sd=sd(MH_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(high_trimmed.df,MHPerc~Number.in.Pool)
kruskal_test(high_trimmed.df,StrainPer~Number.in.Pool)


kruskal_test(high_trimmed.df,Strain_align~Number.in.Pool)
dunn_test(high_trimmed.df, Strain_align~Number.in.Pool, p.adjust.method = "BH")
high_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Strain_align), median=median(Strain_align), n=n(),sd=sd(Strain_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

high_trimmed.df %>%
  reframe(mean=mean(Strain_align), median=median(Strain_align), n=n(),sd=sd(Strain_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(high_trimmed.df,MHPerc~Number.in.Pool)
kruskal_test(high_trimmed.df,StrainPer~Number.in.Pool)

View(high_trimmed.df)

low.df <- as(sample_data(low),"data.frame")
kruskal_test(low.df, InputReads~Number.in.Pool) #NS, P=0.12

low_trimmed.df <- as(sample_data(low_trimmed),"data.frame")
kruskal_test(low_trimmed.df, InputReads~Number.in.Pool) #NS, P=0.141

low_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(InputReads), median=median(InputReads), n=n(),sd=sd(InputReads))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

low_trimmed.df %>%
  reframe(mean=mean(InputReads), median=median(InputReads), n=n(),sd=sd(InputReads))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)


kruskal_test(low_trimmed.df, HostRemoved~Number.in.Pool) #NS, P=0.141

low_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(HostRemoved), median=median(HostRemoved), n=n(),sd=sd(HostRemoved))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

low_trimmed.df %>%
  reframe(mean=mean(HostRemoved), median=median(HostRemoved), n=n(),sd=sd(HostRemoved))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(low_trimmed.df, MH_align~Number.in.Pool) #NS, P=0.141


low_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(MH_align), median=median(MH_align), n=n(),sd=sd(MH_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

low_trimmed.df %>%
  reframe(mean=mean(MH_align), median=median(MH_align), n=n(),sd=sd(MH_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(low_trimmed.df, Strain_align~Number.in.Pool) #NS, P=0.141

low_trimmed.df %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Strain_align), median=median(Strain_align), n=n(),sd=sd(Strain_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

low_trimmed.df %>%
  reframe(mean=mean(Strain_align), median=median(Strain_align), n=n(),sd=sd(Strain_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

pools_amr.df <- as(sample_data(pools_amr),"data.frame")
wilcox_test(pools_amr.df, InputReads~Prevalence) #P=0.502


kruskal_test(low_trimmed.df,MHPerc~Number.in.Pool)
kruskal_test(low_trimmed.df,StrainPer~Number.in.Pool)

#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################
#alpha_div <- estimate_richness(data_trimmed_1k, measures = c("Observed","Shannon","Simpson","InvSimpson"))
#alpha_div.df <- as(sample_data(data_trimmed_1k), "data.frame")
#alpha_div_meta <- cbind(alpha_div, alpha_div.df)
#alpha_div_meta


high_alpha_div <- estimate_richness(high, measures = c("Observed","Shannon","Simpson","InvSimpson"))
high_alpha_div.df <- as(sample_data(high), "data.frame")
high_alpha_div_meta <- cbind(high_alpha_div, high_alpha_div.df)
high_alpha_div_meta

mh_alpha_div <- estimate_richness(data_alpha_mh, measures = c("Observed","Shannon","Simpson","InvSimpson"))
mh_alpha_div.df <- as(sample_data(data_alpha_mh), "data.frame")
mh_alpha_div_meta <- cbind(mh_alpha_div, mh_alpha_div.df)

mh_alpha_div_meta %>%
  reframe(median=median(HostRemoved),mean=mean(HostRemoved), sd=sd(HostRemoved), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

mh_alpha_div_meta %>%
  reframe(median=median(MH_align),mean=mean(MH_align), sd=sd(MH_align), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

mh_alpha_div_meta %>%
  reframe(median=median(Strain_align),mean=mean(Strain_align), sd=sd(Strain_align), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

mh_alpha_div_meta %>%
  reframe(median=median(MHPerc),mean=mean(MHPerc), sd=sd(MHPerc), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

mh_alpha_div_meta %>%
  reframe(median=median(MHPerc),mean=mean(MHPerc), sd=sd(MHPerc), n=n(), min=min(MHPerc), max=max(MHPerc))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

mh_alpha_div_meta %>%
  reframe(median=median(StrainPer),mean=mean(StrainPer), sd=sd(StrainPer), n=n(), min=min(StrainPer), max=max(StrainPer))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

mh_high_metadata <- subset(mh_alpha_div_meta, PrevH=="Y")
mh_low_metadata <-subset(mh_alpha_div_meta, PrevH=="N")

mh_high_metadata%>%
  reframe(median=median(MHPerc),mean=mean(MHPerc), sd=sd(MHPerc), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
mh_high_metadata%>%
  reframe(median=median(StrainPer),mean=mean(StrainPer), sd=sd(StrainPer), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
kruskal_test(mh_high_metadata, MHPerc~Number.in.Pool)
kruskal_test(mh_high_metadata, StrainPer~Number.in.Pool)

kruskal_test(mh_low_metadata, MHPerc~Number.in.Pool)
kruskal_test(mh_low_metadata, StrainPer~Number.in.Pool)


high_mh_per_plot <- ggplot(mh_high_metadata, aes(x= Number.in.Pool, y = MHPerc, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "low_trimmed", y= "Richness") +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  #geom_text(data = ls_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits=c(0,100))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y=element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

high_strain_per_plot <- ggplot(mh_high_metadata, aes(x= Number.in.Pool, y = StrainPer, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "low_trimmed", y= "Richness") +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  #geom_text(data = ls_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits=c(0,1.5))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y=element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

high_mh_strainplot<-high_mh_per_plot + 
  annotation_custom(
    ggplotGrob(high_strain_per_plot), 
    xmin = 2.5, xmax = 4.5, ymin = 5, ymax = 40)

low_mh_per_plot <- ggplot(mh_low_metadata, aes(x= Number.in.Pool, y = MHPerc, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "low_trimmed", y= "Richness") +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  #geom_text(data = ls_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits=c(0,100))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

low_strain_per_plot <- ggplot(mh_low_metadata, aes(x= Number.in.Pool, y = StrainPer, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "low_trimmed", y= "Richness") +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  #geom_text(data = ls_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits=c(0,1.5))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y=element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

low_mh_strainplot<-low_mh_per_plot + 
  annotation_custom(
    ggplotGrob(low_strain_per_plot), 
    xmin = 2, xmax = 3.5, ymin = 5, ymax = 40)

ggarrange(high_mh_strainplot, low_mh_strainplot, nrow = 1)

mh_alpha_div_meta %>%
  group_by(Number.in.Pool) %>%
  wilcox_test(MH_align~MHqPCR)

mh_alpha_div_meta$Number.in.Pool <- factor(mh_alpha_div_meta$Number.in.Pool, levels=c(1,3,6,12))

ggplot(mh_alpha_div_meta, aes(x=Number.in.Pool, y=MH_align, color=MHqPCR))+
  geom_boxplot()+
  scale_color_manual(values=c("seagreen","darkblue"),limits=c("YES","NO"))


high_mh_alpha_div <- estimate_richness(high_mh_alpha, measures = c("Observed","Shannon","Simpson","InvSimpson"))
high_mh_alpha_div.df <- as(sample_data(high_mh_alpha), "data.frame")
high_mh_alpha_div_meta <- cbind(high_mh_alpha_div, high_mh_alpha_div.df)
low_mh_alpha_div <- estimate_richness(low_mh_alpha, measures="Observed")
low_mh_alpha_div.df <- as(sample_data(low_mh_alpha), "data.frame")
low_mh_alpha_div_meta <- cbind(low_mh_alpha_div, low_mh_alpha_div.df)

all_trimmed_ARG_div <- estimate_richness(data_trimmed_middle, measures = c("Observed","Shannon"))
all_trimmed_ARG_div.df <-as(sample_data(data_trimmed_middle),"data.frame")
all_trimmed_ARG_meta <- cbind(all_trimmed_ARG_div,all_trimmed_ARG_div.df)
all_trimmed_ARG_meta

all_trimmed_indiv <- data_trimmed_middle

high_trimmed_alpha_div <- estimate_richness(high_trimmed, measures = c("Observed","Shannon","Simpson","InvSimpson"))
high_trimmed_alpha_div.df <- as(sample_data(high_trimmed), "data.frame")
high_trimmed_alpha_div_meta <- cbind(high_trimmed_alpha_div, high_trimmed_alpha_div.df)
high_trimmed_alpha_div_meta

low_alpha_div <- estimate_richness(low, measures = c("Observed","Shannon","Simpson","InvSimpson"))
low_alpha_div.df <- as(sample_data(low), "data.frame")
low_alpha_div_meta <- cbind(low_alpha_div, low_alpha_div.df)
low_alpha_div_meta

low_trimmed_alpha_div <- estimate_richness(low_trimmed, measures = c("Observed","Shannon","Simpson","InvSimpson"))
low_trimmed_alpha_div.df <- as(sample_data(low_trimmed), "data.frame")
low_trimmed_alpha_div_meta <- cbind(low_trimmed_alpha_div, low_trimmed_alpha_div.df)
low_trimmed_alpha_div_meta

indi_mh_alpha_div <- estimate_richness(indi_mh, measures = c("Observed","Shannon","Simpson","InvSimpson"))
indi_mh_alpha_div.df <-as(sample_data(indi_mh),"data.frame")
indi_mh_alpha_div_meta <- cbind(indi_mh_alpha_div, indi_mh_alpha_div.df)

pools_amr_alpha_div <- estimate_richness(pools_amr,measures = c("Observed","Shannon","Simpson","InvSimpson"))
pools_amr_alpha_div.df <- as(sample_data(pools_amr), "data.frame")
pools_amr_alpha_div_meta <- cbind(pools_amr_alpha_div, pools_amr_alpha_div.df)


#Poolsize-High
high_alpha_div_meta %>%
  wilcox_test(Observed~Contains.CP) #P=0.52
high_trimmed_alpha_div_meta %>%
  kruskal_test(Observed~Number.in.Pool) #P=0.52


ggplot(high_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "High", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

kruskal_test(high_trimmed_alpha_div_meta, Observed~Number.in.Pool) #P=0.238, NS
min(high_trimmed_alpha_div_meta$Observed)
median(high_trimmed_alpha_div_meta$Observed)

hs_observed <- ggplot(high_trimmed_alpha_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Richness") +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

kruskal_test(high_trimmed_alpha_div_meta, Shannon~Number.in.Pool) #P=0.259, NS
dunn_test(high_trimmed_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH")

hs_shannon <- ggplot(high_trimmed_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

#Poolsize-Low
kruskal_test(low_alpha_div_meta, Observed~Number.in.Pool) #P=0.063, NS
kruskal_test(low_trimmed_alpha_div_meta, Observed~Number.in.Pool) #P=0.595, NS
dunn_test(low_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH")

ls_observed <- ggplot(low_alpha_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "No. unique ARGs") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        title = element_blank(),       
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ls_shannon <- ggplot(low_trimmed_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "No. unique ARGs") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggarrange(hs_observed, ls_observed, hs_shannon, ls_shannon, nrow=2, ncol=2)

kruskal_test(low_alpha_div_meta, Shannon~Number.in.Pool) #P=0.080, NS
dunn_test(low_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH")

kruskal_test(low_trimmed_alpha_div_meta, Shannon~Number.in.Pool) #P=0.863



ggplot(low_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggplot(low_trimmed_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

##Compare differences faceted by poolsize
ggplot(alpha_div_meta, aes(x= Prevalence, y = Observed, fill =  Prevalence, colour =  Prevalence)) + 
  theme_bw() + 
  labs(title= "", y= "No. unique ARGs") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  #scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  facet_wrap(~Number.in.Pool,ncol=1)+
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())+
  stat_compare_means() #+
#stat_compare_means(comparisons = poolcomparisons)

ggplot(alpha_div_meta, aes(x= Prevalence, y = Shannon, fill =  Prevalence, colour =  Prevalence)) + 
  theme_bw() + 
  labs(title= "", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  #scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  facet_wrap(~Number.in.Pool,ncol=1)+
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())+
  stat_compare_means()


#Prevalence
wilcox_test(pools_amr_alpha_div_meta, Observed~Prevalence) #P=0.00967
wilcox_test(pools_amr_alpha_div_meta, Shannon~Prevalence) #P=0.00997

ggplot(pools_amr_alpha_div_meta, aes(x= Prevalence, y = Observed, fill =  Prevalence, colour =  Prevalence)) + 
  theme_bw() + 
  labs(title= "", y= "Observed ARGs") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  #scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  #facet_wrap(~Number.in.Pool,ncol=1)+
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggplot(pools_amr_alpha_div_meta, aes(x= Prevalence, y = Shannon, fill =  Prevalence, colour =  Prevalence)) + 
  theme_bw() + 
  labs(title= "", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  #scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  #facet_wrap(~Number.in.Pool,ncol=1)+
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =14, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

#############################################################################################
##############################         BETA DIVERSITY         ###############################
#############################################################################################
#############################################################################################

#high 
sum(taxa_sums(high_trimmed)==0)
sum(sample_sums(high_trimmed)==0)
high.css <- phyloseq_transform_css(high, log = F)
high_trimmed.css <-phyloseq_transform_css(high_trimmed, log = F)
high.css.df <- as(sample_data(high.css), "data.frame")
high_trimmed.css.df <- as(sample_data(high_trimmed.css),"data.frame")


#low
sum(taxa_sums(low)==0)
low.css <- phyloseq_transform_css(low, log = F)
low.css.df <-as(sample_data(low.css),"data.frame")

sum(taxa_sums(low_trimmed)==0)
low_trimmed.css <- phyloseq_transform_css(low_trimmed, log = F)
low_trimmed.css.df <-as(sample_data(low_trimmed.css),"data.frame")

high.dist <- vegdist(t(otu_table(high.css)), method = "bray")
high_trimmed.dist <- vegdist(t(otu_table(high_trimmed.css)),method = "bray")
low.dist <- vegdist(t(otu_table(low.css)),method = "bray")
low_trimmed.dist <- vegdist(t(otu_table(low_trimmed.css)),method = "bray")

high.ord <- vegan::metaMDS(comm = t(otu_table(high.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)
high_trimmed.ord <- vegan::metaMDS(comm = t(otu_table(high_trimmed.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)

low.ord <- vegan::metaMDS(comm = t(otu_table(low.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)
low_trimmed.ord <- vegan::metaMDS(comm = t(otu_table(low_trimmed.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)


###Pools
pools_amr.css <- phyloseq_transform_css(pools_amr,log=F)
pools_amr.css.df <- as(sample_data(pools_amr.css),"data.frame")

pools_amr.dist <- vegdist(t(otu_table(pools_amr.css)), method = "bray")
pools_amr.ord <- vegan::metaMDS(comm = t(otu_table(pools_amr.css)),distance = "bray",try = 20, trymax = 50, autotransform = F)





plot_ordination(high.css,high.ord,  color = "Number.in.Pool") +
  scale_color_manual(values = numberinpoolpalette)+
  stat_ellipse()

plot_ordination(high_trimmed.css,high_trimmed.ord,  color = "Number.in.Pool") +
  scale_color_manual(values = numberinpoolpalette)+
  stat_ellipse()

plot_ordination(low.css, low.ord, color="Number.in.Pool")+
  scale_color_manual(values = numberinpoolpalette)+
  stat_ellipse()

plot_ordination(low_trimmed.css, low_trimmed.ord, color="Number.in.Pool")+
  scale_color_manual(values = poolonly_palette)+
  stat_ellipse()

plot_ordination(pools_amr.css,pools_amr.ord,  color = "Prevalence") +
  #scale_color_manual(values = numberinpoolpalette)+
  stat_ellipse()

##Adding centroids
#### findingcenters
#hs
high_plot1 <- ordiplot(high.ord$points)
high_siteslong <-sites.long(high_plot1,high.css.df)

high_trimmed_plot1 <- ordiplot(high_trimmed.ord$points)
high_trimmed_siteslong <-sites.long(high_trimmed_plot1,high_trimmed.css.df)

high_centroids <- envfit(high.ord~high.css.df$Number.in.Pool)
high_centroids

high_col <- c("1","3","6","12")
#poolonly_col <- c("3","6","12")
high_NMDS_col1 <- c(0.1023,0.1896,-0.2341,-0.2593)
high_NMDS_col2 <- c(0.0731,-0.0275,-0.0704,-0.0867)
high_centroids.df <-data.frame(high_col,high_NMDS_col1,high_NMDS_col2)

high_trimmed_centroids <- envfit(high_trimmed.ord~high_trimmed.css.df$Number.in.Pool)
high_trimmed_centroids

high_trimmed_col <- c("1","3","6","12")
high_trimmed_NMDS_col1 <- c(0.0901,0.1693,-0.1913,-0.2845)
high_trimmed_NMDS_col2 <- c(0.0616,-0.0650,0.0580,-0.1287)
high_trimmed_centroids.df <-data.frame(high_trimmed_col,high_trimmed_NMDS_col1,high_trimmed_NMDS_col2)

#ls
low_plot1 <- ordiplot(low.ord$points)
low_siteslong <- sites.long(low_plot1,low.css.df)
low_centroids <- envfit(low.ord~low.css.df$Number.in.Pool)
low_centroids

low_col <- c("1","3","6","12")
low_NMDS_col1 <- c(-0.4135,0.0162,0.2731,0.1957)
low_NMDS_col2 <- c(-0.1129,0.0262,0.2047,-0.0949)
low_centroids.df <-data.frame(low_col,low_NMDS_col1,low_NMDS_col2)

low_trimmed_plot1 <- ordiplot(low_trimmed.ord$points)
low_trimmed_siteslong <- sites.long(low_trimmed_plot1,low_trimmed.css.df)
low_trimmed_centroids <- envfit(low_trimmed.ord~low_trimmed.css.df$Number.in.Pool)
low_trimmed_centroids

low_trimmed_col <- c("3","6","12")
low_trimmed_NMDS_col1 <- c(0.0344,-0.1344,0.1057)
low_trimmed_NMDS_col2 <- c(-0.1792,0.1114,0.0380)
low_trimmed_centroids.df <-data.frame(low_trimmed_col,low_trimmed_NMDS_col1,low_trimmed_NMDS_col2)

#Pools only
pools_plot1 <- ordiplot(pools_amr.ord$points)
pools_siteslong <- sites.long(pools_plot1, pools_amr.css.df)
pools_centroids <- envfit(pools_amr.ord~pools_amr.css.df$Prevalence)

prevalence_col <- c("H","L")
prev_NMDS_col1 <- c(0.0949,-0.0837)
prev_NMDS_col2 <- c(-0.1021, 0.0901)
prev_centroids.df <- data.frame(prevalence_col, prev_NMDS_col1, prev_NMDS_col2)

## high
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, High Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = high_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool, fill=Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(18,19,19,19))+
  scale_size_manual(values = c(5,5,5,5), guide = "none") +
  stat_ellipse(data = high_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = high_centroids.df, aes(x=high_NMDS_col1, y=high_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(18,19,19,19)) +
  geom_text(data = high_centroids.df, aes(x=high_NMDS_col1, y=high_NMDS_col2, label = high_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, colour = "black"))


# stats
high.adonis <- adonis2(high.dist~Number.in.Pool, high.css.df, perm = 9999)
high.adonis #NS, 0.5
pairwise.adonis2(high.dist~Number.in.Pool, high.css.df, perm=9999, p.adjust.methods="BH")
write.csv(high.adonis, "TE_highprev_adonis.csv")

high.disper <- betadisper(high.dist, high_alpha_div.df$Number.in.Pool)
plot(high.disper)
high.permdisp <- permutest(high.disper, permutations = 9999, pairwise = F) #No dif, 0.1052
high.permdisp_pair <- permutest(high.disper, permutations = 9999, pairwise = T) #pairwise 12 diff from 1 and 3
write.csv(high.permdisp[["pairwise"]][["permuted"]],"highprev_TE_bray_permdisp.csv")

## high_trimmed
high_trimmed_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, High Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = high_trimmed_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool, fill=Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(18,19,19,19))+
  scale_size_manual(values = c(2,2,2,2), guide = "none") +
  stat_ellipse(data = high_trimmed_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = high_trimmed_centroids.df, aes(x=high_trimmed_NMDS_col1, y=high_trimmed_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = high_trimmed_centroids.df, aes(x=high_trimmed_NMDS_col1, y=high_trimmed_NMDS_col2, label = high_trimmed_col), colour = "white", size = 2, fontface = "bold") +
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

ggsave("hs_arg_nmds.tiff", plot = high_trimmed_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


# stats
high_trimmed.adonis <- adonis2(high_trimmed.dist~Number.in.Pool, high_trimmed.css.df, perm = 9999)
high_trimmed.adonis #NS, 0.52
pairwise.adonis2(high_trimmed.dist~Number.in.Pool, high_trimmed.css.df, perm=9999, p.adjust.methods="BH")
write.csv(high_trimmed.adonis, "TE_high_trimmedprev_adonis.csv")

high_trimmed.disper <- betadisper(high_trimmed.dist, high_trimmed_alpha_div.df$Number.in.Pool)
plot(high_trimmed.disper)
high_trimmed.permdisp <- permutest(high_trimmed.disper, permutations = 9999, pairwise = F) #No dif, 0.1139
high_trimmed.permdisp_pair <- permutest(high_trimmed.disper, permutations = 9999, pairwise = T) #pairwise 12 diff from 1 and 3
write.csv(high_trimmed.permdisp[["pairwise"]][["permuted"]],"high_trimmedprev_TE_bray_permdisp.csv")

## low
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, low Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = low_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_size_manual(values = c(5,5,5,5), guide = "none") +
  stat_ellipse(data = low_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = low_centroids.df, aes(x=low_NMDS_col1, y=low_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(18,19,19,19)) +
  geom_text(data = low_centroids.df, aes(x=low_NMDS_col1, y=low_NMDS_col2, label = low_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", size = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, colour = "black"))

## low_trimmed
low_trimmed_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, low_trimmed Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous() +
  scale_y_continuous(breaks = c(-1,0,1)) +
  scale_shape_manual(values = c(19,19,19)) +
  geom_point(data = low_trimmed_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_size_manual(values = c(2,2,2), guide = "none") +
  stat_ellipse(data = low_trimmed_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = low_trimmed_centroids.df, aes(x=low_trimmed_NMDS_col1, y=low_trimmed_NMDS_col2), fill = poolonly_palette, colour = poolonly_palette, size = 3, shape = c(19,19,19)) +
  geom_text(data = low_trimmed_centroids.df, aes(x=low_trimmed_NMDS_col1, y=low_trimmed_NMDS_col2, label = low_trimmed_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = poolonly_palette) +
  scale_fill_manual(values = poolonly_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"),
        plot.title=element_blank())

ggsave("ls_arg_nmds.tiff", plot = low_trimmed_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")

# stats
low.adonis <- adonis2(low.dist~Number.in.Pool, low.css.df,perm = 9999) #NS
low.adonis #NS
write.csv(low.adonis, "TE_lowprev_adonis.csv")

low.disper <- betadisper(low.dist, low_alpha_div.df$Number.in.Pool)
plot(low.disper)
low.permdisp <- permutest(low.disper, permutations = 9999, pairwise = F) #No dif, 0.3251
low.permdisp_pair <- permutest(low.disper, permutations = 9999, pairwise = T) #
write.csv(low.permdisp[["pairwise"]][["permuted"]],"lowprev_TE_bray_permdisp.csv")

# stats
low_trimmed.adonis <- adonis2(low_trimmed.dist~Number.in.Pool, low_trimmed.css.df,perm = 9999) #NS
low_trimmed.adonis #NS,0.8
write.csv(low_trimmed.adonis, "TE_low_trimmedprev_adonis.csv")

low_trimmed.disper <- betadisper(low_trimmed.dist, low_trimmed_alpha_div.df$Number.in.Pool)
plot(low_trimmed.disper)
low_trimmed.permdisp <- permutest(low_trimmed.disper, permutations = 9999, pairwise = F) #No dif, 0.103
low_trimmed.permdisp_pair <- permutest(low_trimmed.disper, permutations = 9999, pairwise = T) #6 dif from 12
write.csv(low_trimmed.permdisp[["pairwise"]][["permuted"]],"low_trimmedprev_TE_bray_permdisp.csv")

####
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = pools_siteslong, aes(x=axis1,y=axis2, colour= Prevalence, shape = Prevalence, alpha = Prevalence, size = Prevalence, fill=Prevalence)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(19,19))+
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Prevalence, fill = Prevalence), alpha = c(0.1), lty=1, level = 0.95, linewidth = 2) +
  geom_point(data = prev_centroids.df, aes(x=prev_NMDS_col1, y=prev_NMDS_col2), fill =c("indianred","skyblue"), colour = c("indianred","skyblue"), size = 10, shape = c(19,19)) +
  geom_text(data = prev_centroids.df, aes(x=prev_NMDS_col1, y=prev_NMDS_col2, label = prevalence_col), colour = "black", size = 6, fontface = "bold") +
  #scale_colour_manual(values = numberinpoolpalette) +
  #scale_fill_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", size = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, colour = "black"))


# stats
prev.adonis <- adonis2(pools_amr.dist~Prevalence, pools_amr.css.df, perm = 9999)
prev.adonis #NS, 0.0685   


prev.disper <- betadisper(pools_amr.dist, pools_amr.css.df$Prevalence)
plot(prev.disper)
prev.permdisp <- permutest(prev.disper, permutations = 9999, pairwise = F) #No dif, 0.4737
prev.permdisp_pair <- permutest(prev.disper, permutations = 9999, pairwise = T) #pairwise 12 diff from 1 and 3
write.csv(prev.permdisp[["pairwise"]][["permuted"]],"prevprev_TE_bray_permdisp.csv")




#### CLUSTERING ON BRAY-CURTIS
sample_data(ARG_high.css)$Number.in.Pool <- as.factor(sample_data(ARG_high.css)$Number.in.Pool)
ARG_high.dist <- vegdist(t(otu_table(ARG_high.css)), method = "bray")
ARG_high.hclust <- hclust(ARG_high.dist, method = "ward.D2")
plot(ARG_high.hclust)
ARG_high.dendro <- as.dendrogram(ARG_high.hclust)
ARG_high.dendro.data <- dendro_data(ARG_high.dendro, type = "rectangle")
ARG_high.dendro.data

ARG_high_dendro_metadata <- as_tibble(ARG_high.css@sam_data)
ARG_high.dendro.data$labels <- ARG_high.dendro.data$labels %>%
  left_join(ARG_high_dendro_metadata , by = c("label" = "samplename"))

ARG_high_dendro_plot<-ggplot(ARG_high.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = ARG_high.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 5, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  #geom_text(data = ARG_high.dendro.data$labels, 
   #         aes(x=x,y=y, label=HSPool),size=3,position = position_nudge(y=-0.3,x=0),color="black",angle=90)+
  #scale_y_continuous(limits = c(-.25,1)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_blank())


ARG_high_order <- ARG_high.dendro.data$labels$label

low.hclust <- hclust(low.dist, method = "ward.D2")
plot(low.hclust)
low.dendro <- as.dendrogram(low.hclust)
low.dendro.data <- dendro_data(low.dendro, type = "rectangle")
low.dendro.data

low_trimmed.hclust <- hclust(low_trimmed.dist, method = "ward.D2")
plot(low_trimmed.hclust)
low_trimmed.dendro <- as.dendrogram(low_trimmed.hclust)
low_trimmed.dendro.data <- dendro_data(low_trimmed.dendro, type = "rectangle")
low_trimmed.dendro.data

low_trimmed_dendro_metadata <- as_tibble(low_trimmed.css@sam_data)
low_trimmed.dendro.data$labels <- low_trimmed.dendro.data$labels %>%
  left_join(low_trimmed_dendro_metadata, by = c("label" = "MEG_ID"))

low_amr_dendro_plot<-ggplot(low_trimmed.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = low_trimmed.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 5, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  #geom_text(data = low_trimmed.dendro.data$labels, 
          #  aes(x=x,y=y, label=HSPool),size=3,position = position_nudge(y=-0.3,x=0),color="black",angle=90)+
  #scale_y_continuous(limits = c(-.25,1)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = poolonly_palette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_blank())


low_amr_order <- low_trimmed.dendro.data$labels$label


#############################################################################################
############################         RELATIVE ABUNDANCE         #############################
#############################################################################################
#############################################################################################


##Doing this for legend
any(sample_sums(data_trimmed_middle)==0)
data.css <- phyloseq_transform_css(data_trimmed_middle, log=F)
ra_all <- transform_sample_counts(data.css, function(x) {x/sum(x)}*100)

ra_Class <- tax_glom(ra_all, taxrank = "Class")
ra_Class_melt <- psmelt(ra_Class)

ra_group <- tax_glom(ra_all, taxrank = "Group")
ra_group_melt <- psmelt(ra_group)

unique(ra_group_melt$Group)
length(unique(ra_Class_melt$Class)) #18


abundance_class_df <- ra_Class_melt %>%
  group_by(Class)%>%
  reframe(class=Class,meanRA=mean(Abundance))

sort(abundance_class_df$meanRA)

abundance_Class_df <- ra_Class_melt %>%
  group_by(Class)%>%
  reframe(meanRA=mean(Abundance))



sort(unique(ra_Class_melt$Class))

View(abundance_Class_df)

class_palette <- distinctColorPalette(k=18)

#[1] "#D0E9BF" "#E2E547" "#80B2DB" "#C79497" "#73E3DC" "#E0A84D" "#DBBB8E" "#DF69DA" "#9147E0" "#75E39E" "#8377C9" "#E0615B" "#D5B4E3" "#7EE557" "#7B9F83" "#D9679F" "#C0D376" "#D2DADF"


all_group_palette <- c("#cb3f54",
                       "#df6f70",
                       "#8c282f",
                       "#d94d4d",
                       "#90281f",
                       "#d0735a",
                       "#b84828",
                       "#e08053",
                       "#cf6c2f",
                       "#7d3f13",
                       "#d38725",
                       "#ca9457",
                       "#e16fbc",
                       "#b53c8f",
                       "#bb689c",
                       "#ace268",
                       "#468027",
                       "#a1dd8a",
                       "#65b652",
                       "#296021",
                       "#569c56",
                       "#5bd581",
                       "#399f69",
                       "#6fdfa2",
                       "#48f2a7",
                       "#37a886",
                       "#55e5c1",
                       "#3be6ea",
                       "#5999e0",
                       "#4663ab",
                       "#a27521",
                       "#dcab56",
                       "#ceaa34",
                       "#8048a5",
                       "#633c7c",
                       "#420e66",
                       "#2c0845",
                       "#3e1957",
                       "#c076db",
                       "#a068b1",
                       "#52115d",
                       "#de9be3",
                       "#792b7c",
                       "#df73d1",
                       "#a73a9a",
                       "#8a3577",
                       "#681552",
                       "#c5c575",
                       "#b6c142",
                       "#828f29",
                       "#c1d168",
                       "#4663ab",
                       "#5b7fee",
                       "#8587de",
                       "#373e89",
                       "#182465",
                       "#5356bb",
                       "#251d77",
                       "#9778e2",
                       "#1c0b53",
                       "#311c66",
                       "#240568")

all_group_order <- c("A16S","ANT3-DPRIME","ANT6","ANT9","APH2-DPRIME","APH3-DPRIME","APH3-PRIME","APH6","RRS","RRSA","RRSC","RRSH",
                     "ROB","CMX","FLOR",
                     "CFR","ERM42","ERMA","ERMB","ERMF","ERMX","LNUA","LNUC","MEFA","MEFE","MLS23S","MPHE","MSRE","MYRA","VATE",
                     "FOLP","SULI","SULII",
                     "TET16S","TET32","TET33","TET40","TET44","TETH","TETM","TETO","TETQ","TETR","TETW","TETX","TETY","TETZ",
                     "ASR",
                     "TCRA",
                     "TCRY",
                     "TCRZ",
                     "QACE",
                     "QACEDELTA1",
                     "TUFAB",
                     "BLE",
                     "BRP",
                     "MERD",
                     "MERE",
                     "CFRA",
                     "EMRE",
                     "OPTRA",
                     "O23S")

ra_group_melt$Group <-factor(ra_group_melt$Group, levels=all_group_order)

ggplot(ra_Class_melt, aes(x= MEG_ID, y= Abundance, fill= Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=class_palette) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 10, face ="bold"),
        axis.text.x = element_blank())

#Need colors for
#A16S, ANT3-DPRIME,ANT6, ANT9, 

write.csv(ra_group_melt, "all_amr_groups.csv")
#Make palettes

viridis_pal(option="E")(63)

viridis_pal(option="E")(63)

#Cool site I found for this http://medialab.github.io/iwanthue/




tetracycline_group_all <- c("#8048a5",
                            "#633c7c",
                            "#420e66",
                            "#2c0845",
                            "#3e1957",
                            "#c076db",
                            "#a068b1",
                            "#52115d",
                            "#de9be3",
                            "#792b7c",
                            "#df73d1",
                            "#a73a9a",
                            "#8a3577",
                            "#681552")

sulfonamide_group_all <- c("#a27521",
                           "#dcab56",
                           "#ceaa34")

aminoglycoside_group_all <- c("#cb3f54",
                              "#df6f70",
                              "#8c282f",
                              "#d94d4d",
                              "#90281f",
                              "#d0735a",
                              "#b84828",
                              "#e08053",
                              "#cf6c2f",
                              "#7d3f13",
                              "#d38725",
                              "#ca9457")
  
mls_group_pal_all <- c("#ace268",
                       "#468027",
                       "#a1dd8a",
                       "#65b652",
                       "#296021",
                       "#569c56",
                       "#5bd581",
                       "#399f69",
                       "#6fdfa2",
                       "#48f2a7",
                       "#37a886",
                       "#55e5c1",
                       "#3be6ea",
                       "#5999e0",
                       "#4663ab")

phenicol_bl_group_pal_all <- c("#e16fbc",
                               "#b53c8f",
                               "#bb689c")

other_group_pal_all <- c("#c5c575",
                         "#b6c142",
                         "#828f29",
                         "#c1d168",
                         "#4663ab",
                         "#5b7fee",
                         "#8587de",
                         "#373e89",
                         "#182465",
                         "#5356bb",
                         "#251d77",
                         "#9778e2",
                         "#1c0b53",
                         "#311c66",
                         "#240568")

####HIGH####

sum(taxa_sums(ARG_resistome_high.ps)==0)
ARG_resistome_high.ps <- prune_taxa(taxa_sums(ARG_resistome_high.ps)>0, ARG_resistome_high.ps)
ARG_high.css <- phyloseq_transform_css(ARG_resistome_high.ps, log = F)
ARG_high.ra <- transform_sample_counts(ARG_high.css, function(x) {x/sum(x)}*100) #157 "taxa", 28 samples
rel_abund_high_trimmed <- transform_sample_counts(high_trimmed.css, function(x) {x/sum(x)}*100)

# agglomerate at different levels
ARG_high_Class.ra<- tax_glom(ARG_high.ra, taxrank = "Class") # 37 groups
ARG_high_Class_melt <- psmelt(ARG_high_Class.ra)


#subsetting from low palette, it has more
ra_group_palette_low



rel_abund_class_high <- tax_glom(rel_abund_high, taxrank = "Class") # 15 classes
ra_class_melt_high <- psmelt(rel_abund_class_high)
#ra_class_palette_high <- distinctColorPalette(15)--done, can redo if don't like

#Mech palette-remove # to see colors
##"#DDDA44" "#79E753" "#80B1DB" "#D8AC7C" "#BEE17B" "#C246E3" "#7BE2DA" "#D9DADD" "#69DA97" "#937ADB" "#CDE7B7" "#D5A0C5" "#809A8A" "#E066BD" "#D46358"

write.csv(otu_table(rel_abund_class_high),"high_class_otus.csv")
write.csv(tax_table(rel_abund_class_high),"high_class_taxa.csv")


abundance_Class_df <- ra_class_melt_high %>%
  group_by(Class)%>%
  reframe(mean=mean(Abundance))
  
View(abundance_Class_df)

sort(unique(ra_class_melt_high$Class))



high_ra_fill_order <- c("Aminoglycosides","betalactams","Copper_resistance","Drug_and_biocide_resistance",
                        "Elfamycins","Glycopeptides","Mercury_resistance","MLS","Multi-drug_resistance",
                        "Nickel_resistance","Nucleosides","Peroxide_resistance","Phenicol",
                        "Sulfonamides","Tetracyclines")



high_ra_fill_palette <- c("#E2E547","#80B2DB","#C79497","#73E3DC","#E0A84D","#DBBB8E","#DF69DA",
                          "#9147E0","#75E39E","#E0615B","#D5B4E3","#7B9F83",
                          "#D9679F","#C0D376","#D2DADF")

ra_class_melt_high$Class <- factor(ra_high_trimmed_group_melt$Class, levels = high_ra_fill_order)

ra_Class_high_plot <- ggplot(ra_class_melt_high, aes(x= MEG_ID, y= Abundance, fill= Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = ARG_high.dendro.data$labels$MEG_ID, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=high_ra_fill_palette) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        #legend.title = element_text(size=8),
        #legend.text = element_text(size = 6, face ="bold"),
        axis.text.x = element_blank())

stackeddendro_group_high_trimmed_plot <- ggarrange(ARG_high_dendro_plot, ra_Class_high_plot, ncol=1, heights = c(1.5,2))
##Tetracyclines
high_tetracylines <- subset_taxa(ra_high_trimmed_group, Class=="Tetracyclines")
high_tetracyclines_ra <- transform_sample_counts(high_tetracylines,function(x) {x/sum(x)}*100)
high_tetracyclines_ra_melt <- psmelt(high_tetracyclines_ra)


##Separating out by pool
#H-1

H12_1_tetracyclines_ra <-  high_tetracyclines_ra_melt[which(high_tetracyclines_ra_melt$H_12.01=="Y"),]
ggplot(H12_1_tetracyclines_ra, aes(x= MEG_ID, y= Abundance, fill= Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6, face ="bold"))

H12_2_tetracyclines_ra <-  high_tetracyclines_ra_melt[which(high_tetracyclines_ra_melt$H_12.02=="Y"),]
ggplot(H12_2_tetracyclines_ra, aes(x= MEG_ID, y= Abundance, fill= Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6, face ="bold"))

H12_3_tetracyclines_ra <-  high_tetracyclines_ra_melt[which(high_tetracyclines_ra_melt$H_12.03=="Y"),]
ggplot(H12_3_tetracyclines_ra, aes(x= MEG_ID, y= Abundance, fill= Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6, face ="bold"))

H12_4_tetracyclines_ra <-  high_tetracyclines_ra_melt[which(high_tetracyclines_ra_melt$H_12.04=="Y"),]
ggplot(H12_4_tetracyclines_ra, aes(x= MEG_ID, y= Abundance, fill= Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6, face ="bold"))

H12_5_tetracyclines_ra <-  high_tetracyclines_ra_melt[which(high_tetracyclines_ra_melt$H_12.05=="Y"),]
ggplot(H12_5_tetracyclines_ra, aes(x= MEG_ID, y= Abundance, fill= Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6, face ="bold"))

H12_6_tetracyclines_ra <-  high_tetracyclines_ra_melt[which(high_tetracyclines_ra_melt$H_12.06=="Y"),]
ggplot(H12_6_tetracyclines_ra, aes(x= MEG_ID, y= Abundance, fill= Group)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  #scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6, face ="bold"))


ggsave("stackeddendroclass_highplot.tiff", path = "~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/TE Pooling/Figures/Layered", plot =stackeddendro_group_high_trimmed_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)
ggsave("grouphigh_plot.tiff",path="~/Desktop",plot=ra_group_high_plot, device="tiff",dpi=600)

ra_high_trimmed_class <- tax_glom(rel_abund_high_trimmed, taxrank = "Class") #9 taxa
ra_high_trimmed_class_melt <- psmelt(ra_high_trimmed_class)

unique(ra_high_trimmed_class_melt$Class)

class_order <- c("Tetracyclines","Sulfonamides","Aminoglycosides","MLS","Phenicol","betalactams","Elfamycins","Drug_and_biocide_resistance","Mercury_resistance")

high_mean_sd <- ra_high_trimmed_class_melt %>%
  group_by(Number.in.Pool,Class)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))


high_class_plot <-ggplot(ra_high_trimmed_class_melt, aes(x= Class, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "ARG Classes", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = high_mean_sd, aes(y= meanRA, ymin=meanRA, ymax=meanRA+sd), linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  scale_x_discrete(limits = class_order, labels=c("Tetracyclines","Sulfonamides","Aminoglycosides","MLS","Phenicol","beta-lactams","Elfamycins","Drug and biocide resistance","Mercury resistance")) +
  scale_y_continuous(breaks=c(0,25,50,75,100))+
  scale_fill_manual(values =numberinpoolpalette) +
  scale_color_manual(values="black")+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle=45, vjust = 1, size=8, hjust =1))

#top groups
ra_high_trimmed_group_top10 <- merge_less_than_top_amr(ra_high_trimmed_group, top = 10)
ra_high_trimmed_group_10_melt <- psmelt(ra_high_trimmed_group_top10)

high_group_meansd <- ra_high_trimmed_group_10_melt %>%
  group_by(Group, Number.in.Pool)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))

high_group_order <- c("SULII","TETH","MLS23S","APH3-DPRIME","RRS")

high_group_plot <- ggplot(high_group_meansd, aes(x= Group, y= meanRA, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "ARG Classes", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = high_group_meansd, aes(y= meanRA, ymin=meanRA, ymax=meanRA+sd), linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  scale_x_discrete(limits = high_group_order) +
  scale_y_continuous(breaks=c(0,25,50,75,100), limits = c(0,50))+
  scale_fill_manual(values =numberinpoolpalette) +
  scale_color_manual(values="black")+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle=45, vjust = 1, size=8, hjust =1))

unique(ra_class_melt_high$Class)

high_amr_class_important <- subset(ra_high_trimmed_class_melt, Class=="Tetracyclines"|Class=="Phenicol"|Class=="betalactams"|Class=="MLS")

high_amr_class <- ggplot(high_amr_class_important, aes(x= Class, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "ARG Classes", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(stat="summary", linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  scale_x_discrete(limits = c("Tetracyclines","MLS","Phenicol","betalactams"), labels=c("Tetracycline","MLS","Phenicols","Beta-Lactams")) +
  scale_y_continuous(breaks=c(0,25,50,75,100), limits = c(0,50))+
  scale_fill_manual(values =numberinpoolpalette) +
  scale_color_manual(values="black")+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_blank())

low_amr_class_important <- subset(ra_low_trimmed_class_melt, Class=="Tetracyclines"|Class=="Phenicol"|Class=="betalactams"|Class=="MLS")

low_amr_class <- ggplot(low_amr_class_important, aes(x= Class, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "ARG Classes", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(stat="summary", linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  scale_x_discrete(limits = c("Tetracyclines","MLS","Phenicol","betalactams"), labels=c("Tetracycline","MLS","Phenicols","Beta-Lactams")) +
  scale_y_continuous(breaks=c(0,25,50,75,100), limits = c(0,50))+
  scale_fill_manual(values =poolonly_palette) +
  scale_color_manual(values="black")+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle=45, vjust = 1, size=20, hjust =1))

ggarrange(high_amr_class,low_amr_class, ncol = 1 )

#stats
###Class
#Tetracyclines
high_trimmed_Class_Tetracyclines <- subset(ra_high_trimmed_class_melt, Class=="Tetracyclines")
kruskal_test(high_trimmed_Class_Tetracyclines,Abundance~Number.in.Pool) #NS, p= 0.633
dunn_test(high_trimmed_Class_Tetracyclines, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_Tetracyclines %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (38.9 23.2), 3(26.2 12.4), 6(31.1 14.0), 12(32.0  5.37)
high_trimmed_Class_Tetracyclines %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #33.80996 17.89951

#Sulfonamides
high_trimmed_Class_Sulfonamides <- subset(ra_high_trimmed_class_melt, Class=="Sulfonamides")
kruskal_test(high_trimmed_Class_Sulfonamides,Abundance~Number.in.Pool) #NS, p= 0.892
dunn_test(high_trimmed_Class_Sulfonamides, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_Sulfonamides %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (26.0 14.2), 3(26.8 17.3), 6(22.4 15.5), 12(25.0  5.17)
high_trimmed_Class_Sulfonamides %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 25.46184 13.31587

#MLS
high_trimmed_Class_MLS <- subset(ra_high_trimmed_class_melt, Class=="MLS")
kruskal_test(high_trimmed_Class_MLS,Abundance~Number.in.Pool) #NS, p= 0.281
dunn_test(high_trimmed_Class_MLS, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_MLS %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (12.2  7.90), 3(16.6  8.20), 6(21.7 10.8), 12( 19.1  1.21)
high_trimmed_Class_MLS %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 15.72894 8.180008

#Aminoglycosides
high_trimmed_Class_Aminoglycosides <- subset(ra_high_trimmed_class_melt, Class=="Aminoglycosides")
kruskal_test(high_trimmed_Class_Aminoglycosides,Abundance~Number.in.Pool) #NS, p= 0.396
dunn_test(high_trimmed_Class_Aminoglycosides, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_Aminoglycosides %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (16.8  7.45), 3(21.4  7.42), 6(17.8  7.93), 12( 19.5  3.74)
high_trimmed_Class_Aminoglycosides %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 18.41293 6.894831

#Phenicol
high_trimmed_Class_Phenicol <- subset(ra_high_trimmed_class_melt, Class=="Phenicol")
kruskal_test(high_trimmed_Class_Phenicol,Abundance~Number.in.Pool) #NS, p= 0.11
dunn_test(high_trimmed_Class_Phenicol, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_Phenicol %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (3.30  6.26), 3(8.48  8.31), 6(6.19  9.16), 12(3.79  3.47)
high_trimmed_Class_Phenicol %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 4.910123 6.775416

#betalactams
high_trimmed_Class_betalactams <- subset(ra_high_trimmed_class_melt, Class=="betalactams")
kruskal_test(high_trimmed_Class_betalactams,Abundance~Number.in.Pool) #NS, p= 0.611
dunn_test(high_trimmed_Class_betalactams, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_betalactams %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.45  6.39), 3(0), 6(0), 12(0.142 0.317)
high_trimmed_Class_betalactams %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.164226 4.435726

#Elfamycins
high_trimmed_Class_Elfamycins <- subset(ra_high_trimmed_class_melt, Class=="Elfamycins")
kruskal_test(high_trimmed_Class_Elfamycins,Abundance~Number.in.Pool) #NS, p= 0.566
dunn_test(high_trimmed_Class_Elfamycins, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_Elfamycins %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.45  6.39), 3(0), 6(0), 12(0.142 0.317)
high_trimmed_Class_Elfamycins %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.164226 4.435726

#Drug_and_biocide_resistance
high_trimmed_Class_Drug_and_biocide_resistance <- subset(ra_high_trimmed_class_melt, Class=="Drug_and_biocide_resistance")
kruskal_test(high_trimmed_Class_Drug_and_biocide_resistance,Abundance~Number.in.Pool) #NS, p= 0.566
dunn_test(high_trimmed_Class_Drug_and_biocide_resistance, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Class_Drug_and_biocide_resistance %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.45  6.39), 3(0), 6(0), 12(0.142 0.317)
high_trimmed_Class_Drug_and_biocide_resistance %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.164226 4.435726

###Group
#APH3-DPRIME-7.5
high_trimmed_Group_APH3DPRIME <- subset(ra_high_trimmed_group_10_melt, Group=="APH3-DPRIME")
kruskal_test(high_trimmed_Group_APH3DPRIME,Abundance~Number.in.Pool) #NS, p= 0.676
dunn_test(high_trimmed_Group_APH3DPRIME, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_APH3DPRIME %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (6.23  6.40), 3(7.06  6.69), 6(8.38  6.38), 12(10.6   2.70)
high_trimmed_Group_APH3DPRIME %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #7.499475 5.902518

#APH6-3.9
high_trimmed_Group_APH6 <- subset(ra_high_trimmed_group_10_melt, Group=="APH6")
kruskal_test(high_trimmed_Group_APH6,Abundance~Number.in.Pool) #NS, p= 0.201
dunn_test(high_trimmed_Group_APH6, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_APH6 %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.66 3.43), 3(3.06 3.42), 6(6.19 5.01), 12(6.26 0.398)
high_trimmed_Group_APH6 %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3.893033 3.595099

#MLS23S-14.3
high_trimmed_Group_MLS23S <- subset(ra_high_trimmed_group_10_melt, Group=="MLS23S")
kruskal_test(high_trimmed_Group_MLS23S,Abundance~Number.in.Pool) #NS, p= 0.284
dunn_test(high_trimmed_Group_MLS23S, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_MLS23S %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (10.6  8.38), 3(15.5  7.67), 6(20.6 11.3), 12(17.4  3.61)
high_trimmed_Group_MLS23S %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #14.28524 8.531505

#ROB-3.5
high_trimmed_Group_ROB <- subset(ra_high_trimmed_group_10_melt, Group=="ROB")
kruskal_test(high_trimmed_Group_ROB,Abundance~Number.in.Pool) #NS, p= 0.611
dunn_test(high_trimmed_Group_ROB, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_ROB %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1(2.45 6.39), 3(0), 6(0), 12(0.142 0.317)
high_trimmed_Group_ROB %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3.497529 6.957579

#RRS-5.8
high_trimmed_Group_RRS <- subset(ra_high_trimmed_group_10_melt, Group=="RRS")
kruskal_test(high_trimmed_Group_RRS,Abundance~Number.in.Pool) #NS, p= 0.708
dunn_test(high_trimmed_Group_RRS, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_RRS %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1(6.51  8.90), 3(9.37 10.3), 6(3.06  4.32), 12(2.11  1.77)
high_trimmed_Group_RRS %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #5.84183 8.003189

#SULII-24.4
high_trimmed_Group_SULII <- subset(ra_high_trimmed_group_10_melt, Group=="SULII")
kruskal_test(high_trimmed_Group_SULII,Abundance~Number.in.Pool) #NS, p= 0.87
dunn_test(high_trimmed_Group_SULII, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_SULII %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (25.1 13.6), 3(26.3 17.3), 6(20.0 14.8), 12(24.1  4.88)
high_trimmed_Group_SULII %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #24.45871 13.02916

#TETH-24.4
high_trimmed_Group_TETH <- subset(ra_high_trimmed_group_10_melt, Group=="TETH")
kruskal_test(high_trimmed_Group_TETH,Abundance~Number.in.Pool) #NS, p= 0.731
dunn_test(high_trimmed_Group_TETH, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_TETH %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (27.5 22.1), 3(17.2 13.4), 6(23.2 14.4), 12(26.0  4.86)
high_trimmed_Group_TETH %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #24.40202 17.14364

#CMX-4.0
high_trimmed_Group_CMX <- subset(ra_high_trimmed_group_10_melt, Group=="CMX")
kruskal_test(high_trimmed_Group_CMX,Abundance~Number.in.Pool) #Sig, p= 0.0302
dunn_test(high_trimmed_Group_CMX, Abundance~Number.in.Pool, p.adjust.method = "BH") #1v3 different, 0.0273
high_trimmed_Group_CMX %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.32 5.63),3(8.48 8.31), 6(5.32 7.45), 12(2.03 0.898)
high_trimmed_Group_CMX %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #4.016333 6.321474

#TETW-3.3
high_trimmed_Group_TETW <- subset(ra_high_trimmed_group_10_melt, Group=="TETW")
kruskal_test(high_trimmed_Group_TETW,Abundance~Number.in.Pool) #NS, p= 0.844
dunn_test(high_trimmed_Group_TETW, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_TETW %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1(5.96  21.5), 3(0.832  2.04), 6(1.54   3.09), 12(0.256  0.573)
high_trimmed_Group_TETW %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3.213225 14.63729

#TETR-4.4
high_trimmed_Group_TETR <- subset(ra_high_trimmed_group_10_melt, Group=="TETR")
kruskal_test(high_trimmed_Group_TETR,Abundance~Number.in.Pool) #NS, p= 0.864
dunn_test(high_trimmed_Group_TETR, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
high_trimmed_Group_TETR %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (4.16 3.49),3(4.62 2.36), 6(4.29 2.15), 12(5.07 0.526)
high_trimmed_Group_TETR %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #4.439633 2.668023

arg_high_kw_stats <- ra_high_trimmed_group_melt %>%
  group_by(Group)%>%
  kruskal_test(Abundance~Number.in.Pool)
write.csv(arg_high_kw_stats,"arg_high_stats.csv")

arg_high_dunn_stats <- ra_high_trimmed_group_melt %>%
  group_by(Group)%>%
  dunn_test(Abundance~Number.in.Pool, p.adjust.method = "BH")

write.csv(arg_high_dunn_stats,"arg_high_stats_dunn.csv")

arg_high_group_number <- ra_high_trimmed_group_melt %>%
  group_by(Group,Number.in.Pool) %>%
  reframe(mean=mean(Abundance), sd=sd(Abundance))

write.csv(arg_high_group_number, "arg_high_group_abundance.csv")

ggplot(ra_high_trimmed_group_melt, aes(x=MEG_ID,y=Abundance, fill=Group))+
  geom_bar(stat="summary")

CMX_high <- subset_taxa(ra_high_trimmed_group, Group=="CMX")
CMX_high_melt <- psmelt(CMX_high)


ggplot(CMX_high_melt, aes(x=MEG_ID,y=Abundance, fill=Group))+
  geom_bar(stat="summary")+
  theme(axis.text.x.bottom = element_text(angle = 90))

####LOW####

rel_abund_low <- transform_sample_counts(low.css, function(x) {x/sum(x)}*100) #251 "taxa", 29 samples
rel_abund_low_trimmed <- transform_sample_counts(low_trimmed.css, function(x) {x/sum(x)}*100) #265 taxa, 17 samples

# agglomerate at different levels
rel_abund_Class_low <- tax_glom(rel_abund_low_trimmed, taxrank = "Class") # 13 classes
ra_class_melt_low <- psmelt(rel_abund_Class_low)

ra_low_trimmed_group <- tax_glom(rel_abund_low_trimmed, taxrank = "Group") #57 groups
ra_low_trimmed_group_melt <- psmelt(ra_low_trimmed_group)

length(unique(ra_class_melt_low$Class))
sort(unique(ra_class_melt_low$Class))
ra_group_palette_low <- distinctColorPalette(57)#--done, can redo if don't like
#ra_group_50_low <- prune_taxa(names(sort(taxa_sums(rel_abund_group_low),T)[1:50]), rel_abund_group_low)


write.csv(ra_low_trimmed_group_melt,"lowgroups.csv")



write.csv(otu_table(rel_abund_class_low),"low_class_otus.csv")
write.csv(tax_table(rel_abund_class_low),"low_class_taxa.csv")

abundance_group_low_df <- ra_low_trimmed_group_melt %>%
  group_by(Class,Group)%>%
  reframe(mean=mean(Abundance))

View(abundance_group_low_df)
  


low_class_order <- c("Acid_resistance","Aminoglycosides","betalactams","Copper_resistance","Drug_and_biocide_resistance",
                     "Elfamycins","Glycopeptides","MLS","Multi-drug_resistance","Oxazolidinone","Phenicol","Sulfonamides","Tetracyclines")


low_class_palette <- c("#D0E9BF","#E2E547","#80B2DB","#C79497","#73E3DC",
                       "#E0A84D","#DBBB8E","#9147E0","#75E39E","#7EE557",
                       "#D9679F","#C0D376","#D2DADF")

ra_class_melt_low$Class <-factor(ra_class_melt_low$Class, levels=low_class_order)

low_arg_abundance <- ra_low_trimmed_group_melt %>%
  group_by(Group, Number.in.Pool)%>%
  reframe(mean=mean(Abundance),sd=sd(Abundance))

write.csv(low_arg_abundance, "low_arg_abundance.csv")

ra_class_low_plot <- ggplot(ra_class_melt_low, aes(x= MEG_ID, y= Abundance, fill= Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = low_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values=low_class_palette) +
  #guides(fill=guide_legend(nrow=5), bycol=T)+
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=6),
        legend.text = element_text(size = 6, face ="bold"),
        axis.text.x = element_blank())

stackeddendro_class_low_trimmed_plot <- ggarrange(low_amr_dendro_plot, ra_class_low_plot, ncol=1,heights = c(30,65))


combined_dendros_ra <- ggarrange(stackeddendro_group_high_trimmed_plot,stackeddendro_group_low_trimmed_plot, ncol=1)

ggsave("stackeddendroclass_lowplot.tiff", path = "~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/TE Pooling/Figures/Layered", plot =stackeddendro_class_low_trimmed_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)
ggsave("group_lowplot.tiff",path = "~/Desktop",plot=ra_group_low_plot, device = "tiff", dpi=600)


ra_low_trimmed_class <- tax_glom(rel_abund_low_trimmed, taxrank = "Class") #13 taxa
ra_low_trimmed_class_melt <- psmelt(ra_low_trimmed_class)

unique(ra_low_trimmed_class_melt$Class)

ra_low_trimmed_gene <- tax_glom(rel_abund_low_trimmed, taxrank = "Gene")
ra_low_trimmed_gene_melt <- psmelt(ra_low_trimmed_gene)

#low_class_order <- c("Tetracyclines","Sulfonamides","MLS","Aminoglycosides","Phenicol","betalactams","Elfamycins","Drug_and_biocide_resistance","Mercury_resistance")

low_mean_sd <- ra_low_trimmed_class_melt %>%
  group_by(Number.in.Pool,Class)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))


low_class_plot <-ggplot(ra_low_trimmed_class_melt, aes(x= Class, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "ARG Classes", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = low_mean_sd, aes(y= meanRA, ymin=meanRA, ymax=meanRA+sd), linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  #scale_x_discrete(limits = class_order, labels=c("Tetracyclines","Sulfonamides","MLS","Aminoglycosides","Phenicol","beta-lactams","Elfamycins","Drug and biocide resistance","Mercury resistance")) +
  scale_y_continuous(breaks=c(0,25,50,75,100))+
  scale_fill_manual(values =poolonly_palette) +
  scale_color_manual(values="black")+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle=45, vjust = 1, size=8, hjust =1))

#top groups
ra_low_trimmed_group_top10 <- merge_less_than_top_amr(ra_low_trimmed_group, top = 10)
ra_low_trimmed_group_10_melt <- psmelt(ra_low_trimmed_group_top10)

unique(ra_low_trimmed_group_10_melt$Group)


low_group_meansd <- ra_low_trimmed_group_10_melt %>%
  group_by(Group, Number.in.Pool)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))

low_group_order <- c("TETH","SULII","MLS23S","APH3-DPRIME","TETM")

lowgroup_plot <- ggplot(low_group_meansd, aes(x= Group, y= meanRA, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "ARG Classes", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = low_group_meansd, aes(y= meanRA, ymin=meanRA, ymax=meanRA+sd), linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  scale_x_discrete(limits = low_group_order) +
  scale_y_continuous(breaks=c(0,25,50,75,100), limits = c(0,50))+
  scale_fill_manual(values =poolonly_palette) +
  scale_color_manual(values="black")+
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle=45, vjust = 1, size=8, hjust =1))



#stats
###Class
#Tetracyclines
low_trimmed_Class_Tetracyclines <- subset(ra_low_trimmed_class_melt, Class=="Tetracyclines")
kruskal_test(low_trimmed_Class_Tetracyclines,Abundance~Number.in.Pool) #NS, p= 0.633
dunn_test(low_trimmed_Class_Tetracyclines, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_Tetracyclines %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(48.1  17.0), 6(43.6  18.2), 12(36.9  11.7)
low_trimmed_Class_Tetracyclines %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #42.56868 15.50709

#Sulfonamides
low_trimmed_Class_Sulfonamides <- subset(ra_low_trimmed_class_melt, Class=="Sulfonamides")
kruskal_test(low_trimmed_Class_Sulfonamides,Abundance~Number.in.Pool) #NS, p= 0.892
dunn_test(low_trimmed_Class_Sulfonamides, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_Sulfonamides %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (26.0 14.2), 3(26.8 17.3), 6(22.4 15.5), 12(25.0  5.17)
low_trimmed_Class_Sulfonamides %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 25.46184 13.31587

#MLS
low_trimmed_Class_MLS <- subset(ra_low_trimmed_class_melt, Class=="MLS")
kruskal_test(low_trimmed_Class_MLS,Abundance~Number.in.Pool) #NS, p= 0.281
dunn_test(low_trimmed_Class_MLS, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_MLS %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (12.2  7.90), 3(16.6  8.20), 6(21.7 10.8), 12( 19.1  1.21)
low_trimmed_Class_MLS %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 15.72894 8.180008

#Aminoglycosides
low_trimmed_Class_Aminoglycosides <- subset(ra_low_trimmed_class_melt, Class=="Aminoglycosides")
kruskal_test(low_trimmed_Class_Aminoglycosides,Abundance~Number.in.Pool) #NS, p= 0.215
dunn_test(low_trimmed_Class_Aminoglycosides, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_Aminoglycosides %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(16.0  5.48), 6(10.8  4.80), 12(14.4  4.34)
low_trimmed_Class_Aminoglycosides %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 13.61443 5.061697

#Phenicol
low_trimmed_Class_Phenicol <- subset(ra_low_trimmed_class_melt, Class=="Phenicol")
kruskal_test(low_trimmed_Class_Phenicol,Abundance~Number.in.Pool) #NS, p= 0.538
dunn_test(low_trimmed_Class_Phenicol, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_Phenicol %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(1.32   1.53), 6(0.689  1.34), 12(3.14   3.90)
low_trimmed_Class_Phenicol %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.741743 2.664016

#betalactams
low_trimmed_Class_betalactams <- subset(ra_low_trimmed_class_melt, Class=="betalactams")
kruskal_test(low_trimmed_Class_betalactams,Abundance~Number.in.Pool) #NS, p= 0.0673
dunn_test(low_trimmed_Class_betalactams, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_betalactams %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(0), 6(8.15 10.3), 12(1.76  2.73)
low_trimmed_Class_betalactams %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.164226 4.435726

#Elfamycins
low_trimmed_Class_Elfamycins <- subset(ra_low_trimmed_class_melt, Class=="Elfamycins")
kruskal_test(low_trimmed_Class_Elfamycins,Abundance~Number.in.Pool) #NS, p= 0.566
dunn_test(low_trimmed_Class_Elfamycins, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_Elfamycins %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.45  6.39), 3(0), 6(0), 12(0.142 0.317)
low_trimmed_Class_Elfamycins %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.164226 4.435726

#Drug_and_biocide_resistance
low_trimmed_Class_Drug_and_biocide_resistance <- subset(ra_low_trimmed_class_melt, Class=="Drug_and_biocide_resistance")
kruskal_test(low_trimmed_Class_Drug_and_biocide_resistance,Abundance~Number.in.Pool) #NS, p= 0.566
dunn_test(low_trimmed_Class_Drug_and_biocide_resistance, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Class_Drug_and_biocide_resistance %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (2.45  6.39), 3(0), 6(0), 12(0.142 0.317)
low_trimmed_Class_Drug_and_biocide_resistance %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) # 1.164226 4.435726

###Group
#APH3-DPRIME-6.3
low_trimmed_Group_APH3DPRIME <- subset(ra_low_trimmed_group_10_melt, Group=="APH3-DPRIME")
kruskal_test(low_trimmed_Group_APH3DPRIME,Abundance~Number.in.Pool) #NS, p= 0.539
dunn_test(low_trimmed_Group_APH3DPRIME, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_APH3DPRIME %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(7.15  4.15), 6(4.81  3.84), 12(7.08  3.04)
low_trimmed_Group_APH3DPRIME %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #6.299662 3.618852

#APH6-3.1
low_trimmed_Group_APH6 <- subset(ra_low_trimmed_group_10_melt, Group=="APH6")
kruskal_test(low_trimmed_Group_APH6,Abundance~Number.in.Pool) #NS, p= 0.46
dunn_test(low_trimmed_Group_APH6, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_APH6 %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(2.81  2.54), 6(2.44  2.36), 12(4.12  2.56)
low_trimmed_Group_APH6 %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3.139477 2.443875

#MLS23S-12.8
low_trimmed_Group_MLS23S <- subset(ra_low_trimmed_group_10_melt, Group=="MLS23S")
kruskal_test(low_trimmed_Group_MLS23S,Abundance~Number.in.Pool) #NS, p= 0.375
dunn_test(low_trimmed_Group_MLS23S, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_MLS23S %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(9.96  6.73), 6(14.7   8.09), 12(13.4   6.68)
low_trimmed_Group_MLS23S %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #12.84701 7.05049

#ROB-3.5
low_trimmed_Group_ROB <- subset(ra_low_trimmed_group_10_melt, Group=="ROB")
kruskal_test(low_trimmed_Group_ROB,Abundance~Number.in.Pool) #NS, p= 0.0673
dunn_test(low_trimmed_Group_ROB, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_ROB %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(0), 6(8.15 10.3 ), 12(1.76  2.73)
low_trimmed_Group_ROB %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3.497529 6.957579

#RRS-2.3
low_trimmed_Group_RRS <- subset(ra_low_trimmed_group_10_melt, Group=="RRS")
kruskal_test(low_trimmed_Group_RRS,Abundance~Number.in.Pool) #NS, p= 0.551
dunn_test(low_trimmed_Group_RRS, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_RRS %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(3.10  2.13), 6(1.95  2.57), 12(1.95  1.68)
low_trimmed_Group_RRS %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #2.287171 2.090505

#SULII-21.3
low_trimmed_Group_SULII <- subset(ra_low_trimmed_group_10_melt, Group=="SULII")
kruskal_test(low_trimmed_Group_SULII,Abundance~Number.in.Pool) #NS, p= 0.551
dunn_test(low_trimmed_Group_SULII, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_SULII %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(20.0 14.2 ), 6(17.3 13.3), 12(26.3  8.55)
low_trimmed_Group_SULII %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #21.27242 12.03348

#TETH-24.7
low_trimmed_Group_TETH <- subset(ra_low_trimmed_group_10_melt, Group=="TETH")
kruskal_test(low_trimmed_Group_TETH,Abundance~Number.in.Pool) #NS, p= 0.71
dunn_test(low_trimmed_Group_TETH, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_TETH %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(23.8  9.51), 6(25.2 13.2), 12(24.9  5.59)
low_trimmed_Group_TETH %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #24.67101 9.347528

#TETM-5.0
low_trimmed_Group_TETM <- subset(ra_low_trimmed_group_10_melt, Group=="TETM")
kruskal_test(low_trimmed_Group_TETM,Abundance~Number.in.Pool) #NS, p= 0.424
dunn_test(low_trimmed_Group_TETM, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_TETM %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(6.86  7.39), 6(3.55  8.69), 12(4.94 12.1)
low_trimmed_Group_TETM %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #5.012132 9.211883

#TETO-3.3
low_trimmed_Group_TETO <- subset(ra_low_trimmed_group_10_melt, Group=="TETO")
kruskal_test(low_trimmed_Group_TETO,Abundance~Number.in.Pool) #NS, p= 0.932
dunn_test(low_trimmed_Group_TETO, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_TETO %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(5.28  10.7), 6(4.46  10.8), 12(0.609  0.968)
low_trimmed_Group_TETO %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3.340319 8.364782

#TETR-4.9
low_trimmed_Group_TETR <- subset(ra_low_trimmed_group_10_melt, Group=="TETR")
kruskal_test(low_trimmed_Group_TETR,Abundance~Number.in.Pool) #NS, p= 0.401
dunn_test(low_trimmed_Group_TETR, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
low_trimmed_Group_TETR %>%
  group_by(Number.in.Pool)%>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #3(4.62  1.06), 6(5.46  2.32), 12(4.62  2.32)
low_trimmed_Group_TETR %>%
  reframe(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #4.915973 1.955736

majorgroups_combined <- ggarrange(high_group_plot, lowgroup_plot,  ncol = 1, nrow = 2, labels = "AUTO")

ggsave("amr_groups_plot.tiff",plot=majorgroups_combined,device = "tiff", path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/ForPublication/Figures", units = "mm", height=200, width = 85, dpi=600)


low_ra_stats_kw <- ra_low_trimmed_group_melt %>%
  group_by(Group)%>%
  kruskal_test(Abundance~Number.in.Pool)

write.csv(low_ra_stats_kw, "low_arg_stats_kw.csv")

#REL ABUND MH####



high_mh_alpha_div <- estimate_richness(high_mh, measures = c("Observed","Shannon","Simpson","InvSimpson"))
high_alpha_div.df <- as(sample_data(high), "data.frame")
high_alpha_div_meta <- cbind(high_alpha_div, high_alpha_div.df)
high_alpha_div_meta


sample_data(data_mh)
str(sample_data(data_mh))
high_mh <- subset_samples(data_mh, PrevH=="Y"&MEG_ID!="SC_193")
sum(sample_sums(high_mh)==0)
sum(taxa_sums(high_mh)==0) #No taxa just found in high

high_mh_melt <- psmelt(high_mh)

high_mh_species <- tax_glom(high_mh, taxrank = "Species")
high_mh_species_melt <- psmelt(high_mh_species)
high_mh_strain <- subset_taxa(high_mh, Strain!="unclassified")
high_mh_strain_melt <- tax_glom(high_mh_strain,taxrank = "Strain") %>%
  psmelt()
high_mh_strain_species <- tax_glom(high_mh_strain, taxrank = "Species") %>%
  psmelt()

high_strain_abundance_table <- high_mh_strain_species %>%
  group_by(Number.in.Pool,Species)%>%
  reframe(mean=mean(Abundance), sd=sd(Abundance))

mh_strain_palette <- c("#E39DDE","#CCE7E0","#DA855C","#CFC794","#D8BFCE",
                       "#66D0D0","#913EE0","#84AADE","#E15480","#966794",
                       "#8CDA96","#7CE656","#D4DA55","#DC59D7","#8575DC")

ggplot(high_mh_species_melt, aes(x= MEG_ID, y= Abundance, fill= Species)) +
  labs(y= "Normalized Abundance", title = "Species") +
  #geom_errorbar(data = drug_ra_class_low_windiv, aes(x= Class, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary") +
  #scale_x_discrete(limits = phyla_abundance_order) +
  #scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  #scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, size=14, hjust =1))

ggplot(high_mh_species_melt, aes(x= Number.in.Pool, y= Abundance, fill= Species)) +
  labs(y= "Normalized Abundance", title = "Classification") +
  #geom_errorbar(data = drug_ra_class_low_windiv, aes(x= Class, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary") +
  #scale_x_discrete(limits = phyla_abundance_order) +
  #scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  #scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme(#legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_text(size=24),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, size=14, hjust =1))
sort(unique(high_mh_strain_melt$MEG_ID))

unique(high_mh_strain_melt$MEG_ID) 

high_indi_order <- c("SC_142","SC_147","SC_148","SC_151","SC_162","SC_163","SC_171","SC_180",
                     "SC_181","SC_183","SC_187","SC_194","SC_197","SC_202","SC_207","SC_209",
                     "SC_215","SC_218","SC_220","SC_230","SC_233","SC_234","SC_237",
                     "HP_03-01","HP_03-02","HP_03-03","HP_03-04","HP_03-05","HP_03-06","HP_06-01",
                     "HP_06-02","HP_06-03","HP_06-04","HP_06-05","HP_06-06","H_12-01","H_12-02","H_12-03","H_12-04","H_12-05","H_12-06")

high_indi <- ggplot(high_mh_strain_melt, aes(x= MEG_ID, y= Abundance, fill= Strain)) +
  theme_bw()+
  labs(y= "Normalized Abundance", title = "MH Strains-High") +
  #geom_errorbar(data = drug_ra_class_low_windiv, aes(x= Class, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary", color="black") +
  scale_fill_manual(values =mh_strain_palette) +
  scale_x_discrete(limits=high_indi_order)+
  scale_y_continuous(limits = c(0,250))+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text( angle=90 , size=8, face="bold"))

high_group <- ggplot(high_mh_strain_melt, aes(x= Number.in.Pool)) +
  labs(y= "Normalized Abundance", x="Pool Size") +
  geom_errorbar(data = high_strain_abundance_table, aes(x= Number.in.Pool, y=mean,ymin= mean, ymax=mean+sd),  width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary", aes(fill=Strain, y=Abundance), color="black") +
  #scale_x_discrete(limits = c(1,3,6,12)) +
  scale_fill_manual(values =mh_strain_palette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  scale_y_continuous(limits = c(0,250), )+
  theme_bw()+
  theme(#legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8))



high_mh_plot <- high_indi + 
  annotation_custom(
    ggplotGrob(high_group), 
    xmin = 24, xmax = 36, ymin = 100, ymax = 250
  )

ggsave("high_mh_strain_indi.tiff", high_indi, device = "tiff", path = "~/Desktop", units = "mm", width = 180, height = 100, dpi = 600)
ggsave("high_mh_strain_group.tiff", high_group, device = "tiff", path = "~/Desktop", units = "mm", width = 40, height = 80, dpi = 600)

sort(unique(high_mh_strain_melt$Strain))

high_mh_Brain2012 <- subset(high_mh_strain_melt, Strain=="MhBrain2012")

strainabundance_df <- high_mh_strain_melt %>%
  group_by(Strain, Number.in.Pool)%>%
  reframe(mean=mean(Abundance), sd=sd(Abundance))

write.csv(strainabundance_df, "mh_abundance.csv")

kruskal_test_mh_strain_df <- high_mh_strain_melt %>%
  group_by(Strain)%>%
  kruskal_test(Abundance~Number.in.Pool)
write.csv(kruskal_test_mh_strain_df,"strain_stats.csv")

dunn_test_mh_strain_df <- high_mh_strain_melt %>%
  group_by(Strain)%>%
  dunn_test(Abundance~Number.in.Pool,p.adjust.method = "BH")
write.csv(dunn_test_mh_strain_df,"strain_stats_dunn.csv")

kruskal_test(high_mh_Brain2012, Abundance~Number.in.Pool)
dunn_test(high_mh_Brain2012, Abundance~Number.in.Pool,p.adjust.method = "BH")

high_mh_D193 <- subset(high_mh_strain_melt, Strain=="D193")

kruskal_test(high_mh_D193, Abundance~Number.in.Pool) #P=0.06
dunn_test(high_mh_D193, Abundance~Number.in.Pool,p.adjust.method = "BH") 

high_mh_Brain2012 %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Abundance))

#REL ABUND MH####
sample_data(data_mh)
str(sample_data(data_mh))
tax_table(data_mh)
low_mh <- subset_samples(data_mh, PrevL=="Y")
low_mh <- subset_samples(low_mh, Number.in.Pool!=1)
sum(sample_sums(low_mh)==0)
sum(taxa_sums(low_mh)==0) #No taxa just found in low


low_mh_species <- tax_glom(low_mh, taxrank = "Species")
low_mh_strain <- subset_taxa(low_mh, Strain!="unclassified")
low_mh_strain_melt <- tax_glom(low_mh_strain,taxrank = "Strain") %>%
  psmelt()
low_mh_species_melt <- psmelt(low_mh_species)

ggplot(low_mh_species_melt, aes(x= MEG_ID, y= Abundance, fill= Species)) +
  labs(y= "Normalized Abundance") +
  #geom_errorbar(data = drug_ra_class_low_windiv, aes(x= Class, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary") +
  #scale_x_discrete(limits = phyla_abundance_order) +
  #scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  #scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_text(size=24),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, size=14, hjust =1))

ggplot(low_mh_species_melt, aes(x= Number.in.Pool, y= Abundance, fill= Species)) +
  labs(y= "Normalized Abundance", title = "Classification") +
  #geom_errorbar(data = drug_ra_class_low_windiv, aes(x= Class, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary") +
  #scale_x_discrete(limits = phyla_abundance_order) +
  #scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  #scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme(#legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    title = element_text(size=24),
    axis.line.y = element_line(size = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", size = 0.75),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, vjust = 1, size=14, hjust =1))

low_mh_strain_species <- tax_glom(low_mh_strain, taxrank = "Species") %>%
  psmelt()

low_strain_abundance_table <- low_mh_strain_species %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Abundance), sd=sd(Abundance))

sort(unique(low_mh_strain_melt$MEG_ID))
low_indi_order <- c("LP_03-01","LP_03-02","LP_03-03","LP_03-04","LP_03-05","LP_03-06",
                    "LP_06-01","LP_06-02","LP_06-03","LP_06-04","LP_06-05","LP_06-06",
                    "L_12-01","L_12-02", "L_12-03","L_12-04" , "L_12-05" , "L_12-06")

low_indi <- ggplot(low_mh_strain_melt, aes(x= MEG_ID, y= Abundance, fill= Strain)) +
  theme_bw()+
  labs(y= "Normalized Abundance", title = "MH Strains-low") +
  #geom_errorbar(data = drug_ra_class_low_windiv, aes(x= Class, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary", color="black") +
  scale_x_discrete(limits = low_indi_order) +
  scale_fill_manual(values =mh_strain_palette) +
  scale_y_continuous(limits = c(0,250))+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

low_group <- ggplot(low_mh_strain_melt, aes(x= Number.in.Pool)) +
  labs(y= "Normalized Abundance", title = "MH Strains-low") +
  geom_errorbar(data = low_strain_abundance_table, aes(x= Number.in.Pool, y=mean,ymin= mean, ymax=mean+sd),  width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary", aes(fill=Strain, y=Abundance), color="black") +
  #scale_x_discrete(limits = c(1,3,6,12)) +
  scale_fill_manual(values =mh_strain_palette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  scale_y_continuous(limits = c(0,250))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8))

ggsave("low_mh_strain_indi.tiff", low_indi, device = "tiff", path = "~/Desktop", units = "mm", width = 100, height = 80, dpi = 600)
ggsave("low_mh_strain_group.tiff", low_group, device = "tiff", path = "~/Desktop", units = "mm", width = 20, height = 80, dpi = 600)


ggarrange(low_indi, low_group,nrow=1)

MHRAplots

lowMhBrain <- subset(low_mh_strain_melt, Strain=="MhBrain2012")

kruskal_test(lowMhBrain, Abundance~Number.in.Pool)
dunn_test(lowMhBrain, Abundance~Number.in.Pool, p.adjust.method = "BH")

lowMhBrain %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Abundance),sd=sd(Abundance))


MHStrain_plot <- ggarrange(high_indi,high_group,low_indi,low_group, nrow = 2, ncol = 2)

low_mh_strain_abundance_df <- low_mh_strain_melt %>%
  group_by(Strain,Number.in.Pool)%>%
  reframe(mean=mean(Abundance), sd=sd(Abundance))

write.csv(low_mh_strain_abundance_df, "mh_abundance_low.csv")

kruskal_test_mh_strain_low_df <- low_mh_strain_melt %>%
  group_by(Strain)%>%
  kruskal_test(Abundance~Number.in.Pool)
write.csv(kruskal_test_mh_strain_low_df,"strain_stats_low.csv")

dunn_test_mh_strain_low_df <- low_mh_strain_melt %>%
  group_by(Strain)%>%
  dunn_test(Abundance~Number.in.Pool,p.adjust.method = "BH")
write.csv(dunn_test_mh_strain_low_df,"strain_stats_dunn_low.csv")

ggsave("mhstrain_plot.tif", plot=MHStrain_plot, device = "tiff", units="mm", height = 100, width=180, dpi = 600)

##Individuals only:
indi_mh_ra <- transform_sample_counts(indi_mh, function(x) {x/sum(x)}*100)
indi_mh_strain <- tax_glom(indi_mh, taxrank = "Strain")
indi_mh_strain_only <- subset_taxa(indi_mh_strain, Strain!="unclassified")
indi_mh_strain_only <- transform_sample_counts(indi_mh_strain_only, function(x) {x/sum(x)}*100)
indi_mh_strain_melt <- psmelt(indi_mh_strain_only)

indi_mh_ra_strain <- tax_glom(indi_mh_ra, taxrank = "Strain")
indi_mh_ra_strain_only <- subset_taxa(indi_mh_ra_strain, Strain!="unclassified")
indi_mh_ra_strain_melt <- psmelt(indi_mh_ra_strain_only)

indi_mh_strain_melt$Mh.Cult <- factor(indi_mh_strain_melt$Mh.Cult, levels=c("Y","N"))

unique(indi_mh_strain_melt$Strain)

strain_abundance <- indi_mh_strain_melt %>%
  group_by(Strain) %>%
  reframe(mean=mean(Abundance))

sort(strain_abundance$mean)

indi_mh_strain_melt$Strain <- factor(indi_mh_strain_melt$Strain, levels=c("USDA-ARS-USMARC-184","serotype A2 str. BOVINE", "M42548" , "35","D193", "serotype 6 str. H23" ,"USMARC_2286", "USDA-ARS-USMARC-185","MhSwine2000" ,"D153","serotype A2 str. OVINE","D38","D174","D171", "MhBrain2012" ))

indi_mh_norm <- ggplot(indi_mh_strain_melt, aes(x= MEG_ID)) +
  labs(y= "Relative Abundance", title = "MH Strains-Individuals", shape="MH Culture", fill="MH Culture") +
  geom_bar(stat = "summary", aes( y= Abundance, fill= Strain), color="black") +
  geom_point(data=indi_mh_alpha_div_meta, aes(y=-6, shape=Mh.Cult, color=Mh.Cult), size=4)+
  geom_text(data=indi_mh_alpha_div_meta, aes(y=-6, label=Observed), size=3)+
  scale_fill_manual(values =mh_strain_palette, limits = c("USDA-ARS-USMARC-184","serotype A2 str. BOVINE", "M42548" , "35","D193", "serotype 6 str. H23" ,"USMARC_2286", "USDA-ARS-USMARC-185","MhSwine2000" ,"D153","serotype A2 str. OVINE","D38","D174","D171", "MhBrain2012" )) +
  #facet_wrap(~Mh.Cult,scales = "free_x")+
  #scale_y_continuous(limits = c(-10,300))+
  theme_bw()+
  theme(
    legend.text = element_text(size=6),
    legend.key.size = unit(4,'mm'),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    title = element_blank(),
    axis.line.y = element_line(size = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", size = 0.75),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, vjust = 1, size=8, hjust =1))

indi_mh_ra <- ggplot(indi_mh_ra_strain_melt, aes(x= MEG_ID)) +
  labs(y= "Relative Abundance (%)", title = "MH Strains-Individuals", shape="MH Culture", fill="MH Culture") +
  geom_bar(stat = "summary", aes( y= Abundance, fill= Strain), color="black") +
  geom_point(data=indi_mh_alpha_div_meta, aes(y=3, shape=Mh.Cult, color=Mh.Cult), size=4)+
  geom_text(data=indi_mh_alpha_div_meta, aes(y=3, label=Observed), size=3)+
  scale_fill_manual(values =mh_strain_palette) +
  #facet_wrap(~Mh.Cult,scales = "free_x")+
  #scale_y_continuous(limits = c(0,0.5,1,1.5,2, 2.5), labels = c(0,0.5,1,1.5,2, 2.5))+
  theme_bw()+
  theme(
    legend.text = element_text(size=6),
    legend.key.size = unit(4,'mm'),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    title = element_blank(),
    axis.line.y = element_line(size = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", size = 0.75),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, vjust = 1, size=8, hjust =1))


wilcox_test(indi_mh_alpha_div_meta, Observed~Mh.Cult)
indi_mh_alpha_div_meta$Mh.Cult <- factor(indi_mh_alpha_div_meta$Mh.Cult, levels = c("Y","N"))

observedplot <- ggplot(data = indi_mh_alpha_div_meta, aes(x=Mh.Cult, y=Observed, color=Mh.Cult, fill=Mh.Cult, alpha=Mh.Cult))+
  theme_bw()+
  labs(y="Observed MH Strains", x="MH Culture", color="MH Culture", fill="MH Culture", alpha="MH Culture", shape="MH Culture")+
  geom_boxplot()+
  geom_point(aes(shape=Mh.Cult))+
  scale_alpha_manual(values=c(0.7,0.7))+
  scale_y_continuous(limits=c(0,17))+
  stat_compare_means(label = "p.signif", label.x.npc = 0.5)+
  theme(legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    title = element_blank(),
    axis.line.y = element_line(size = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", size = 0.75),
    axis.title.y = element_text(size = 10),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text())

indi_mh_alpha_div_meta %>%
  group_by(Mh.Cult) %>%
  reframe(medianOb=median(Observed),meanOb=mean(Observed), sdOb=sd(Observed), nOb=n())%>%
  mutate(seOb=sdOb/sqrt(nOb), lower_ci=meanOb-qt(1-(0.05/2),nOb-1)*seOb,
         upper_ci=meanOb+qt(1-(0.05/2),nOb-1)*seOb)

hs_alpha_data <-high_mh_alpha_div_meta %>%
  group_by(Number.in.Pool)%>%
  reframe(median=median(Observed))

hs_mh <- ggplot(high_mh_alpha_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Richness") +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  geom_text(data = hs_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits=c(3,20))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

kruskal_test(high_mh_alpha_div_meta, Observed~Number.in.Pool) #P=0.0833, NS
dunn_test(high_mh_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH")

mh_alpha_div_meta %>%
  group_by(Number.in.Pool)%>%
  wilcox_test(Observed~MHqPCR)

mh_alpha_div_meta %>%
  group_by(Number.in.Pool)%>%
  wilcox_test(MH_align~ MHqPCR)

ggplot(mh_alpha_div_meta, aes(x= MHqPCR, y = Observed, fill =  MHqPCR, colour =  MHqPCR)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Richness") +
  #scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  #geom_text(data = hs_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  #scale_y_continuous(limits=c(3,20))+
  facet_wrap(~Number.in.Pool)+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggplot(mh_alpha_div_meta, aes(x= MHqPCR, y = Observed, fill =  MHqPCR, colour =  MHqPCR)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Richness") +
  #scale_fill_manual(values = numberinpoolpalette) +
  #scale_colour_manual(values = numberinpoolpalette) +
  #geom_text(data = hs_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  #scale_y_continuous(limits=c(3,20))+
  #facet_wrap(~Number.in.Pool)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggplot(mh_alpha_div_meta, aes(x= MHqPCR, y = MH_align, fill =  MHqPCR, colour =  MHqPCR)) + 
    theme_bw() + 
    labs(title= "high_trimmed", y= "Richness") +
    #scale_fill_manual(values = numberinpoolpalette) +
    #scale_colour_manual(values = numberinpoolpalette) +
    #geom_text(data = hs_alpha_data, aes(y=18,label=median), color="black", size=6)+
    geom_boxplot(alpha = 0.5, size = 1) +
    geom_point(size = 2.5) +
    #scale_y_continuous(limits=c(3,20))+
    facet_wrap(~Number.in.Pool)+
    theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())


kruskal_test(low_mh_alpha_div_meta, Observed~Number.in.Pool) #NS, P=0.607
dunn_test(low_mh_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH")

ls_alpha_data <-low_mh_alpha_div_meta %>%
  group_by(Number.in.Pool)%>%
  reframe(median=median(Observed))

ls_mh <- ggplot(low_mh_alpha_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "low_trimmed", y= "Richness") +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  geom_text(data = ls_alpha_data, aes(y=18,label=median), color="black", size=6)+
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits=c(3,20))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 18, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggarrange(hs_mh, ls_mh,  ncol=2)


ggarrange(indi_mh_ra,observedplot, ncol = 1)

inset_plot<-indi_mh_ra + 
  annotation_custom(
    ggplotGrob(observedplot), 
    xmin = 18, xmax = 24, ymin = 0.8, ymax = 2.5
  )
ggsave("mhstrain_insetl.tif",plot=inset_plot,device = "tiff", path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/ForPublication/Figures", units = "mm", height=150, width = 180, dpi=300)

##PCR status
data_mh_strain <- tax_glom(data_mh, taxrank = "Strain")
data_mh_strain <- subset_taxa(data_mh_strain, Strain!="unclassified")
data_mh_strain_melt <- psmelt(data_mh_strain)

data_mh_strain_melt$MHqPCR <- factor(data_mh_strain_melt$MHqPCR, levels = c("YES","NO"))

ggplot(data_mh_strain_melt, aes(x= Strain, y=Abundance, fill=MHqPCR, alpha=MHqPCR)) +
  labs(y= "Normalized Abundance", title = "MH Strains-low") +
  geom_errorbar(stat = "summary" , width = 0.55, position = position_dodge(0.92))+
  geom_bar(stat = "summary", color="black", position=position_dodge(0.92)) +
  #scale_x_discrete(limits = c("YES","NO")) +
  scale_fill_manual(values =mh_strain_palette, limits=c("YES","NO")) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  facet_wrap(~Number.in.Pool)+
  #scale_y_continuous(limits = c(0,250))+
  theme_bw()+
  theme(#legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, angle = 90, vjust=-0.5))

data_mh_strain_melt %>%
  group_by(Strain,Number.in.Pool) %>%
  wilcox_test(Abundance~MHqPCR)


