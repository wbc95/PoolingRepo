###LoadLibraries###
setwd("~/path/to/working/directory")


### #load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr);
library(ggplot2); library(btools); library(randomcoloR); library(cowplot);
library(pairwiseAdonis); library(picante); library(gridExtra); library(grid); library(wrapr);
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(rstatix); library(viridis); library(ggprism); library(ggpubr);
library(UpSetR); library(vegan);library(BiodiversityR); library(randomForest);
library(MicEco);library(ggbreak); library(tidyr); library(DAtest); library(ANCOMBC); library(knitr); library(ggrepel)

#### source some functions ####
source("~/change16STaxaNames.R")
source("~/w_unifrac.R")
source("~/g_unifrac.R")
source("~/uw_unifrac.R")
source("~/MergeLowAbund.R")
source("~/removeNARows.R")


####ImportData####
###16S####
qiimedata <- import_biom("~/16S/table-with-taxonomy.biom", "~/16S/tree.nwk", "~/16S/dna-sequences.fasta")
metadata_16S <- import_qiime_sample_data("~/16Smetadata_complete.txt") 

data_16S <- merge_phyloseq(qiimedata,metadata_16S) #163 samples, 86100 ASVs
sum(taxa_sums(data_16S)==0) # 22504 taxa with no counts
data_16S <- prune_taxa(taxa_sums(data_16S) > 0, data_16S) #63596 taxa and 163 samples
sum(sample_sums(data_16S)==0)##None

### # check the names of our ranks
rank_names(data_16S) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data_16S)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data_16S) # beauty, now they are named properly

### # changing the GG style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data_16S)) # extract the taxonomy table as a data frame
tax.data.names <- change16Staxa(tax.data) # this gets rid of the GG format

### # now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

### # now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- kingdom
  } else if (tax.data.names[i,3] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- phylum
  } else if (tax.data.names[i,4] == ""){
    class <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- class
  } else if (tax.data.names[i,5] == ""){
    order <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- order
  } else if (tax.data.names[i,6] == ""){
    family <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- family
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Species[i] <- paste("unclassified ",tax.data.names$Genus[i], sep = "")
  }
}

head(tax.data.names) # great, no more NAs and no more k__
tax_table(data_16S) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(data_16S), 20) # "unclassified Unassigned" is that what we want?
#heading names changed

sample_names(data_16S) <- sample_data(data_16S)$MEG_ID
sample_names(data_16S)
sample_data(data_16S)$sample.id <- sample_data(data_16S)$MEG_ID

data_16S

data_16S <- subset_taxa(data_16S, Kingdom=="Bacteria"|Kingdom=="Unassigned"|Kingdom=="Archaea")
sum(taxa_sums(data_16S)==0)
sum(sample_sums(data_16S)==0)
min(sample_sums(data_16S))
max(sample_sums(data_16S))

##Pruning####
data_16S_50k <- subset_samples(data_16S, sample_sums(data_16S)>50000)
data_16S_50k #158 samples, 5 dropped; 63496 Taxa
sample_data(data_16S_50k)$Number.in.Pool <- as.factor(sample_data(data_16S_50k)$Number.in.Pool)
sort(sample_sums(data_16S_50k))


sum(taxa_sums(data_16S_50k)==0) #42 taxa
data_16S_50k <- prune_taxa(taxa_sums(data_16S_50k)>0, data_16S_50k)
any(sample_sums(data_16S_50k)==0)
data_16S_50k <- subset_taxa(data_16S_50k, Kingdom=="Bacteria"|Kingdom=="Unassigned"|Kingdom=="Archaea")

data_16S_50k_kingdom <- tax_glom(data_16S_50k, taxrank = "Kingdom", NArm = F)
kingdom_melt <- psmelt(data_16S_50k_kingdom)
unique(kingdom_melt$Kingdom)
##

indiv_16S <- subset_samples(data_16S_50k, Type=="Individual")
sum(taxa_sums(indiv_16S)==0)
indiv_16S <- prune_taxa(taxa_sums(indiv_16S)>0, indiv_16S)
sample_data(indiv_16S)
min(sample_sums(indiv_16S))



pool_16S <- subset_samples(data_16S_50k, Type=="Pool")
sum(taxa_sums(pool_16S)==0)
pool_16S <- prune_taxa(taxa_sums(pool_16S)>0, pool_16S)

pool3_16S <- subset_samples(pool_16S,Number.in.Pool==3)
sum(taxa_sums(pool3_16S)==0)
pool3_16S <- prune_taxa(taxa_sums(pool3_16S)>0, pool3_16S)
min(sample_sums(pool3_16S))

pool6_16S <-subset_samples(pool_16S, Number.in.Pool==6)
min(sample_sums(pool6_16S))

pool12_16S <- subset_samples(pool_16S, Number.in.Pool==12)
min(sample_sums(pool12_16S))

sample_sum_df_50k <- data.frame(ASVs = sample_sums(data_16S_50k))
sample_sum_meta <- as(sample_data(data_16S_50k), "data.frame") # making into DF for metadata
sample_sum_50k_stats <- cbind(sample_sum_df_50k, sample_sum_meta)
kruskal_test(sample_sum_50k_stats, ASVs~Number.in.Pool)
dunn_test(sample_sum_50k_stats, ASVs~Number.in.Pool, p.adjust.method = "BH")

sample_sum_50k_stats %>%
  group_by(Number.in.Pool) %>%
  summarise(min=min(ASVs), med=median(ASVs), max=max(ASVs))

sample_sum_50k_stats$Number.in.Pool <- as.factor(sample_sum_50k_stats$Number.in.Pool)

###AMR####
TE_data <- import_biom("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/AMR/metaSNV.biom")
sample_names(TE_data)
write.csv(sample_names(TE_data),"~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/samplenamesformetadata_allTE.csv")

te_mapfile <- import_qiime_sample_data("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/TEmetadata_SNV.txt")

ARGdata <- merge_phyloseq(TE_data,te_mapfile) #17265 taxa, 60 samples
sample_data(ARGdata)
sample_names(ARGdata) <- sample_data(ARGdata)$samplename

rank_names(ARGdata)
head(tax_table(ARGdata))
colnames(tax_table(ARGdata)) <- c("Type","Class","Mechanism","Group","SNV","SNP")
rank_names(ARGdata)
ARGdata_noSNP <- subset_taxa(ARGdata, SNP=="N")
ARGdata_noSNP
sum(taxa_sums(ARGdata_noSNP)==0)
sum(sample_sums(ARGdata_noSNP)==0)

min(sample_sums(ARGdata_noSNP)) #203
max(sample_sums(ARGdata_noSNP)) #37893
mean(sample_sums(ARGdata_noSNP)) #4934.3
median(sample_sums(ARGdata_noSNP)) #2255.5
sort(sample_sums(ARGdata_noSNP)) #Not going to remove any?


ARG_trimmed <- subset_samples(ARGdata_noSNP, sample_sums(ARGdata_noSNP)>1000)
ARG_trimmed <- prune_taxa(taxa_sums(ARG_trimmed)>0,ARG_trimmed)

sample_data(ARG_trimmed)$Number.in.Pool <- factor(sample_data(ARG_trimmed)$Number.in.Pool, levels = c("1","3","6","12"))

ARG_trimmed_data <- as(sample_data(ARG_trimmed),"data.frame")

t.test(ARG_trimmed_data$InputReads)

median(ARG_trimmed_data$InputReads)

IQR(ARG_trimmed_data$InputReads)

IQR(ARG_trimmed_data$MH_align)
IQR(ARG_trimmed_data$MHPerc)


##FInidng overall median and 95% CI of generated reads
reads <- ARG_trimmed_data$InputReads
n <- length(reads)
sample_median <- median(reads)
se <- sd(reads)/sqrt(n)

t_crit <- qt(0.975, df = n-1)
low_bound <- sample_median - (t_crit * se)
up_bound <- sample_median + (t_crit * se)

MHreads <- ARG_trimmed_data$MH_align
n <- length(MHreads)
sample_median <- median(MHreads)
se <- sd(MHreads)/sqrt(n)

t_crit <- qt(0.975, df = n-1)
low_bound <- sample_median - (t_crit * se)
up_bound <- sample_median + (t_crit * se)

MHreads_perc <- ARG_trimmed_data$MHPerc
n <- length(MHreads_perc)
sample_median <- median(MHreads_perc)
se <- sd(MHreads_perc)/sqrt(n)

t_crit <- qt(0.975, df = n-1)
low_bound <- sample_median - (t_crit * se)
up_bound <- sample_median + (t_crit * se)

ARGcounts <- ARGcount_data$ARGs
n=length(ARGcounts)
sample_median <- median(ARGcounts)
se <- sd(ARGcounts)/sqrt(n)
t_crit <- qt(0.975, df = n-1)
low_bound <- sample_median - (t_crit*se)
up_bound <- sample_median+(t_crit*se)

aggregate(InputReads ~ Number.in.Pool, data = ARG_trimmed_data, FUN = function(x) {
  n <- length(x)
  mean_val <- mean(x)
  med_val <- median(x)
  se <- sd(x) / sqrt(n)
  
  # qt() finds the critical t-value for a 95% two-tailed test
  t_val <- qt(0.975, df = n - 1) 
  
  c(
    Mean = mean_val,
    Median = med_val,
    Lower = med_val - (t_val * se),
    Upper = med_val + (t_val * se)
  )
})

head(sample_data(ARGdata_noSNP),5)
sample_data(ARGdata_noSNP)$Number.in.Pool <- factor(sample_data(ARGdata_noSNP)$Number.in.Pool, levels = c("1","3","6","12"))

ARG_indiv <- subset_samples(ARGdata_noSNP, SampleType=="Individual")
any(taxa_sums(ARG_indiv)==0)
ARG_indiv <- prune_taxa(taxa_sums(ARG_indiv)>0, ARG_indiv)

ARG_pool <- subset_samples(ARGdata_noSNP, SampleType=="Pool")
any(taxa_sums(ARG_pool)==0)
ARG_pool <- prune_taxa(taxa_sums(ARG_pool)>0, ARG_pool)

ARG_high <- subset_samples(ARGdata_noSNP,PrevH=="Y")
any(taxa_sums(ARG_high)==0)
any(sample_sums(ARG_high)==0)
ARG_high <- prune_taxa(taxa_sums(ARG_high)>0,ARG_high)
ARG_high_trimmed <- subset_samples(ARG_high, samplename!="HP_06-03")
sum(taxa_sums(ARG_high_trimmed)==0)#131 SNVs only seen in this sample
ARG_high_trimmed <- prune_taxa(taxa_sums(ARG_high_trimmed)>0, ARG_high_trimmed)
sample_names(ARG_high_trimmed)

sort(sample_sums(ARG_high_trimmed))


ARGcount_data <- data.frame(ARGs=sample_sums(ARG_trimmed))

IQR(ARGcount_data$ARGs)

View(sample_data(ARG_high))
ARG_high_63 <-subset_samples(ARGdata_noSNP, HP_06.03=="Y" & Number.in.Pool !=3)
ARG_high_63 <- prune_taxa(taxa_sums(ARG_high_63)>0, ARG_high_63)
high_alpha_div_63 <- estimate_richness(ARG_high_63, measures=c("Observed"))
high_63_alpha.df <- as(sample_data(ARG_high_63),"data.frame")
alpha_div_high_63 <- cbind(high_alpha_div_63,high_63_alpha.df)

wilcox_test(alpha_div_high_63, Observed~Number.in.Pool)


ARG_low <- subset_samples(ARGdata_noSNP,Prevalence=="Low")
any(taxa_sums(ARG_low)==0)
any(sample_sums(ARG_low)==0)
ARG_low <- prune_taxa(taxa_sums(ARG_low)>0,ARG_low)


head(tax_table(ARGdata_noSNP))
ARG_noSNP_ANCOMBC <- ARGdata_noSNP
colnames(tax_table(ARG_noSNP_ANCOMBC)) <-c("Phylum","Class","Order","Genus","Species","SNP")




###CHeck that one was resequenced
#ARG_6_3 <- subset_samples(ARGdata_noSNP, HP_06.03=="Y")
#ARG_6_3 <- prune_taxa(taxa_sums(ARG_6_3)>0,ARG_6_3)

#ARG_63.css <- phyloseq_transform_css(ARG_6_3, norm = T, log = F)
#ARG_63.ra <- transform_sample_counts(ARG_63.css, function(x) {x/sum(x)}*100)
#ARG_63_Group <- tax_glom(ARG_63.ra,taxrank = "Group")

#ARG_63_Group_melt <- psmelt(ARG_63_Group)




#####MH####
library(ggsci) # if not using this package, make sure to switch out the palatte in the figures
library(metagMisc)
library(dplyr)
library(ggdendro)
library(cowplot)

MH_data <- import_biom("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/MH/MHclass.biom")
write.csv(sample_names(MH_data),"~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/samplenamesformetadata_MH.csv")

MH_mapfile <- import_qiime_sample_data("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/MH/TEmetadata_MH.txt")

MHdata <- merge_phyloseq(MH_data,MH_mapfile) #561 taxa, 60 samples
sample_data(MHdata)
sample_names(MHdata) <- sample_data(MHdata)$samplename

rank_names(MHdata)
head(tax_table(MHdata))
colnames(tax_table(MHdata)) <- c("Genus","Species","PSV","All")

PSVdata <- subset_taxa(MHdata, PSV!="")
PSVdata <-subset_taxa(PSVdata, Species=="haemolytica")
PSVdata <- subset_samples(PSVdata, samplename!="HP_06-03")
sum(taxa_sums(PSVdata)==0)
sample_data(mh_TE_PSV_transformed.ps)$Number.in.Pool <- factor(sample_data(mh_TE_PSV_transformed.ps)$Number.in.Pool, levels=c("1","3","6","12"))
View(tax_table(mh_TE_PSV_transformed.ps))

MH_high <- subset_samples(mh_TE_PSV_transformed.ps, PrevH=="Y")
any(taxa_sums(MH_high)==0)
MH_high <- prune_taxa(taxa_sums(MH_high)>0,MH_high)
MH_low <- subset_samples(mh_TE_PSV_transformed.ps, Prevalence=="Low")
any(taxa_sums(MH_low)==0)
MH_low <- prune_taxa(taxa_sums(MH_low)>0,MH_low)

MH_indiv <- subset_samples(PSVdata,SampleType=="Individual")
any(taxa_sums(MH_indiv)==0)
MH_pool <- subset_samples(PSVdata,SampleType!="Individual")
any(taxa_sums(MH_pool)==0)


numberinpoolpalette <- c("#30123BFF","#1AE4B6FF" ,"#FABA39FF" ,"#7A0403FF")
poolonly_palette<- c("#1AE4B6FF" ,"#FABA39FF" ,"#7A0403FF")

poolcomparisons <- list(c("1","3"),c("1","6"),c("1","12"),c("3","6"),c("3","12"),c("6","12"))
poolonly_comparisons <- list(c("3","6"),c("3","12"),c("6","12"))

MHmetadata.df <- as(sample_data(PSVdata), "data.frame")
MHmetadata.df %>%
  reframe(mean=mean(MH_align), median=median(MH_align), n=n(),sd=sd(MH_align), min=min(MH_align),max=max(MH_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
MHmetadata.df %>%
  reframe(mean=mean(PSV_align), median=median(PSV_align), n=n(),sd=sd(PSV_align), min=min(PSV_align),max=max(PSV_align))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
kruskal_test(MHhigh.df,MH_align~Number.in.Pool)
kruskal_test(MHhigh.df,PSV_align~Number.in.Pool)
kruskal_test(MHlow.df, MH_align~Number.in.Pool)
kruskal_test(MHlow.df, PSV_align~Number.in.Pool)

MHmetadata.df %>%
  reframe(mean=mean(MH_per), median=median(MH_per), n=n(),sd=sd(MH_per), min=min(MH_per),max=max(MH_per))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
MHmetadata.df %>%
  reframe(mean=mean(PSV_per), median=median(PSV_per), n=n(),sd=sd(PSV_per), min=min(PSV_per),max=max(PSV_per))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
kruskal_test(MHhigh.df,MH_per~Number.in.Pool)
kruskal_test(MHhigh.df,PSV_per~Number.in.Pool)
kruskal_test(MHlow.df, MH_per~Number.in.Pool)
kruskal_test(MHlow.df, PSV_per~Number.in.Pool)


