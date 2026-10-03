###Script from from LivAbs_16S (Lee)--indicated by ### before

### #load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr);
library(ggplot2); library(btools); library(randomcoloR); library(cowplot);
library(pairwiseAdonis); library(picante); library(gridExtra); library(grid); library(wrapr);
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(rstatix); library(viridis); library(ggprism); library(ggpubr);
library(UpSetR); library(vegan);library(BiodiversityR); library(randomForest); #library(MicrobiotaProcess)
library(MicEco);library(ggbreak)


#Set working directory
setwd("~/path/working/directory")

### # source some stuff
source("g_unifrac.R")
source("w_unifrac.R")
source("uw_unifrac.R")
source("change16STaxaNames.R")
source("~/MergeLowAbund.R")
source("~/removeNARows.R")

###Scripts from Lee

### # import data
qiimedata <- import_biom("table-with-taxonomy.biom", "tree.nwk", "dna-sequences.fasta")
map_file <- import_qiime_sample_data("16Smetadata_complete.txt") # need to convert the date to date type

View(map_file)

sample_names(qiimedata)

write.csv(sample_names(qiimedata), "samplenamesformetadata.csv")

# combining sample data with the rest
data <- merge_phyloseq(qiimedata,map_file)



sample_names(data) <- sample_data(data)$MEG_ID
sample_names(data)

##Data exploration
data #163 samples with 86100 ASVS
sum(taxa_sums(data)==0)  # 22504 taxa with no counts
data <- prune_taxa(taxa_sums(data) > 0, data)
sum(sample_sums(data)==0)


##### # check the names of our ranks
View(tax_table(data))
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) # beauty, now they are named properly
head(tax_table(data))


### # changing the GG style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- change16Staxa(tax.data) # this gets rid of the GG format
head(tax_table(data))

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
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
tail(tax_table(data), 20) # sweet, lookin good!
head(tax_table(data), 20)

data # 163 samples, 63596 ASVs
any(taxa_sums(data)==0)
any(sample_sums(data)==0)

data <- subset_taxa(data, Kingdom!="Eukaryota")
data #163 samples, 63590 ASVs, lost 6 taxa

data <- subset_taxa(data,Kingdom!="Stramenopiles") #163 samples, 63584 ASVs, lost 6 more taxa

any(taxa_sums(data)==0)



##Read count visualization
#all
sample_sum_df <- data.frame(sum = sample_sums(data))
readplot_all <- ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill="skyblue", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) #  Talk to Lee about distribution
readplot_all

###Trimming
#### some QC checks
#all\
data
min(sample_sums(data)) #357
max(sample_sums(data)) #7194819
mean(sample_sums(data)) #422551.6
median(sample_sums(data)) #273121
sort(sample_sums(data)) #357, 4555, 5222, 12130,... trim at 50000 or 10000

write.csv(sample_sums(data),"ASVcount.csv")

data_50k <- subset_samples(data, sample_sums(data)>50000)
data_50k #158 samples, 5 dropped
sample_data(data_50k)$Number.in.Pool <- as.factor(sample_data(data_50k)$Number.in.Pool)
sort(sample_sums(data_50k))


sum(taxa_sums(data_50k)==0) #42 taxa
data_50k <- prune_taxa(taxa_sums(data_50k)>0, data_50k)
any(sample_sums(data_50k)==0)
data_50k_kingdom <- tax_glom(data_50k, taxrank ="Kingdom")

data_50k_kingdom_melt <- psmelt(data_50k_kingdom)
unique(data_50k_kingdom_melt$Kingdom) # "Bacteria"   "Unassigned" "Archaea"   OK, good

sample_sum_df_50k <- data.frame(ASVs = sample_sums(data_50k))
readplot_all_50k <- ggplot(sample_sum_df_50k, aes(x = ASVs)) + 
  geom_histogram(color = "black", fill="skyblue", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) #  Talk to Lee about distribution
readplot_all_50k

numberinpoolpalette <- viridis_pal(option = "H")(4) # "#30123BFF" "#1AE4B6FF" "#FABA39FF" "#7A0403FF"
poolonlypalette <- c("#1AE4B6FF", "#FABA39FF", "#7A0403FF")

sample_data(data_50k)$Number.in.Pool <- as.factor(sample_data(data_50k)$Number.in.Pool)


##Split all data by prevalence (thinking about this as just pools, not worried about how they were selected)
#high
high_all<-subset_samples(data_50k, PrevH=="Y")
sum(taxa_sums(high_all)==0) #38262
high_all <-prune_taxa(taxa_sums(high_all)>0,high_all)
sum(sample_sums(high_all)==0)

#low
low_all<-subset_samples(data_50k, PrevL=="Y")
sum(taxa_sums(low_all)==0) #9285
low_all <-prune_taxa(taxa_sums(low_all)>0,low_all)
sum(sample_sums(low_all)==0)


##Split up pooling strategies

#subset
subsetpooling <- subset_samples(data_50k, PoolingType!="Random") #137 samples
sum(taxa_sums(subsetpooling)==0)
subsetpooling <- prune_taxa(taxa_sums(subsetpooling)>0, subsetpooling)


##HP
subsethigh <- subset_samples(subsetpooling,PrevH=="Y") #88 samples
sum(taxa_sums(subsethigh)==0)
subsethigh <- prune_taxa(taxa_sums(subsethigh)>0, subsethigh)
sample_data(subsethigh)$Number.in.Pool <- as.factor(sample_data(subsethigh)$Number.in.Pool)



##LP
subsetlow <- subset_samples(subsetpooling,PrevL=="Y") #88 samples
sum(taxa_sums(subsetlow)==0)
subsetlow <- prune_taxa(taxa_sums(subsetlow)>0, subsetlow)

#random
randompooling <- subset_samples(data_50k,PoolingType!="Subset") #135 samples
sum(taxa_sums(randompooling)==0)
randompooling <- prune_taxa(taxa_sums(randompooling)>0, randompooling)

##HP
randomhigh <- subset_samples(randompooling,PrevH=="Y") #88 samples
sum(taxa_sums(randomhigh)==0)
randomhigh <- prune_taxa(taxa_sums(randomhigh)>0, randomhigh)

##LP
randomlow <- subset_samples(randompooling,PrevL=="Y") #88 samples
sum(taxa_sums(randomlow)==0)
randomlow <- prune_taxa(taxa_sums(randomlow)>0, randomlow)

#both for comparison
pools_3v6 <- subset_samples(data_50k, PoolingType=="Subset"|PoolingType=="Random") #44 samples
pools_3v6 <- prune_taxa(taxa_sums(pools_3v6)>0, pools_3v6)

##HP
pools_3v6high <- subset_samples(pools_3v6,PrevH=="Y") #20 samples
pools_3v6high <- prune_taxa(taxa_sums(pools_3v6high)>0,pools_3v6high)

##LP
pools_3v6low <- subset_samples(pools_3v6,PrevL=="Y") #24 samples
pools_3v6low <- prune_taxa(taxa_sums(pools_3v6low)>0,pools_3v6low)





##Just compare high and low individuals
View(sample_data(data_50k))

##high##
high_individual <- subset_samples(data_50k, PrevH=="Y"&PrevL=="N"&Number.in.Pool=="1") #32, 10964 taxa
sum(taxa_sums(high_individual)==0)


##low##
low_individual <- subset_samples(data_50k, PrevH=="N"&PrevL=="Y"&Number.in.Pool=="1")
sum(taxa_sums(low_individual)==0)

individual_hl<- merge_phyloseq(high_individual,low_individual)
sum(taxa_sums(individual_hl)==0)
individual_hl <- prune_taxa(taxa_sums(individual_hl)>0, individual_hl) #63 samples, 31415 taxa

high_individual<- prune_taxa(taxa_sums(high_individual)>0,high_individual)
low_individual<- prune_taxa(taxa_sums(low_individual)>0,low_individual) #31 sample, 25575 taxa






poolcomparisons = list(c("1","3"),c("1","6"),c("1","12"),c("3","6"),c("3","12"),c("6","12"))

#For Mannheimia culture status
#all
any(taxa_sums(data_50k)==0)

data_50k_metadata <- as(sample_data(data_50k),"data.frame")
seqdepth_50k <- cbind(sample_sum_df_50k,data_50k_metadata)
View(seqdepth_50k)

seqdepth_50k$Number.in.Pool <- as.factor(seqdepth_50k$Number.in.Pool)

kruskal_test(seqdepth_50k, ASVs~Number.in.Pool)
wilcox_test(seqdepth_50k,ASVs~Contains.CP)
seqdepth_50k %>%
  group_by(Number.in.Pool)%>%
  wilcox_test(ASVs~Contains.CP)

t.test(seqdepth_50k$ASVs)

ggplot(data=seqdepth_50k, aes(x=Number.in.Pool, y=ASVs, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  facet_wrap(~Contains.CP)+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)


aggregate(ASVs ~ Number.in.Pool, data = seqdepth_50k, FUN = function(x) {
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

#all
#high
sample_sum_df_high <- data.frame(ASVs = sample_sums(high_all))
high_all_metadata <- as(sample_data(high_all),"data.frame")
seqdepth_high <- cbind(sample_sum_df_high, high_all_metadata)

View(seqdepth_high)

kruskal_test(seqdepth_high,ASVs~Number.in.Pool) #0.00000398, Sig
dunn_test(seqdepth_high,ASVs~Number.in.Pool, p.adjust.method="BH") #1 different from 3 and 12

ggplot(data=seqdepth_high, aes(x=Number.in.Pool, y=ASVs, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)

#low
sample_sum_df_low <- data.frame(ASVs = sample_sums(low_all))
low_all_metadata <- as(sample_data(low_all),"data.frame")
seqdepth_low <- cbind(sample_sum_df_low, low_all_metadata)
View(seqdepth_low)

kruskal_test(seqdepth_low,ASVs~Number.in.Pool) #P<0.00001, Sig
dunn_test(seqdepth_low,ASVs~Number.in.Pool, p.adjust.method="BH") #1 different from all

ggplot(data=seqdepth_low, aes(x=Number.in.Pool, y=ASVs, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)


#subset
##high
sample_sum_subhigh_df <- data.frame(ASV_count = sample_sums(subsethigh))
metadata_subhigh_df <- as(sample_data(subsethigh), "data.frame" )
seqdepth_subhigh_df <- cbind(sample_sum_subhigh_df,metadata_subhigh_df)

seqdepth_subhigh_df$Number.in.Pool <- as.factor(seqdepth_subhigh_df$Number.in.Pool)
View(seqdepth_subhigh_df)

seqdepth_subhigh_df %>%
  group_by(Number.in.Pool) %>%
  reframe(mean=mean(ASV_count),median=median(ASV_count),sd=sd(ASV_count), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
       upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(seqdepth_subhigh_df,ASV_count~Number.in.Pool)
dunn_test(seqdepth_subhigh_df, ASV_count~Number.in.Pool, p.adjust.method = "BH")
ggplot(data=seqdepth_subhigh_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)
#Significant difference in read depth, can't evaluate alpha diversity, can compare pools to each other
 
##low
sample_sum_sublow_df <- data.frame(ASV_count = sample_sums(subsetlow))
metadata_sublow_df <- as(sample_data(subsetlow), "data.frame" )
seqdepth_sublow_df <- cbind(sample_sum_sublow_df,metadata_sublow_df)

seqdepth_sublow_df$Number.in.Pool <- as.factor(seqdepth_sublow_df$Number.in.Pool)

seqdepth_sublow_df %>%
  group_by(Number.in.Pool) %>%
  reframe(mean=mean(ASV_count),median=median(ASV_count),sd=sd(ASV_count), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(seqdepth_sublow_df,ASV_count~Number.in.Pool)
dunn_test(seqdepth_sublow_df, ASV_count~Number.in.Pool, p.adjust.method = "BH")
ggplot(data=seqdepth_sublow_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))+
  scale_color_manual(values = numberinpoolpalette)  #Significant difference in read depth, can't evaluate alpha diversity, can compare pools to each other

#random#
##high
sample_sum_randhigh_df <- data.frame(ASV_count = sample_sums(randomhigh))
metadata_randhigh_df <- as(sample_data(randomhigh), "data.frame" )
seqdepth_randhigh_df <- cbind(sample_sum_randhigh_df,metadata_randhigh_df)

seqdepth_randhigh_df$Number.in.Pool <- as.factor(seqdepth_randhigh_df$Number.in.Pool)

seqdepth_randhigh_df %>%
  group_by(Number.in.Pool) %>%
  reframe(mean=mean(ASV_count),median=median(ASV_count),sd=sd(ASV_count), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(seqdepth_randhigh_df, ASV_count~Number.in.Pool)
dunn_test(seqdepth_randhigh_df, ASV_count~Number.in.Pool, p.adjust.method = "BH")



ggplot(data=seqdepth_randhigh_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)+
  stat_compare_means()+
  stat_compare_means(comparisons = poolcomparisons) #Significant difference in read depth, can't evaluate alpha diversity, can compare pools to each other

##low
sample_sum_randlow_df <- data.frame(ASV_count = sample_sums(randomlow))
metadata_randlow_df <- as(sample_data(randomlow), "data.frame" )
seqdepth_randlow_df <- cbind(sample_sum_randlow_df,metadata_randlow_df)

seqdepth_randlow_df$Number.in.Pool <- as.factor(seqdepth_randlow_df$Number.in.Pool)

seqdepth_randlow_df %>%
  group_by(Number.in.Pool) %>%
  reframe(mean=mean(ASV_count),median=median(ASV_count),sd=sd(ASV_count), n=n())%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

kruskal_test(seqdepth_randlow_df, ASV_count~Number.in.Pool)
dunn_test(seqdepth_randlow_df, ASV_count~Number.in.Pool, p.adjust.method = "BH")

ggplot(data=seqdepth_randlow_df, aes(x=Number.in.Pool, y=ASV_count, fill=Number.in.Pool,color=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))+
  scale_color_manual(values = numberinpoolpalette)+
 stat_compare_means()+
 stat_compare_means(comparisons = poolcomparisons) ##Basically individual sequencing depth lower than pools; can compare pools to each other



###Pull out just pools
##all
#high
high_pools <- subset_samples(high_all,Number.in.Pool!=1)
sum(taxa_sums(high_pools)==0)
high_pools <- prune_taxa(taxa_sums(high_pools)>0,high_pools)

sample_sum_high_pool_df <- data.frame(ASV_count = sample_sums(high_pools))
metadata_high_pool_df <- as(sample_data(high_pools), "data.frame" )
seqdepth_high_pool_df <- cbind(sample_sum_high_pool_df,metadata_high_pool_df)

kruskal_test(seqdepth_high_pool_df,ASV_count~Number.in.Pool)
ggplot(data=seqdepth_high_pool_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = poolonlypalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values =poolonlypalette)

#low
low_pools <- subset_samples(low_all,Number.in.Pool!=1)
sum(taxa_sums(low_pools)==0)
low_pools <- prune_taxa(taxa_sums(low_pools)>0,low_pools)

sample_sum_low_pool_df <- data.frame(ASV_count = sample_sums(low_pools))
metadata_low_pool_df <- as(sample_data(low_pools), "data.frame" )
seqdepth_low_pool_df <- cbind(sample_sum_low_pool_df,metadata_low_pool_df)

kruskal_test(seqdepth_low_pool_df,ASV_count~Number.in.Pool)
ggplot(data=seqdepth_low_pool_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = poolonlypalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values =poolonlypalette)

#low
low_pools <- subset_samples(low_all,Number.in.Pool!=1)
sum(taxa_sums(low_pools)==0)
low_pools <- prune_taxa(taxa_sums(low_pools)>0,low_pools)

##subset
#high
hs_pools <- subset_samples(subsethigh, Number.in.Pool!="1")
sum(taxa_sums(hs_pools)==0)
hs_pools <- prune_taxa(taxa_sums(hs_pools)>0, hs_pools)

sample_sum_hspool_df <- data.frame(ASV_count = sample_sums(hs_pools))
metadata_hspool_df <- as(sample_data(hs_pools), "data.frame" )
seqdepth_hspool_df <- cbind(sample_sum_hspool_df,metadata_hspool_df)

seqdepth_hspool_df$Number.in.Pool <- as.factor(seqdepth_hspool_df$Number.in.Pool)

ggplot(data=seqdepth_hspool_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = poolonlypalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values =poolonlypalette)+
  stat_compare_means() ##no difference

#low
ls_pools <- subset_samples(subsetlow, Number.in.Pool!="1")
sum(taxa_sums(ls_pools)==0)
ls_pools <- prune_taxa(taxa_sums(ls_pools)>0, ls_pools)

sample_sum_lspool_df <- data.frame(ASV_count = sample_sums(ls_pools))
metadata_lspool_df <- as(sample_data(ls_pools), "data.frame" )
seqdepth_lspool_df <- cbind(sample_sum_lspool_df,metadata_lspool_df)

seqdepth_lspool_df$Number.in.Pool <- as.factor(seqdepth_lspool_df$Number.in.Pool)

ggplot(data=seqdepth_lspool_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = poolonlypalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values =poolonlypalette)+
  stat_compare_means()

##random
#high
hr_pools <- subset_samples(randomhigh, Number.in.Pool!="1")
sum(taxa_sums(hr_pools)==0)
hr_pools <- prune_taxa(taxa_sums(hr_pools)>0, hr_pools)

sample_sum_hrpool_df <- data.frame(ASV_count = sample_sums(hr_pools))
metadata_hrpool_df <- as(sample_data(hr_pools), "data.frame" )
seqdepth_hrpool_df <- cbind(sample_sum_hrpool_df,metadata_hrpool_df)

seqdepth_hrpool_df$Number.in.Pool <- as.factor(seqdepth_hrpool_df$Number.in.Pool)

ggplot(data=seqdepth_hrpool_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = poolonlypalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values =poolonlypalette)+
  stat_compare_means() ##no difference

#low
lr_pools <- subset_samples(randomlow, Number.in.Pool!="1")
sum(taxa_sums(lr_pools)==0)
lr_pools <- prune_taxa(taxa_sums(lr_pools)>0, lr_pools)

sample_sum_lrpool_df <- data.frame(ASV_count = sample_sums(lr_pools))
metadata_lrpool_df <- as(sample_data(lr_pools), "data.frame" )
seqdepth_lrpool_df <- cbind(sample_sum_lrpool_df,metadata_lrpool_df)

seqdepth_lrpool_df$Number.in.Pool <- as.factor(seqdepth_lrpool_df$Number.in.Pool)

ggplot(data=seqdepth_lrpool_df, aes(x=Number.in.Pool, y=ASV_count, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  scale_fill_manual(values = poolonlypalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values =poolonlypalette)+
  stat_compare_means()

##No difference in any pools in sequencing depth


##individual
sample_sum_indivhl_df <- data.frame(ASV_count = sample_sums(individual_hl))
metadata_indivhl_df <- as(sample_data(individual_hl), "data.frame" )
seqdepth_indivhl_df <- cbind(sample_sum_indivhl_df,metadata_indivhl_df)

seqdepth_indivhl_df$Number.in.Pool <- as.factor(seqdepth_indivhl_df$Number.in.Pool)

ggplot(data=seqdepth_indivhl_df, aes(x=PrevH, y=ASV_count, color=PrevH))+
  geom_boxplot()+
  stat_compare_means() #also different looking at inviduals only in one "population"


###randomvsubset
#high
sample_sum_pool3v6high_df <- data.frame(ASV_count = sample_sums(pools_3v6high))
metadata_pool3v6high_df <- as(sample_data(pools_3v6high), "data.frame" )
seqdepth_pool3v6high_df <- cbind(sample_sum_pool3v6high_df,metadata_pool3v6high_df)

seqdepth_pool3v6high_df$Number.in.Pool <- as.factor(seqdepth_pool3v6high_df$Number.in.Pool)

ggplot(data=seqdepth_pool3v6high_df, aes(x=PoolingType, y=ASV_count, color=PoolingType))+
  geom_boxplot()+
  stat_compare_means() #No difference

#low
sample_sum_pool3v6low_df <- data.frame(ASV_count = sample_sums(pools_3v6low))
metadata_pool3v6low_df <- as(sample_data(pools_3v6low), "data.frame" )
seqdepth_pool3v6low_df <- cbind(sample_sum_pool3v6low_df,metadata_pool3v6low_df)

seqdepth_pool3v6low_df$Number.in.Pool <- as.factor(seqdepth_pool3v6low_df$Number.in.Pool)

ggplot(data=seqdepth_pool3v6low_df, aes(x=PoolingType, y=ASV_count, color=PoolingType))+
  geom_boxplot()+
  facet_wrap(~Number.in.Pool)+
  stat_compare_means() #No difference

##Getting Number at each level

data_kingdom <- tax_glom(data_50k, taxrank = "Kingdom",NArm = F) #3 Kingdoms
data_kingdom_melt <- psmelt(data_kingdom)
unique(data_kingdom_melt$Kingdom) #"Bacteria"   "Unassigned" "Archaea"
write.csv(unique(data_kindom_melt$Kingdom),"kindomsindata.csv")

data_phylum <- tax_glom(data_50k, taxrank = "Phylum",NArm=F) #60 Phyla
data_phylum_melt <- psmelt(data_phylum)
unique(data_phylum_melt$Phylum)
write.csv(unique(data_phylum_melt$Phylum), "phylaindata.csv")

# [1] "Firmicutes"                   "Proteobacteria"               "Actinobacteriota"             "Bacteroidota"                 "unclassified Bacteria"        "unclassified Unassigned"     
#[7] "Deinococcota"                 "Spirochaetota"                "Fibrobacterota"               "Chloroflexi"                  "Verrucomicrobiota"            "Patescibacteria"             
#[13] "Fusobacteriota"               "Cyanobacteria"                "Planctomycetota"              "Acidobacteriota"              "Desulfobacterota"             "Gemmatimonadota"             
#[19] "Campylobacterota"             "Myxococcota"                  "Euryarchaeota"                "Cloacimonadota"               "Halobacterota"                "Synergistota"                
#[25] "Bdellovibrionota"             "Sumerlaeota"                  "Armatimonadota"               "SAR324_clade(Marine_group_B)" "Elusimicrobiota"              "Thermoplasmatota"            
#[31] "Nitrospirota"                 "Hydrogenedentes"              "Dependentiae"                 "Crenarchaeota"                "MBNT15"                       "Abditibacteriota"            
#[37] "Sva0485"                      "Entotheonellaeota"            "NB1-j"                        "Methylomirabilota"            "Latescibacterota"             "Caldatribacteriota"          
#[43] "WOR-1"                        "Aquificota"                   "RCP2-54"                      "GAL15"                        "WPS-2"                        "WS2"                         
#[49] "Dadabacteria"                 "Chrysiogenetota"              "WS4"                          "Halanaerobiaeota"             "uncultured"                   "Nanoarchaeota"               
#[55] "Fermentibacterota"            "Calditrichota"                "FW113"                        "TX1A-33"                      "Margulisbacteria"             "Modulibacteria"              

data_class <- tax_glom(data_50k, taxrank = "Class", NArm = F) #166 classes, 158 samples
data_class_melt <- psmelt(data_class) 
unique(data_class_melt$Class) #162 unique? differences from NArm?
write.csv(unique(data_class_melt$Class),"classesindata.csv") 

data_order <- tax_glom(data_50k, taxrank = "Order", NArm = F) #401 orders
data_order_melt <- psmelt(data_order)
unique(data_order_melt$Order) #389 unique? differences from NArm?
write.csv(unique(data_order_melt$Order), "ordersindata.csv")

data_family <- tax_glom(data_50k, taxrank = "Family", NArm = F) #767 families
data_family_melt <- psmelt(data_family)
unique(data_family_melt$Family) #721 unique? difference from NArm?
write.csv(unique(data_family_melt$Family), "familiesindata.csv")

data_genus <- tax_glom(data_50k, taxrank = "Genus", NArm = F) # 2024 genera
data_genus_melt <- psmelt(data_genus)
unique(data_genus_melt$Genus) #1826 unique? difference from NArm?
write.csv(unique(data_genus_melt$Genus), "generaindata.csv")

#data_species <- tax_glom(data, taxrank = "Species", NArm = F) #
#data_species_melt <- psmelt(data_species)
#unique(data_species_melt$Species)
#write.csv(unique(data_species_melt$Species), "speciesindata.csv")

###Writing csvs
write.csv(tax_table(data_kingdom),"kingdom_taxa.csv")
write.csv(otu_table(data_kingdom),"kingdom_otu.csv")

write.csv(tax_table(data_phylum),"phylum_taxa.csv")
write.csv(otu_table(data_phylum),"phylum_otu.csv")

write.csv(tax_table(data_class),"class_taxa.csv")
write.csv(otu_table(data_class),"class_otu.csv")

write.csv(tax_table(data_order),"order_taxa.csv")
write.csv(otu_table(data_order),"order_otu.csv")

write.csv(tax_table(data_family),"family_taxa.csv")
write.csv(otu_table(data_family),"family_otu.csv")

write.csv(tax_table(data_genus),"genus_taxa.csv")
write.csv(otu_table(data_genus),"genus_otu.csv")

#write.csv(tax_table(data_species),"species_taxa.csv")
#write.csv(otu_table(data_species),"species_otu.csv")

#write.csv(tax_table(data), "data_taxa.csv")
#write.csv(otu_table(data),"data_otu.csv")






