#Rarefaction for Reviewers

any(sample_sums(data_50k)==0)
data_50k


data_50k_rarefiedNR <- rarefy_even_depth(data_50k, rngseed = 888, sample.size = min(sample_sums(data_50k)), replace = F)

data_50k_rarefiedR <- rarefy_even_depth(data_50k, rngseed = 888, sample.size = min(sample_sums(data_50k)), replace = T)

any(sample_sums(data_50k_rarefiedNR)==0)

any(sample_sums(data_50k_rarefiedR)==0)

####WITH REPLACEMENT####

##Split all data by prevalence (thinking about this as just pools, not worried about how they were selected)
#high
high_all_rareR<-subset_samples(data_50k_rarefiedR, PrevH=="Y")
sum(taxa_sums(high_all_rareR)==0) #20363
high_all_rareR <-prune_taxa(taxa_sums(high_all_rareR)>0,high_all_rareR)
sum(sample_sums(high_all_rareR)==0)

#low
low_all_rareR<-subset_samples(data_50k_rarefiedR, PrevL=="Y")
sum(taxa_sums(low_all_rareR)==0) #5972
low_all_rareR<-prune_taxa(taxa_sums(low_all_rareR)>0,low_all_rareR)
sum(sample_sums(low_all_rareR)==0)


##Split up pooling strategies

#subset
subsetpooling_rareR <- subset_samples(data_50k_rarefiedR, PoolingType!="Random") #137 samples
sum(taxa_sums(subsetpooling_rareR)==0)#2413
subsetpooling_rareR <- prune_taxa(taxa_sums(subsetpooling_rareR)>0, subsetpooling_rareR)


##HP
subsethigh_rareR <- subset_samples(subsetpooling_rareR,PrevH=="Y") #88 samples
sum(taxa_sums(subsethigh_rareR)==0) #19051
subsethigh_rareR <- prune_taxa(taxa_sums(subsethigh_rareR)>0, subsethigh_rareR)
sample_data(subsethigh_rareR)$Number.in.Pool <- as.factor(sample_data(subsethigh_rareR)$Number.in.Pool)



##LP
subsetlow_rareR <- subset_samples(subsetpooling_rareR,PrevL=="Y") #88 samples
sum(taxa_sums(subsetlow_rareR)==0) #5309
subsetlow_rareR <- prune_taxa(taxa_sums(subsetlow_rareR)>0, subsetlow_rareR)

#random
randompooling_rareR <- subset_samples(data_50k_rarefiedR,PoolingType!="Subset") #135 samples
sum(taxa_sums(randompooling_rareR)==0) #3447
randompooling_rareR <- prune_taxa(taxa_sums(randompooling_rareR)>0, randompooling_rareR)

##HP
randomhigh_rareR <- subset_samples(randompooling_rareR,PrevH=="Y") #88 samples
sum(taxa_sums(randomhigh_rareR)==0) #18480
randomhigh_rareR <- prune_taxa(taxa_sums(randomhigh_rareR)>0, randomhigh_rareR)

##LP
randomlow_rareR <- subset_samples(randompooling_rareR,PrevL=="Y") #88 samples
sum(taxa_sums(randomlow_rareR)==0) #5061
randomlow_rareR <- prune_taxa(taxa_sums(randomlow_rareR)>0, randomlow_rareR)

#both for comparison
pools_3v6_rareR <- subset_samples(data_50k_rarefiedR, PoolingType=="Subset"|PoolingType=="Random") #44 samples
pools_3v6_rareR <- prune_taxa(taxa_sums(pools_3v6_rareR)>0, pools_3v6_rareR)

##HP
pools_3v6high <- subset_samples(pools_3v6,PrevH=="Y") #20 samples
pools_3v6high <- prune_taxa(taxa_sums(pools_3v6high)>0,pools_3v6high)

##LP
pools_3v6low <- subset_samples(pools_3v6,PrevL=="Y") #24 samples
pools_3v6low <- prune_taxa(taxa_sums(pools_3v6low)>0,pools_3v6low)





##Just compare high and low individuals
View(sample_data(data_50k_rarefiedR))

##high##
high_individual_rareR <- subset_samples(data_50k_rarefiedR, PrevH=="Y"&PrevL=="N"&Number.in.Pool=="1") #32, 10964 taxa
sum(taxa_sums(high_individual_rareR)==0) #28889


##low##
low_individual_rareR <- subset_samples(data_50k_rarefiedR, PrevH=="N"&PrevL=="Y"&Number.in.Pool=="1")
sum(taxa_sums(low_individual_rareR)==0)#21541

individual_hl_rareR<- merge_phyloseq(high_individual_rareR,low_individual_rareR)
sum(taxa_sums(individual_hl_rareR)==0)#16425
individual_hl_rareR <- prune_taxa(taxa_sums(individual_hl_rareR)>0, individual_hl_rareR) #63 samples, 31415 taxa

high_individual_rareR<- prune_taxa(taxa_sums(high_individual_rareR)>0,high_individual_rareR)
low_individual_rareR<- prune_taxa(taxa_sums(low_individual_rareR)>0,low_individual_rareR) #31 sample, 25575 taxa






poolcomparisons = list(c("1","3"),c("1","6"),c("1","12"),c("3","6"),c("3","12"),c("6","12"))

#For Mannheimia culture status
#all
sample_sum_df_50k_rareR <- data.frame(ASVs = sample_sums(data_50k_rarefiedR))
readplot_all_50k_rareR <- ggplot(sample_sum_df_50k_rareR, aes(x = ASVs)) + 
  geom_histogram(color = "black", fill="skyblue", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) #  Talk to Lee about distribution
readplot_all_50k_rareR

any(taxa_sums(data_50k_rarefiedR)==0)

data_50k_metadata_rareR <- as(sample_data(data_50k_rarefiedR),"data.frame")
seqdepth_50k_rareR <- cbind(sample_sum_df_50k_rareR,data_50k_metadata_rareR)
View(seqdepth_50k_rareR)

seqdepth_50k$Number.in.Pool <- as.factor(seqdepth_50k$Number.in.Pool)

kruskal_test(seqdepth_50k_rareR, ASVs~Number.in.Pool)

ggplot(data=seqdepth_50k_rareR, aes(x=Number.in.Pool, y=ASVs, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  facet_wrap(~Contains.CP)+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)

####Without replacement####
#high
high_all_rareNR<-subset_samples(data_50k_rarefiedNR, PrevH=="Y")
sum(taxa_sums(high_all_rareNR)==0) #20821
high_all_rareNR <-prune_taxa(taxa_sums(high_all_rareNR)>0,high_all_rareNR)
sum(sample_sums(high_all_rareNR)==0)

#low
low_all_rareNR<-subset_samples(data_50k_rarefiedNR, PrevL=="Y")
sum(taxa_sums(low_all_rareNR)==0) #6125
low_all_rareNR<-prune_taxa(taxa_sums(low_all_rareNR)>0,low_all_rareNR)
sum(sample_sums(low_all_rareNR)==0)


##Split up pooling strategies

#subset
subsetpooling_rareNR <- subset_samples(data_50k_rarefiedNR, PoolingType!="Random") #137 samples
sum(taxa_sums(subsetpooling_rareNR)==0)#2483
subsetpooling_rareNR <- prune_taxa(taxa_sums(subsetpooling_rareNR)>0, subsetpooling_rareNR)


##HP
subsethigh_rareNR <- subset_samples(subsetpooling_rareNR,PrevH=="Y") #88 samples
sum(taxa_sums(subsethigh_rareNR)==0) #19441
subsethigh_rareNR <- prune_taxa(taxa_sums(subsethigh_rareNR)>0, subsethigh_rareNR)
sample_data(subsethigh_rareNR)$Number.in.Pool <- as.factor(sample_data(subsethigh_rareNR)$Number.in.Pool)



##LP
subsetlow_rareNR <- subset_samples(subsetpooling_rareNR,PrevL=="Y") #88 samples
sum(taxa_sums(subsetlow_rareNR)==0) #5470
subsetlow_rareNR <- prune_taxa(taxa_sums(subsetlow_rareNR)>0, subsetlow_rareNR)

#random
randompooling_rareNR <- subset_samples(data_50k_rarefiedNR,PoolingType!="Subset") #135 samples
sum(taxa_sums(randompooling_rareNR)==0) #3508
randompooling_rareNR <- prune_taxa(taxa_sums(randompooling_rareNR)>0, randompooling_rareNR)

##HP
randomhigh_rareNR <- subset_samples(randompooling_rareNR,PrevH=="Y") #88 samples
sum(taxa_sums(randomhigh_rareNR)==0) #19026
randomhigh_rareNR <- prune_taxa(taxa_sums(randomhigh_rareNR)>0, randomhigh_rareNR)

##LP
randomlow_rareNR <- subset_samples(randompooling_rareNR,PrevL=="Y") #88 samples
sum(taxa_sums(randomlow_rareNR)==0) #5154
randomlow_rareNR <- prune_taxa(taxa_sums(randomlow_rareNR)>0, randomlow_rareNR)

#both for comparison
pools_3v6_rareNR <- subset_samples(data_50k_rarefiedNR, PoolingType=="Subset"|PoolingType=="Random") #44 samples
pools_3v6_rareNR <- prune_taxa(taxa_sums(pools_3v6_rareNR)>0, pools_3v6_rareNR)

##HP
pools_3v6high <- subset_samples(pools_3v6,PrevH=="Y") #20 samples
pools_3v6high <- prune_taxa(taxa_sums(pools_3v6high)>0,pools_3v6high)

##LP
pools_3v6low <- subset_samples(pools_3v6,PrevL=="Y") #24 samples
pools_3v6low <- prune_taxa(taxa_sums(pools_3v6low)>0,pools_3v6low)


##Just compare high and low individuals
View(sample_data(data_50k_rarefiedNR))

##high##
high_individual_rareNR <- subset_samples(data_50k_rarefiedNR, PrevH=="Y"&PrevL=="N"&Number.in.Pool=="1") #32, 10964 taxa
sum(taxa_sums(high_individual_rareNR)==0) #29597


##low##
low_individual_rareNR <- subset_samples(data_50k_rarefiedNR, PrevH=="N"&PrevL=="Y"&Number.in.Pool=="1")
sum(taxa_sums(low_individual_rareNR)==0)#22031

individual_hl_rareNR<- merge_phyloseq(high_individual_rareNR,low_individual_rareNR)
sum(taxa_sums(individual_hl_rareNR)==0)#16425
individual_hl_rareNR <- prune_taxa(taxa_sums(individual_hl_rareNR)>0, individual_hl_rareNR) #63 samples, 31415 taxa

high_individual_rareNR<- prune_taxa(taxa_sums(high_individual_rareNR)>0,high_individual_rareNR)
low_individual_rareNR<- prune_taxa(taxa_sums(low_individual_rareNR)>0,low_individual_rareNR) #31 sample, 25575 taxa


#For Mannheimia culture status
#all
sample_sum_df_50k_rareNR <- data.frame(ASVs = sample_sums(data_50k_rarefiedNR))
readplot_all_50k_rareNR <- ggplot(sample_sum_df_50k_rareNR, aes(x = ASVs)) + 
  geom_histogram(color = "black", fill="skyblue", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) #  Talk to Lee about distribution
readplot_all_50k_rareNR

any(taxa_sums(data_50k_rarefiedNR)==0)

data_50k_metadata_rareNR <- as(sample_data(data_50k_rarefiedNR),"data.frame")
seqdepth_50k_rareNR <- cbind(sample_sum_df_50k_rareNR,data_50k_metadata_rareNR)
View(seqdepth_50k_rareNR)

seqdepth_50k$Number.in.Pool <- as.factor(seqdepth_50k$Number.in.Pool)

kruskal_test(seqdepth_50k_rareNR, ASVs~Number.in.Pool)

ggplot(data=seqdepth_50k_rareNR, aes(x=Number.in.Pool, y=ASVs, color=Number.in.Pool, fill=Number.in.Pool, alpha=Number.in.Pool))+
  geom_boxplot()+
  facet_wrap(~Contains.CP)+
  scale_fill_manual(values = numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5)) +
  scale_color_manual(values = numberinpoolpalette)


####DON'T NEED THIS SECTION####
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

####Getting Number at each level####

data_kingdom_rareR <- tax_glom(data_50k_rarefiedR, taxrank = "Kingdom",NArm = F) #3 Kingdoms
data_kingdom_melt_rareR <- psmelt(data_kingdom_rareR)
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

####ALPHA DIVERSITY-With Replacement####

##double check no samples with 0 reads
any(sample_sums(data_50k_rarefiedR)==0)
any(sample_sums(subsethigh_rareR)==0)
any(sample_sums(subsetlow_rareR)==0)
any(sample_sums(randomhigh_rareR)==0)
any(sample_sums(randomlow_rareR)==0)
any(sample_sums(pools_3v6high)==0)
any(sample_sums(pools_3v6low)==0)
any(sample_sums(individual_hl_rareR)==0)
any(sample_sums(high_all_rareR)==0)
any(sample_sums(low_all_rareR)==0)


##None
##double check taxa sums
any(taxa_sums(data_50k_rarefiedR)==0)
any(taxa_sums(subsethigh_rareR)==0)
any(taxa_sums(subsetlow_rareR)==0)
any(taxa_sums(randomhigh_rareR)==0)
any(taxa_sums(randomlow_rareR)==0)
any(taxa_sums(pools_3v6high)==0)
any(taxa_sums(pools_3v6low)==0)
any(taxa_sums(individual_hl_rareR)==0)
any(taxa_sums(high_all_rareR)==0)
any(taxa_sums(low_all_rareR)==0)

#also false


##make a palette

numberinpoolpalette <- viridis_pal( option = "H")(4)

#make comparison list
poolonly_comparisons <- list(c("3","6"),c("3","12"),c("6","12"))

##all high_pools
high_alpha_div1_rareR <- estimate_richness(high_all_rareR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
high_alpha_div2_rareR <- estimate_pd(high_all_rareR) # calculating Faith's PD

high_alpha_div_rareR <- cbind(high_alpha_div1_rareR, high_alpha_div2_rareR) # combining Faith's and other metrics
high_alpha_div_rareR # looks good, but have some duplicates so lets trim it a bit
high_alpha_div_rareR <- high_alpha_div_rareR[,c(1:5)]
high_alpha_div_rareR
high_alpha_div_rareR.df <- as(sample_data(high_all_rareR), "data.frame") # making into DF for metadata
high_alpha_div_rareR_meta <- cbind(high_alpha_div_rareR, high_alpha_div_rareR.df)
high_alpha_div_rareR_meta 

ggplot(high_alpha_div_rareR_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(high_alpha_div_rareR_meta, Observed~Number.in.Pool) #P=0.0623
dunn_test(high_alpha_div_rareR_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #NS,1-12 P.adj=0.261

ggplot(high_alpha_div_rareR_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's Diversity", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(high_alpha_div_rareR_meta, Shannon~Number.in.Pool) #P=0.0092
dunn_test(high_alpha_div_rareR_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-12 P.adj=0.0432


##all low_pools
low_alpha_div1_rareR <- estimate_richness(low_all_rareR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
low_alpha_div2_rareR <- estimate_pd(low_all_rareR) # calculating Faith's PD

low_alpha_div_rareR <- cbind(low_alpha_div1_rareR, low_alpha_div2_rareR) # combining Faith's and other metrics
low_alpha_div_rareR # looks good, but have some duplicates so lets trim it a bit
low_alpha_div_rareR <- low_alpha_div_rareR[,c(1:5)]
low_alpha_div_rareR
low_alpha_div_rareR.df <- as(sample_data(low_all_rareR), "data.frame") # making into DF for metadata
low_alpha_div_rareR_meta <- cbind(low_alpha_div_rareR, low_alpha_div_rareR.df)
low_alpha_div_rareR_meta 


ggplot(low_alpha_div_rareR_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(low_alpha_div_rareR_meta, Observed~Number.in.Pool) #P=0.000159
dunn_test(low_alpha_div_rareR_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #1-6, P=0.0102; 1-12 P=0.00225

ggplot(low_alpha_div_rareR_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(low_alpha_div_rareR_meta, Shannon~Number.in.Pool) #NS, 0.000666
dunn_test(low_alpha_div_rareR_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-6, P=0.0174; 1-12, P=0.00678

##subset_highprev
subseth_rareR_alpha_div1 <- estimate_richness(subsethigh_rareR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
subseth_rareR_alpha_div2 <- estimate_pd(subsethigh_rareR) # calculating Faith's PD

subseth_rareR_alpha_div <- cbind(subseth_rareR_alpha_div1, subseth_rareR_alpha_div2) # combining Faith's and other metrics
subseth_rareR_alpha_div # looks good, but have some duplicates so lets trim it a bit
subseth_rareR_alpha_div <- subseth_rareR_alpha_div[,c(1:5)]
subseth_rareR_alpha_div
subseth_rareR_alpha_div.df <- as(sample_data(subsethigh_rareR), "data.frame") # making into DF for metadata
str(subseth_rareR_alpha_div.df)
subseth_rareR_alpha_div_meta <- cbind(subseth_rareR_alpha_div, subseth_rareR_alpha_div.df)
subseth_rareR_alpha_div_meta # great now we've got alpha_div values and metadata
subseth_rareR_alpha_div_meta$Number.in.Pool <- as.factor(subseth_rareR_alpha_div_meta$Number.in.Pool)


hs_rareR_observed <- ggplot(subseth_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))
hs_rareR_observed

kruskal_test(subseth_rareR_alpha_div_meta, Observed~Number.in.Pool)#P=0.0322
dunn_test(subseth_rareR_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #NS, 1-3 and 1-12 p-adj=0.128



hs_rareR_shannon <- ggplot(subseth_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's Diversity Index", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))#ND

hs_rareR_shannon

kruskal_test(subseth_rareR_alpha_div_meta, Shannon~Number.in.Pool) #P=0.0152
dunn_test(subseth_rareR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #NS, 1-12 P.adj.=0.0563; 1-3 =0.137



##subset_lowprev
subsetl_rareR_alpha_div1 <- estimate_richness(subsetlow_rareR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
subsetl_rareR_alpha_div2 <- estimate_pd(subsetlow_rareR) # calculating Faith's PD

subsetl_rareR_alpha_div <- cbind(subsetl_rareR_alpha_div1, subsetl_rareR_alpha_div2) # combining Faith's and other metrics
subsetl_rareR_alpha_div # looks good, but have some duplicates so lets trim it a bit
subsetl_rareR_alpha_div <- subsetl_rareR_alpha_div[,c(1:5)]
subsetl_rareR_alpha_div
subsetl_rareR_alpha_div.df <- as(sample_data(subsetlow_rareR), "data.frame") # making into DF for metadata
str(subsetl_rareR_alpha_div.df)
subsetl_rareR_alpha_div_meta <- cbind(subsetl_rareR_alpha_div, subsetl_rareR_alpha_div.df)
subsetl_rareR_alpha_div_meta # great now we've got alpha_div values and metadata
subsetl_rareR_alpha_div_meta$Number.in.Pool <- as.factor(subsetl_rareR_alpha_div_meta$Number.in.Pool)

ls_observed_rareR <- ggplot(subsetl_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))
ls_observed_rareR

kruskal_test(subsetl_rareR_alpha_div_meta, Observed~Number.in.Pool) #NS, P=0.000782
dunn_test(subsetl_rareR_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #1-12, p,adj=0.00361

ls_shannon_rareR <- ggplot(subsetl_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon,fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(y= "Shannon's Diversity Index", x= "Number in Pool", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))

kruskal_test(subsetl_rareR_alpha_div_meta, Shannon~Number.in.Pool) #p=0.00317
dunn_test(subsetl_rareR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-12, P.adj.=0.0122

###Random
##highrand
randomh_rareR_alpha_div1 <- estimate_richness(randomhigh_rareR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
randomh_rareR_alpha_div2 <- estimate_pd(randomhigh_rareR) # calculating Faith's PD

randomh_rareR_alpha_div <- cbind(randomh_rareR_alpha_div1, randomh_rareR_alpha_div2) # combining Faith's and other metrics
randomh_rareR_alpha_div # looks good, but have some duplicates so lets trim it a bit
randomh_rareR_alpha_div <- randomh_rareR_alpha_div[,c(1:5)]
randomh_rareR_alpha_div
randomh_rareR_alpha_div.df <- as(sample_data(randomhigh_rareR), "data.frame") # making into DF for metadata
str(randomh_rareR_alpha_div.df)
randomh_rareR_alpha_div_meta <- cbind(randomh_rareR_alpha_div, randomh_rareR_alpha_div.df)
randomh_rareR_alpha_div_meta # great now we've got alpha_div values and metadata
randomh_rareR_alpha_div_meta$Number.in.Pool <- as.factor(randomh_rareR_alpha_div_meta$Number.in.Pool)


hr_rareR_observed <- ggplot(randomh_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))
hr_rareR_observed

kruskal_test(randomh_rareR_alpha_div_meta, Observed~Number.in.Pool)#P=0.146, NS

hr_rareR_shannon <- ggplot(randomh_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's Diversity Index", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))#ND

hr_rareR_shannon

kruskal_test(randomh_rareR_alpha_div_meta, Shannon~Number.in.Pool) #P=0.0.0375
dunn_test(randomh_rareR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #NS, 1-12 P.adj.=0.0529


##random_lowprev
randoml_rareR_alpha_div1 <- estimate_richness(randomlow_rareR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
randoml_rareR_alpha_div2 <- estimate_pd(randomlow_rareR) # calculating Faith's PD

randoml_rareR_alpha_div <- cbind(randoml_rareR_alpha_div1, randoml_rareR_alpha_div2) # combining Faith's and other metrics
randoml_rareR_alpha_div # looks good, but have some duplicates so lets trim it a bit
randoml_rareR_alpha_div <- randoml_rareR_alpha_div[,c(1:5)]
randoml_rareR_alpha_div
randoml_rareR_alpha_div.df <- as(sample_data(randomlow_rareR), "data.frame") # making into DF for metadata
str(randoml_rareR_alpha_div.df)
randoml_rareR_alpha_div_meta <- cbind(randoml_rareR_alpha_div, randoml_rareR_alpha_div.df)
randoml_rareR_alpha_div_meta # great now we've got alpha_div values and metadata
randoml_rareR_alpha_div_meta$Number.in.Pool <- as.factor(randoml_rareR_alpha_div_meta$Number.in.Pool)

lr_observed_rareR <- ggplot(randoml_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))
lr_observed_rareR

kruskal_test(randoml_rareR_alpha_div_meta, Observed~Number.in.Pool) #NS, P=0.00142
dunn_test(randoml_rareR_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #1-12, p,adj=0.00408

lr_shannon_rareR <- ggplot(randoml_rareR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon,fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(y= "Shannon's Diversity Index", x= "Number in Pool", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))
lr_shannon_rareR

kruskal_test(randoml_rareR_alpha_div_meta, Shannon~Number.in.Pool) #p=0.00364
dunn_test(randoml_rareR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-12, P.adj.=0.0950

ggarrange(hs_rareR_observed,ls_observed_rareR,hr_rareR_observed,lr_observed_rareR, hs_rareR_shannon, ls_shannon_rareR, hr_rareR_shannon, lr_shannon_rareR, nrow = 2, ncol = 4, common.legend = F)

####NoReplacemnt####
##double check no samples with 0 reads
any(sample_sums(data_50k_rarefiedR)==0)
any(sample_sums(subsethigh_rareR)==0)
any(sample_sums(subsetlow_rareR)==0)
any(sample_sums(randomhigh_rareR)==0)
any(sample_sums(randomlow_rareR)==0)
any(sample_sums(pools_3v6high)==0)
any(sample_sums(pools_3v6low)==0)
any(sample_sums(individual_hl_rareR)==0)
any(sample_sums(high_all_rareR)==0)
any(sample_sums(low_all_rareR)==0)


##None
##double check taxa sums
any(taxa_sums(data_50k_rarefiedR)==0)
any(taxa_sums(subsethigh_rareR)==0)
any(taxa_sums(subsetlow_rareR)==0)
any(taxa_sums(randomhigh_rareR)==0)
any(taxa_sums(randomlow_rareR)==0)
any(taxa_sums(pools_3v6high)==0)
any(taxa_sums(pools_3v6low)==0)
any(taxa_sums(individual_hl_rareR)==0)
any(taxa_sums(high_all_rareR)==0)
any(taxa_sums(low_all_rareR)==0)

#also false


##make a palette

numberinpoolpalette <- viridis_pal( option = "H")(4)

#make comparison list
poolonly_comparisons <- list(c("3","6"),c("3","12"),c("6","12"))

##all high_pools
high_alpha_div1_rareNR <- estimate_richness(high_all_rareNR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
high_alpha_div2_rareNR <- estimate_pd(high_all_rareNR) # calculating Faith's PD

high_alpha_div_rareNR <- cbind(high_alpha_div1_rareNR, high_alpha_div2_rareNR) # combining Faith's and other metrics
high_alpha_div_rareNR # looks good, but have some duplicates so lets trim it a bit
high_alpha_div_rareNR <- high_alpha_div_rareNR[,c(1:5)]
high_alpha_div_rareNR
high_alpha_div_rareNR.df <- as(sample_data(high_all_rareNR), "data.frame") # making into DF for metadata
high_alpha_div_rareNR_meta <- cbind(high_alpha_div_rareNR, high_alpha_div_rareNR.df)
high_alpha_div_rareNR_meta 


ggplot(high_alpha_div_rareNR_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(high_alpha_div_rareNR_meta, Observed~Number.in.Pool) #P=0.0704
dunn_test(high_alpha_div_rareNR_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #NS,1-12 P.adj=0.266

ggplot(high_alpha_div_rareNR_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's Diversity", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(high_alpha_div_rareNR_meta, Shannon~Number.in.Pool) #P=0.00986
dunn_test(high_alpha_div_rareNR_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-12 P.adj=0.0473


##all low_pools
low_alpha_div1_rareNR <- estimate_richness(low_all_rareNR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
low_alpha_div2_rareNR <- estimate_pd(low_all_rareNR) # calculating Faith's PD

low_alpha_div_rareNR <- cbind(low_alpha_div1_rareNR, low_alpha_div2_rareNR) # combining Faith's and other metrics
low_alpha_div_rareNR # looks good, but have some duplicates so lets trim it a bit
low_alpha_div_rareNR <- low_alpha_div_rareNR[,c(1:5)]
low_alpha_div_rareNR
low_alpha_div_rareNR.df <- as(sample_data(low_all_rareNR), "data.frame") # making into DF for metadata
low_alpha_div_rareNR_meta <- cbind(low_alpha_div_rareNR, low_alpha_div_rareNR.df)
low_alpha_div_rareNR_meta 


ggplot(low_alpha_div_rareNR_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(low_alpha_div_rareNR_meta, Observed~Number.in.Pool) #P=0.000153
dunn_test(low_alpha_div_rareNR_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #1-6, P=0.0106; 1-12 P=0.00200

ggplot(low_alpha_div_rareNR_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))

kruskal_test(low_alpha_div_rareNR_meta, Shannon~Number.in.Pool) #NS, 0.000705
dunn_test(low_alpha_div_rareNR_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-6, P=0.0185; 1-12, P=0.00684

##subset_highprev
subseth_rareNR_alpha_div1 <- estimate_richness(subsethigh_rareNR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
subseth_rareNR_alpha_div2 <- estimate_pd(subsethigh_rareNR) # calculating Faith's PD

subseth_rareNR_alpha_div <- cbind(subseth_rareNR_alpha_div1, subseth_rareNR_alpha_div2) # combining Faith's and other metrics
subseth_rareNR_alpha_div # looks good, but have some duplicates so lets trim it a bit
subseth_rareNR_alpha_div <- subseth_rareNR_alpha_div[,c(1:5)]
subseth_rareNR_alpha_div
subseth_rareNR_alpha_div.df <- as(sample_data(subsethigh_rareNR), "data.frame") # making into DF for metadata
str(subseth_rareNR_alpha_div.df)
subseth_rareNR_alpha_div_meta <- cbind(subseth_rareNR_alpha_div, subseth_rareNR_alpha_div.df)
subseth_rareNR_alpha_div_meta # great now we've got alpha_div values and metadata
subseth_rareNR_alpha_div_meta$Number.in.Pool <- as.factor(subseth_rareNR_alpha_div_meta$Number.in.Pool)


hs_rareNR_observed <- ggplot(subseth_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=10,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))
hs_rareNR_observed

kruskal_test(subseth_rareNR_alpha_div_meta, Observed~Number.in.Pool)#P=0.0378
dunn_test(subseth_rareNR_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #NS, 1-3 and 1-12 p-adj=0.134



hs_rareNR_shannon <- ggplot(subseth_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's Diversity Index", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))#ND

hs_rareNR_shannon

kruskal_test(subseth_rareNR_alpha_div_meta, Shannon~Number.in.Pool) #P=0.0159
dunn_test(subseth_rareNR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #NS, 1-12 P.adj.=0.0591; 1-3 =0.137



##subset_lowprev
subsetl_rareNR_alpha_div1 <- estimate_richness(subsetlow_rareNR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
subsetl_rareNR_alpha_div2 <- estimate_pd(subsetlow_rareNR) # calculating Faith's PD

subsetl_rareNR_alpha_div <- cbind(subsetl_rareNR_alpha_div1, subsetl_rareNR_alpha_div2) # combining Faith's and other metrics
subsetl_rareNR_alpha_div # looks good, but have some duplicates so lets trim it a bit
subsetl_rareNR_alpha_div <- subsetl_rareNR_alpha_div[,c(1:5)]
subsetl_rareNR_alpha_div
subsetl_rareNR_alpha_div.df <- as(sample_data(subsetlow_rareNR), "data.frame") # making into DF for metadata
str(subsetl_rareNR_alpha_div.df)
subsetl_rareNR_alpha_div_meta <- cbind(subsetl_rareNR_alpha_div, subsetl_rareNR_alpha_div.df)
subsetl_rareNR_alpha_div_meta # great now we've got alpha_div values and metadata
subsetl_rareNR_alpha_div_meta$Number.in.Pool <- as.factor(subsetl_rareNR_alpha_div_meta$Number.in.Pool)

ls_observed_rareNR <- ggplot(subsetl_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=10,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))
ls_observed_rareNR

kruskal_test(subsetl_rareNR_alpha_div_meta, Observed~Number.in.Pool) # P=0.00078
dunn_test(subsetl_rareNR_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #1-12, p,adj=0.00317

ls_shannon_rareNR <- ggplot(subsetl_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon,fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(y= "Shannon's Diversity Index", x= "Number in Pool", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))

kruskal_test(subsetl_rareNR_alpha_div_meta, Shannon~Number.in.Pool) #p=0.00315
dunn_test(subsetl_rareNR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-12, P.adj.=0.0122

###Random
##highrand
randomh_rareNR_alpha_div1 <- estimate_richness(randomhigh_rareNR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
randomh_rareNR_alpha_div2 <- estimate_pd(randomhigh_rareNR) # calculating Faith's PD

randomh_rareNR_alpha_div <- cbind(randomh_rareNR_alpha_div1, randomh_rareNR_alpha_div2) # combining Faith's and other metrics
randomh_rareNR_alpha_div # looks good, but have some duplicates so lets trim it a bit
randomh_rareNR_alpha_div <- randomh_rareNR_alpha_div[,c(1:5)]
randomh_rareNR_alpha_div
randomh_rareNR_alpha_div.df <- as(sample_data(randomhigh_rareNR), "data.frame") # making into DF for metadata
str(randomh_rareNR_alpha_div.df)
randomh_rareNR_alpha_div_meta <- cbind(randomh_rareNR_alpha_div, randomh_rareNR_alpha_div.df)
randomh_rareNR_alpha_div_meta # great now we've got alpha_div values and metadata
randomh_rareNR_alpha_div_meta$Number.in.Pool <- as.factor(randomh_rareNR_alpha_div_meta$Number.in.Pool)


hr_rareNR_observed <- ggplot(randomh_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=10,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))
hr_rareNR_observed

kruskal_test(randomh_rareNR_alpha_div_meta, Observed~Number.in.Pool)#P=0.153, NS

hr_rareNR_shannon <- ggplot(randomh_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon, fill = Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Shannon's Diversity Index", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5))#ND

hr_rareNR_shannon

kruskal_test(randomh_rareNR_alpha_div_meta, Shannon~Number.in.Pool) #P=0.0.0408
dunn_test(randomh_rareNR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #NS, 1-12 P.adj.=0.0584


##random_lowprev
randoml_rareNR_alpha_div1 <- estimate_richness(randomlow_rareNR, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
randoml_rareNR_alpha_div2 <- estimate_pd(randomlow_rareNR) # calculating Faith's PD

randoml_rareNR_alpha_div <- cbind(randoml_rareNR_alpha_div1, randoml_rareNR_alpha_div2) # combining Faith's and other metrics
randoml_rareNR_alpha_div # looks good, but have some duplicates so lets trim it a bit
randoml_rareNR_alpha_div <- randoml_rareNR_alpha_div[,c(1:5)]
randoml_rareNR_alpha_div
randoml_rareNR_alpha_div.df <- as(sample_data(randomlow_rareNR), "data.frame") # making into DF for metadata
str(randoml_rareNR_alpha_div.df)
randoml_rareNR_alpha_div_meta <- cbind(randoml_rareNR_alpha_div, randoml_rareNR_alpha_div.df)
randoml_rareNR_alpha_div_meta # great now we've got alpha_div values and metadata
randoml_rareNR_alpha_div_meta$Number.in.Pool <- as.factor(randoml_rareNR_alpha_div_meta$Number.in.Pool)

lr_observed_rareNR <- ggplot(randoml_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(x="Number in Pool", y="Richness")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=10,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))
lr_observed_rareNR

kruskal_test(randoml_rareNR_alpha_div_meta, Observed~Number.in.Pool) #NS, P=0.00124
dunn_test(randoml_rareNR_alpha_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #1-12, p,adj=0.00350

lr_shannon_rareNR <- ggplot(randoml_rareNR_alpha_div_meta, aes(x= Number.in.Pool, y= Shannon,fill = Number.in.Pool, color=Number.in.Pool ,alpha=Number.in.Pool)) +
  theme_bw() + 
  labs(y= "Shannon's Diversity Index", x= "Number in Pool", fill="Pool Size", color="Pool Size", alpha="Pool Size")+
  geom_boxplot()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title= element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    axis.text.y=element_text(size=20,face = "bold"))+
  scale_color_manual(values=numberinpoolpalette)+
  scale_fill_manual(values=numberinpoolpalette)+
  scale_alpha_manual(values=c(0.5,0.5, 0.5, 0.5))
lr_shannon_rareNR

kruskal_test(randoml_rareNR_alpha_div_meta, Shannon~Number.in.Pool) #p=0.00398
dunn_test(randoml_rareNR_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #1-12, P.adj.=0.0959

ggarrange(hs_rareNR_observed,ls_observed_rareNR,hr_rareNR_observed,lr_observed_rareNR, hs_rareNR_shannon, ls_shannon_rareNR, hr_rareNR_shannon, lr_shannon_rareNR, nrow = 2, ncol = 4, common.legend = F)





####BETA DIVERSITY-Without Replacement####
#overall
data_16S_rareR.css <- phyloseq_transform_css(data_50k_rarefiedR, log = F)
data_16S_rareR.css
any(taxa_sums(data_16S_rareR.css)==0)
data_16S_rareR.css.ra <- transform_sample_counts(data_16S_rareR.css, function(x) {x/sum(x)} * 100)

data_16S_rareR.css.df <- as(sample_data(data_16S_rareR.css), "data.frame")

#all pools
##high
high_rareR.css <- phyloseq_transform_css(high_all_rareR,log = F)
high_rareR.css
high_rareR.css.df <- as(sample_data(high_rareR.css),"data.frame")

##low
low_rareR.css <- phyloseq_transform_css(low_all_rareR,log = F)
low_rareR.css
low_rareR.css.df <- as(sample_data(low_rareR.css),"data.frame")



#subset
subset_rareNR.css  <-phyloseq_transform_css(subsetpooling)

sh_16S_rareNR.css <- phyloseq_transform_css(subsethigh_rareNR, log = F)
sh_16S_rareNR.css.ra <- transform_sample_counts(sh_16S_rareNR.css, function(x) {x/sum(x)}*100)
sh_16S_rareNR.css.df <- as(sample_data(sh_16S_rareNR.css),"data.frame")

sl_16S_rareNR.css <- phyloseq_transform_css(subsetlow_rareNR, log = F)
sl_16S_rareNR.css.ra <- transform_sample_counts(sl_16S_rareNR.css, function(x) {x/sum(x)}*100)
sl_16s_rareNR.css.df <- as(sample_data(sl_16S_rareNR.css),"data.frame")


#random
##high
HR_16S_rareNR.css <- phyloseq_transform_css(randomhigh_rareNR,log = F)
HR_16S_rareNR.css.ra <-transform_sample_counts(HR_16S_rareNR.css, function(x) {x/sum(x)}*100)
HR_16S_rareNR.css.df <- as(sample_data(HR_16S_rareNR.css), "data.frame")

LR_16S_rareNR.css <- phyloseq_transform_css(randomlow_rareNR,log = F)
LR_16S_rareNR.css.ra <-transform_sample_counts(LR_16S_rareNR.css, function(x) {x/sum(x)}*100)
LR_16S_rareNR.css.df <- as(sample_data(LR_16S_rareNR.css), "data.frame")


#subsetpooling
#high
subseth_rareNR_wunifrac.dist <- wunifrac(sh_16S_rareNR.css)
subseth_rareNR_gunifrac.dist <- gunifrac(sh_16S_rareNR.css)
subseth_rareNR_uwunifrac.dist <- uwunifrac(sh_16S_rareNR.css)

#low
subsetl_rareNR_wunifrac.dist <- wunifrac(sl_16S_rareNR.css)
subsetl_rareNR_gunifrac.dist <- gunifrac(sl_16S_rareNR.css)
subsetl_rareNR_uwunifrac.dist <- uwunifrac(sl_16S_rareNR.css)

#randompooling
#high
randomh_rareNR_wunifrac.dist <- wunifrac(HR_16S_rareNR.css)
randomh_rareNR_gunifrac.dist <- gunifrac(HR_16S_rareNR.css)
randomh_rareNR_uwunifrac.dist <- uwunifrac(HR_16S_rareNR.css)

#low
randoml_rareNR_wunifrac.dist <- wunifrac(LR_16S_rareNR.css)
randoml_rareNR_gunifrac.dist <- gunifrac(LR_16S_rareNR.css)
randoml_rareNR_uwunifrac.dist <- uwunifrac(LR_16S_rareNR.css)




# ordinate
#subsetpooling
#high
subseth_wunifrac_rareNR.ord <- ordinate(sh_16S_rareNR.css,method = "NMDS",distance = subseth_rareNR_wunifrac.dist)
subseth_gunifrac_rareNR.ord <- ordinate(sh_16S_rareNR.css,method = "NMDS",distance = subseth_rareNR_gunifrac.dist)
subseth_uwunifrac_rareNR.ord <- ordinate(sh_16S_rareNR.css,method = "NMDS",distance = subseth_rareNR_uwunifrac.dist)

#low
subsetl_wunifrac_rareNR.ord <- ordinate(sl_16S_rareNR.css,method = "NMDS",distance = subsetl_rareNR_wunifrac.dist)
subsetl_gunifrac_rareNR.ord <- ordinate(sl_16S_rareNR.css,method = "NMDS",distance = subsetl_rareNR_gunifrac.dist)
subsetl_uwunifrac_rareNR.ord <- ordinate(sl_16S_rareNR.css,method = "NMDS",distance = subsetl_rareNR_uwunifrac.dist)

#randompooling
#high
randomh_wunifrac_rareNR.ord <- ordinate(HR_16S_rareNR.css,method = "NMDS",distance = randomh_rareNR_wunifrac.dist)
randomh_gunifrac_rareNR.ord <- ordinate(HR_16S_rareNR.css,method = "NMDS",distance = randomh_rareNR_gunifrac.dist)
randomh_uwunifrac_rareNR.ord <- ordinate(HR_16S_rareNR.css,method = "NMDS",distance = randomh_rareNR_uwunifrac.dist)

#low
randoml_wunifrac_rareNR.ord <- ordinate(LR_16S_rareNR.css,method = "NMDS",distance = randoml_rareNR_wunifrac.dist)
randoml_gunifrac_rareNR.ord <- ordinate(LR_16S_rareNR.css,method = "NMDS",distance = randoml_rareNR_gunifrac.dist)
randoml_uwunifrac_rareNR.ord <- ordinate(LR_16S_rareNR.css,method = "NMDS",distance = randoml_rareNR_uwunifrac.dist)


#PCoA-wunifrac
subseth_wunifrac_rareNR.pcoa <- ordinate(sh_16S_rareNR.css,method = "PCoA",distance = subseth_rareNR_wunifrac.dist)
subsetl_wunifrac_rareNR.pcoa <- ordinate(sl_16S_rareNR.css,method = "PCoA",distance = subsetl_rareNR_wunifrac.dist)
randomh_wunifrac_rareNR.pcoa <- ordinate(HR_16S_rareNR.css,method = "PCoA",distance = randomh_rareNR_wunifrac.dist)
randoml_wunifrac_rareNR.pcoa <- ordinate(LR_16S_rareNR.css,method = "PCoA",distance = randoml_rareNR_wunifrac.dist)

#PCoA-uwunifrac
subseth_uwunifrac_rareNR.pcoa <- ordinate(sh_16S_rareNR.css,method = "PCoA",distance = subseth_uwunifrac.dist)
subsetl_uwunifrac_rareNR.pcoa <- ordinate(sl_16S_rareNR.css,method = "PCoA",distance = subsetl_uwunifrac.dist)
randomh_uwunifrac_rareNR.pcoa <- ordinate(HR_16S_rareNR.css,method = "PCoA",distance = randomh_uwunifrac.dist)
randoml_uwunifrac_rareNR.pcoa <- ordinate(LR_16S_rareNR.css,method = "PCoA",distance = randoml_uwunifrac.dist)


# plot weighted unifrac

pool_col <- c("1","3","6","12")
Number.in.Pool <- c("1","3","6","12")

#### findingcenters
#hs

subseth_wunifrac_rareNR_plot1 <- ordiplot(subseth_wunifrac_rareNR.ord$points)
subseth_wunifrac_rareNR_pcoaplot <- ordiplot(subseth_wunifrac_rareNR.pcoa$vectors)
subseth_wunifrac_rareNR_siteslong <-sites.long(subseth_wunifrac_rareNR_plot1,sh_16S_rareNR.css)
subseth_wunifrac_rareNR_pcoasites <- sites.long(subseth_wunifrac_rareNR_pcoaplot, sh_16S_rareNR.css)


subseth_wunifrac_rareNR_centroids <- envfit(subseth_wunifrac_rareNR.ord~sh_16S_rareNR.css$Number.in.Pool)
subseth_wunifrac_rareNR_centroids
subseth_wunifrac_rareNR_pcoa_centroids <- envfit(subseth_wunifrac_rareNR.pcoa$vectors~sh_16S_rareNR.css.df$Number.in.Pool)
subseth_wunifrac_rareNR_pcoa_centroids


subseth_wunifrac_rareNR_NMDS_col1 <- c(-0.0042,0.0101,0.0422,0.0040)
subseth_wunifrac_rareNR_NMDS_col2 <- c(-0.0008,0.0134,0.0053,-0.0089)
subseth_wunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,subseth_wunifrac_rareNR_NMDS_col1,subseth_wunifrac_rareNR_NMDS_col2)


##ls
subsetl_wunifrac_rareNR_plot1 <- ordiplot(subsetl_wunifrac_rareNR.ord$points)
subsetl_wunifrac_rareNR_siteslong <- sites.long(subsetl_wunifrac_rareNR_plot1,sl_16S_rareNR.css.df)
sl_16S_rareNR.css.df$Number.in.Pool <- as.factor(sl_16S_rareNR.css.df$Number.in.Pool)
subsetl_wunifrac_rareNR_centroids <- envfit(subsetl_wunifrac_rareNR.ord~sl_16s_rareNR.css.df$Number.in.Pool)
subsetl_wunifrac_rareNR_centroids

subsetl_wunifrac_rareNR_NMDS_col1 <- c(0.0019,0.0240,-0.0028,-0.0434)
subsetl_wunifrac_rareNR_NMDS_col2 <- c(0.0010,0.0298,0.0080,-0.0491)
subsetl_wunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,subsetl_wunifrac_rareNR_NMDS_col1,subsetl_wunifrac_rareNR_NMDS_col2)

##hr
HR_16S_rareNR.css.df$Number.in.Pool <- as.factor(HR_16S_rareNR.css.df$Number.in.Pool)
randomh_wunifrac_rareNR_plot1 <- ordiplot(randomh_wunifrac_rareNR.ord$points)
randomh_wunifrac_rareNR_siteslong <- sites.long(randomh_wunifrac_rareNR_plot1,HR_16S_rareNR.css.df)

randomh_wunifrac_rareNR_centroids <- envfit(randomh_wunifrac_rareNR.ord~HR_16S_rareNR.css.df$Number.in.Pool)
randomh_wunifrac_rareNR_centroids

randomh_wunifrac_rareNR_NMDS_col1 <- c( 0.0011,-0.0094,-0.0050,-0.0018)
randomh_wunifrac_rareNR_NMDS_col2 <- c(0.0015,-0.0193,-0.0008,-0.0010)
randomh_wunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,randomh_wunifrac_rareNR_NMDS_col1,randomh_wunifrac_rareNR_NMDS_col2)

##lr
LR_16S_rareNR.css.df$Number.in.Pool <- as.factor(LR_16S_rareNR.css.df$Number.in.Pool)
randoml_wunifrac_rareNR_plot1 <- ordiplot(randoml_wunifrac_rareNR.ord$points)
randoml_wunifrac_rareNR_siteslong <- sites.long(randoml_wunifrac_rareNR_plot1,LR_16S_rareNR.css.df)

randoml_wunifrac_rareNR_centroids <- envfit(randoml_wunifrac_rareNR.ord~LR_16S_rareNR.css.df$Number.in.Pool)
randoml_wunifrac_rareNR_centroids

randoml_wunifrac_rareNR_NMDS_col1 <- c(0.0027,-0.0083,0.0031,-0.0264)
randoml_wunifrac_rareNR_NMDS_col2 <- c(0.0056,0.0033,-0.0135,-0.0556)
randoml_wunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,randoml_wunifrac_rareNR_NMDS_col1,randoml_wunifrac_rareNR_NMDS_col2)



##plot+stats
###hs
hs_wunifrac_rareNR_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Weighted Unifrac") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subseth_wunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2)) +
  stat_ellipse(data = subseth_wunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subseth_wunifrac_rareNR_centroids.df, aes(x=subseth_wunifrac_rareNR_NMDS_col1, y=subseth_wunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subseth_wunifrac_rareNR_centroids.df, aes(x=subseth_wunifrac_rareNR_NMDS_col1, y=subseth_wunifrac_rareNR_NMDS_col2, label = subseth_wunifrac_rareNR_col), colour = "white", size = 2, fontface = "bold") +
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



adonis2(subseth_rareNR_wunifrac.dist ~ Number.in.Pool, sh_16S_rareNR.css.df, nperm = 9999)

hs_wunifrac_rareNR_pa <- pairwise.adonis2(subseth_rareNR_wunifrac.dist ~ Number.in.Pool, sh_16S_rareNR.css.df, nperm = 9999, p.adjust.methods="BH") # NS

write.csv(hs_wunifrac_rareNR_pa,"hs_wunifrac_rareNR_permanova.csv")


subseth_wunifrac_rareNR.disper <- betadisper(subseth_rareNR_wunifrac.dist, sh_16S_rareNR.css.df$Number.in.Pool)
plot(subseth_wunifrac_rareNR.disper)
boxplot(subseth_wunifrac_rareNR.disper)
TukeyHSD(subseth_wunifrac_rareNR.disper)
anova(subseth_wunifrac_rareNR.disper)
permutest(subseth_wunifrac_rareNR.disper, permutations = 9999, pairwise=F)
subseth_wunifrac_rareNR.permdisp <- permutest(subseth_wunifrac_rareNR.disper, permutations = 9999, pairwise = T)
subseth_wunifrac_rareNR.permdisp #NS
write.csv(subseth_wunifrac_rareNR.permdisp[["pairwise"]][["permuted"]],"hs_wunifrac_rareNR_permdisp.csv")

#PCoA
plot_ordination(subseth.css, subseth_wunifrac_rareNR.pcoa, type = "samples", color = "Number.in.Pool") +
  theme_bw() +
  #labs(title = "RUMEN", x= "Axis 1 (30.2% variation explained)", y= "Axis 2 (21.1% variation explained)") +
  geom_point(size = 5, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill = Number.in.Pool), alpha = 0.5, lty = 2, size = 1) +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour= "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))


###ls
subsetl_wunifrac_rareNR_siteslong$Number.in.Pool <- as.factor(subsetl_wunifrac_rareNR_siteslong$Number.in.Pool)
ls_wunifrac_rareNR_nmdsplot <-ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subsetl_wunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = subsetl_wunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subsetl_wunifrac_rareNR_centroids.df, aes(x=subsetl_wunifrac_rareNR_NMDS_col1, y=subsetl_wunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subsetl_wunifrac_rareNR_centroids.df, aes(x=subsetl_wunifrac_rareNR_NMDS_col1, y=subsetl_wunifrac_rareNR_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6, colour = "black"),
        plot.title=element_blank())

ggsave("ls_wunifrac_rareNR.tiff", plot = ls_wunifrac_rareNR_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(subsetl_rareNR_wunifrac.dist ~ Number.in.Pool, sl_16s_rareNR.css.df, nperm = 9999) # NS
ls_wunifrac_rareNR_pa <- pairwise.adonis2(subsetl_rareNR_wunifrac.dist ~ Number.in.Pool, sl_16s_rareNR.css.df, p.adjust.methods = "BH", nperm = 9999)
write.csv(ls_wunifrac_rareNR_pa,"ls_wunifrac_rareNR_permanova.csv")


subsetl_wunifrac_rareNR.disper <- betadisper(subsetl_rareNR_wunifrac.dist, sl_16s_rareNR.css.df$Number.in.Pool)
plot(subsetl_wunifrac_rareNR.disper)
boxplot(subsetl_wunifrac_rareNR.disper)
TukeyHSD(subsetl_wunifrac_rareNR.disper)
anova(subsetl_wunifrac_rareNR.disper)
permutest(subsetl_wunifrac_rareNR.disper, permutations = 9999, pairwise = F)
subsetl_wunifrac_rareNR.permdisp <- permutest(subsetl_wunifrac_rareNR.disper, permutations = 9999, pairwise = T)
subsetl_wunifrac_rareNR.permdisp #no significant difference in homogeneity of dispersion (overall 0.22, pairwise 6 and 1 different)

###hr
randomh_wunifrac_rareNR_siteslong$Number.in.Pool <- as.factor(randomh_wunifrac_rareNR_siteslong$Number.in.Pool)
hr_wunifrac_rareNR_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randomh_wunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randomh_wunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randomh_wunifrac_rareNR_centroids.df, aes(x=randomh_wunifrac_rareNR_NMDS_col1, y=randomh_wunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = randomh_wunifrac_rareNR_centroids.df, aes(x=randomh_wunifrac_rareNR_NMDS_col1, y=randomh_wunifrac_rareNR_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.title.x= element_blank(),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("hr_wunifrac_rareNR.tiff", plot = hr_wunifrac_rareNR_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(randomh_rareNR_wunifrac.dist ~ Number.in.Pool, HR_16S_rareNR.css.df, p.adjust.methods = "BH", nperm = 9999) # NS

randomh_wunifrac_rareNR.disper <- betadisper(randomh_rareNR_wunifrac.dist, HR_16S_rareNR.css.df$Number.in.Pool)
plot(randomh_wunifrac_rareNR.disper)
boxplot(randomh_wunifrac_rareNR.disper)
TukeyHSD(randomh_wunifrac_rareNR.disper)
anova(randomh_wunifrac_rareNR.disper)
permutest(randomh_wunifrac_rareNR.disper, permutations = 9999, pairwise = F)
randomh_wunifrac_rareNR.permdisp <- permutest(randomh_wunifrac_rareNR.disper, permutations = 9999, pairwise = T)
randomh_wunifrac_rareNR.permdisp #NS, 0.07 overall, 1v12 p=0.052

###lr
randoml_wunifrac_rareNR_siteslong$Number.in.Pool <- as.factor(randoml_wunifrac_rareNR_siteslong$Number.in.Pool)
lr_wunifrac_rareNR_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randoml_wunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randoml_wunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randoml_wunifrac_rareNR_centroids.df, aes(x=randoml_wunifrac_rareNR_NMDS_col1, y=randoml_wunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = randoml_wunifrac_rareNR_centroids.df, aes(x=randoml_wunifrac_rareNR_NMDS_col1, y=randoml_wunifrac_rareNR_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("lr_wunifrac_rareNR.tiff", plot = lr_wunifrac_rareNR_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm", dpi=600)


adonis2(randoml_rareNR_wunifrac.dist ~ Number.in.Pool, LR_16S_rareNR.css.df, nperm = 9999) # NS, p=0.072

randoml_wunifrac_rareNR.disper <-betadisper(randoml_rareNR_wunifrac.dist, LR_16S_rareNR.css.df$Number.in.Pool)
plot(randoml_wunifrac_rareNR.disper)
boxplot(randoml_wunifrac_rareNR.disper)
TukeyHSD(randoml_wunifrac_rareNR.disper)
anova(randoml_wunifrac_rareNR.disper)
permutest(randoml_wunifrac_rareNR.disper, permutations = 9999, pairwise = F)
randoml_wunifrac_rareNR.permdisp <- permutest(randoml_wunifrac_rareNR.disper, permutations = 9999, pairwise = T)
randoml_wunifrac_rareNR.permdisp #NS overall p=0.22

##Unweighted Unifrac
#hs

subseth_uwunifrac_rareNR_plot1 <- ordiplot(subseth_uwunifrac_rareNR.ord$points)
subseth_uwunifrac_rareNR_pcoaplot <- ordiplot(subseth_uwunifrac_rareNR.pcoa$vectors)
subseth_uwunifrac_rareNR_siteslong <-sites.long(subseth_uwunifrac_rareNR_plot1,sh_16S_rareNR.css)
subseth_uwunifrac_rareNR_pcoasites <- sites.long(subseth_uwunifrac_rareNR_pcoaplot, sh_16S_rareNR.css)


subseth_uwunifrac_rareNR_centroids <- envfit(subseth_uwunifrac_rareNR.ord~sh_16S_rareNR.css.df$Number.in.Pool)
subseth_uwunifrac_rareNR_centroids



subseth_uwunifrac_rareNR_NMDS_col1 <- c(-0.0042,0.0101,0.0422,0.0040)
subseth_uwunifrac_rareNR_NMDS_col2 <- c(-0.0008,0.0134,0.0053,-0.0089)
subseth_uwunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,subseth_uwunifrac_rareNR_NMDS_col1,subseth_uwunifrac_rareNR_NMDS_col2)


##ls
subsetl_uwunifrac_rareNR_plot1 <- ordiplot(subsetl_uwunifrac_rareNR.ord$points)
subsetl_uwunifrac_rareNR_siteslong <- sites.long(subsetl_uwunifrac_rareNR_plot1,sl_16S_rareNR.css.df)
sl_16S_rareNR.css.df$Number.in.Pool <- as.factor(sl_16S_rareNR.css.df$Number.in.Pool)
subsetl_uwunifrac_rareNR_centroids <- envfit(subsetl_uwunifrac_rareNR.ord~sl_16s_rareNR.css.df$Number.in.Pool)
subsetl_uwunifrac_rareNR_centroids

subsetl_uwunifrac_rareNR_NMDS_col1 <- c(0.0019,0.0240,-0.0028,-0.0434)
subsetl_uwunifrac_rareNR_NMDS_col2 <- c(0.0010,0.0298,0.0080,-0.0491)
subsetl_uwunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,subsetl_uwunifrac_rareNR_NMDS_col1,subsetl_uwunifrac_rareNR_NMDS_col2)

##hr
HR_16S_rareNR.css.df$Number.in.Pool <- as.factor(HR_16S_rareNR.css.df$Number.in.Pool)
randomh_uwunifrac_rareNR_plot1 <- ordiplot(randomh_uwunifrac_rareNR.ord$points)
randomh_uwunifrac_rareNR_siteslong <- sites.long(randomh_uwunifrac_rareNR_plot1,HR_16S_rareNR.css.df)

randomh_uwunifrac_rareNR_centroids <- envfit(randomh_uwunifrac_rareNR.ord~HR_16S_rareNR.css.df$Number.in.Pool)
randomh_uwunifrac_rareNR_centroids

randomh_uwunifrac_rareNR_NMDS_col1 <- c( 0.0011,-0.0094,-0.0050,-0.0018)
randomh_uwunifrac_rareNR_NMDS_col2 <- c(0.0015,-0.0193,-0.0008,-0.0010)
randomh_uwunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,randomh_uwunifrac_rareNR_NMDS_col1,randomh_uwunifrac_rareNR_NMDS_col2)

##lr
LR_16S_rareNR.css.df$Number.in.Pool <- as.factor(LR_16S_rareNR.css.df$Number.in.Pool)
randoml_uwunifrac_rareNR_plot1 <- ordiplot(randoml_uwunifrac_rareNR.ord$points)
randoml_uwunifrac_rareNR_siteslong <- sites.long(randoml_uwunifrac_rareNR_plot1,LR_16S_rareNR.css.df)

randoml_uwunifrac_rareNR_centroids <- envfit(randoml_uwunifrac_rareNR.ord~LR_16S_rareNR.css.df$Number.in.Pool)
randoml_uwunifrac_rareNR_centroids

randoml_uwunifrac_rareNR_NMDS_col1 <- c(0.0027,-0.0083,0.0031,-0.0264)
randoml_uwunifrac_rareNR_NMDS_col2 <- c(0.0056,0.0033,-0.0135,-0.0556)
randoml_uwunifrac_rareNR_centroids.df <-data.frame(Number.in.Pool,randoml_uwunifrac_rareNR_NMDS_col1,randoml_uwunifrac_rareNR_NMDS_col2)


adonis2(subseth_rareNR_uwunifrac.dist ~ Number.in.Pool, sh_16S_rareNR.css.df, nperm = 9999)

hs_uwunifrac_rareNR_pa <- pairwise.adonis2(subseth_rareNR_uwunifrac.dist ~ Number.in.Pool, sh_16S_rareNR.css.df, nperm = 9999, p.adjust.methods="BH") # NS



subseth_uwunifrac_rareNR.disper <- betadisper(subseth_rareNR_uwunifrac.dist, sh_16S_rareNR.css.df$Number.in.Pool)
plot(subseth_uwunifrac_rareNR.disper)
boxplot(subseth_uwunifrac_rareNR.disper)
TukeyHSD(subseth_uwunifrac_rareNR.disper)
anova(subseth_uwunifrac_rareNR.disper)
permutest(subseth_uwunifrac_rareNR.disper, permutations = 9999, pairwise=F)
subseth_uwunifrac_rareNR.permdisp <- permutest(subseth_uwunifrac_rareNR.disper, permutations = 9999, pairwise = T)
subseth_uwunifrac_rareNR.permdisp #NS
write.csv(subseth_uwunifrac_rareNR.permdisp[["pairwise"]][["permuted"]],"hs_uwunifrac_rareNR_permdisp.csv")

#PCoA
plot_ordination(subseth.css, subseth_uwunifrac_rareNR.pcoa, type = "samples", color = "Number.in.Pool") +
  theme_bw() +
  #labs(title = "RUMEN", x= "Axis 1 (30.2% variation explained)", y= "Axis 2 (21.1% variation explained)") +
  geom_point(size = 5, shape = 18) +
  stat_ellipse(geom = "polygon", aes(fill = Number.in.Pool), alpha = 0.5, lty = 2, size = 1) +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour= "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))


###ls
subsetl_uwunifrac_rareNR_siteslong$Number.in.Pool <- as.factor(subsetl_uwunifrac_rareNR_siteslong$Number.in.Pool)
ls_uwunifrac_rareNR_nmdsplot <-ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subsetl_uwunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = subsetl_uwunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subsetl_uwunifrac_rareNR_centroids.df, aes(x=subsetl_uwunifrac_rareNR_NMDS_col1, y=subsetl_uwunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subsetl_uwunifrac_rareNR_centroids.df, aes(x=subsetl_uwunifrac_rareNR_NMDS_col1, y=subsetl_uwunifrac_rareNR_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6, colour = "black"),
        plot.title=element_blank())

ggsave("ls_uwunifrac_rareNR.tiff", plot = ls_uwunifrac_rareNR_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(subsetl_rareNR_uwunifrac.dist ~ Number.in.Pool, sl_16s_rareNR.css.df, nperm = 9999) # NS
ls_uwunifrac_rareNR_pa <- pairwise.adonis2(subsetl_rareNR_uwunifrac.dist ~ Number.in.Pool, sl_16s_rareNR.css.df, p.adjust.methods = "BH", nperm = 9999)
write.csv(ls_uwunifrac_rareNR_pa,"ls_uwunifrac_rareNR_permanova.csv")


subsetl_uwunifrac_rareNR.disper <- betadisper(subsetl_rareNR_uwunifrac.dist, sl_16s_rareNR.css.df$Number.in.Pool)
plot(subsetl_uwunifrac_rareNR.disper)
boxplot(subsetl_uwunifrac_rareNR.disper)
TukeyHSD(subsetl_uwunifrac_rareNR.disper)
anova(subsetl_uwunifrac_rareNR.disper)
permutest(subsetl_uwunifrac_rareNR.disper, permutations = 9999, pairwise = F)
subsetl_uwunifrac_rareNR.permdisp <- permutest(subsetl_uwunifrac_rareNR.disper, permutations = 9999, pairwise = T)
subsetl_uwunifrac_rareNR.permdisp #no significant difference in homogeneity of dispersion (overall 0.22, pairwise 6 and 1 different)

###hr
randomh_uwunifrac_rareNR_siteslong$Number.in.Pool <- as.factor(randomh_uwunifrac_rareNR_siteslong$Number.in.Pool)
hr_uwunifrac_rareNR_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randomh_uwunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randomh_uwunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randomh_uwunifrac_rareNR_centroids.df, aes(x=randomh_uwunifrac_rareNR_NMDS_col1, y=randomh_uwunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = randomh_uwunifrac_rareNR_centroids.df, aes(x=randomh_uwunifrac_rareNR_NMDS_col1, y=randomh_uwunifrac_rareNR_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.title.x= element_blank(),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("hr_uwunifrac_rareNR.tiff", plot = hr_uwunifrac_rareNR_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(randomh_rareNR_uwunifrac.dist ~ Number.in.Pool, HR_16S_rareNR.css.df, p.adjust.methods = "BH", nperm = 9999) # NS
pairwise.adonis2(randomh_rareNR_uwunifrac.dist ~ Number.in.Pool, HR_16S_rareNR.css.df, p.adjust.methods = "BH", nperm = 9999)

randomh_uwunifrac_rareNR.disper <- betadisper(randomh_rareNR_uwunifrac.dist, HR_16S_rareNR.css.df$Number.in.Pool)
plot(randomh_uwunifrac_rareNR.disper)
boxplot(randomh_uwunifrac_rareNR.disper)
TukeyHSD(randomh_uwunifrac_rareNR.disper)
anova(randomh_uwunifrac_rareNR.disper)
permutest(randomh_uwunifrac_rareNR.disper, permutations = 9999, pairwise = F)
randomh_uwunifrac_rareNR.permdisp <- permutest(randomh_uwunifrac_rareNR.disper, permutations = 9999, pairwise = T)
randomh_uwunifrac_rareNR.permdisp #NS, 0.07 overall, 1v12 p=0.052

###lr
randoml_uwunifrac_rareNR_siteslong$Number.in.Pool <- as.factor(randoml_uwunifrac_rareNR_siteslong$Number.in.Pool)
lr_uwunifrac_rareNR_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randoml_uwunifrac_rareNR_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randoml_uwunifrac_rareNR_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randoml_uwunifrac_rareNR_centroids.df, aes(x=randoml_uwunifrac_rareNR_NMDS_col1, y=randoml_uwunifrac_rareNR_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = randoml_uwunifrac_rareNR_centroids.df, aes(x=randoml_uwunifrac_rareNR_NMDS_col1, y=randoml_uwunifrac_rareNR_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("lr_uwunifrac_rareNR.tiff", plot = lr_uwunifrac_rareNR_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm", dpi=600)


adonis2(randoml_rareNR_uwunifrac.dist ~ Number.in.Pool, LR_16S_rareNR.css.df, nperm = 9999) # NS, p=0.072
pairwise.adonis2(randoml_rareNR_uwunifrac.dist ~ Number.in.Pool, LR_16S_rareNR.css.df, nperm = 9999, p.adjust.methods="BH")

randoml_uwunifrac_rareNR.disper <-betadisper(randoml_rareNR_uwunifrac.dist, LR_16S_rareNR.css.df$Number.in.Pool)
plot(randoml_uwunifrac_rareNR.disper)
boxplot(randoml_uwunifrac_rareNR.disper)
TukeyHSD(randoml_uwunifrac_rareNR.disper)
anova(randoml_uwunifrac_rareNR.disper)
permutest(randoml_uwunifrac_rareNR.disper, permutations = 9999, pairwise = F)
randoml_uwunifrac_rareNR.permdisp <- permutest(randoml_uwunifrac_rareNR.disper, permutations = 9999, pairwise = T)
randoml_uwunifrac_rareNR.permdisp #NS overall p=0.22




####ANCOM-BC####
##DataSetup
BiocManager::install("ANCOMBC")
library(ANCOMBC)

###High subset
View(sample_data(subsethigh_rareNR))

sample_data(subsethigh_rareNR)$HSPool

sample_data(subsethigh_rareNR)$StatsPool<-c("A","B","C","D","E","F","A","B","C","D","E","F","A","B","C","D",
                                            "F","B","E","A","E","B","E","E","F","C","F","A","F","D","D",
                                            "E","C","F","F","F","D","C","A","B","A","F","C","B","F",
                                            "B","F","F","D","D","B","D","E","B","C","E","A","D","D",
                                            "C","E","A","A","B","C","F","D","E","C","A","D","B",
                                            "C","D","F","A","E","E","A","B","A","C","C","C","E","A","D","B")

sample_data(subsethigh_rareNR)$StatsPool <- factor(sample_data(subsethigh_rareNR)$StatsPool, levels = c("A","B","C","D","E","F"))

sample_data(subsethigh_rareNR)$Number.in.Pool <- factor(sample_data(subsethigh_rareNR)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(subsethigh_rareNR)


ancom_Family_16S_HS_rareNR = ancombc2(data = subsethigh_rareNR, assay_name = "counts", tax_level = "Family",
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


dunnet_16S_HS_rareNR = ancom_Family_16S_HS_rareNR$res_dunn


write.csv(dunnet_16S_HS_rareNR,"dunnet_16S_HS_rareNR.csv")
write.csv(BCabundancetable_16S_HS, "BCLogAbundance_16S_HS_rareNR.csv")


df_fig_dunnet_16S_HS1 = dunnet_16S_HS_rareNR %>%
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

diff_16S_HS_3 = dunnet_16S_HS_rareNR %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_16S_HS2 = dunnet_16S_HS_rareNR %>%
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


####Low Subset
sample_data(subsetlow_rareNR)$LSPool

sample_data(subsetlow_rareNR)$StatsPool<-c("A","B","C","D","E","F","A","B","C","D","E","F","A","B","C","D",
                                           "E","F","D","F","E","E","F","A","B","D","F","D","E","C","C","F",
                                           "D","F","C","B","A","E","A","A","E","F","F","F","A","D","D","A",
                                           "E","F","B","A","B","E","A","E","C","A","F","B","C","A","C","C",
                                           "B","C","B","D","C","D","E","E","A","B","D","B","B","F","F","A",
                                           "B","D","E","C","B","D","C","C")

sample_data(subsetlow_rareNR)$StatsPool <- factor(sample_data(subsetlow_rareNR)$StatsPool, levels = c("A","B","C","D","E","F"))

sample_data(subsetlow_rareNR)$Number.in.Pool <- factor(sample_data(subsetlow_rareNR)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(subsetlow_rareNR)


ancom_Family_16S_LS_rareNR = ancombc2(data = subsetlow_rareNR, assay_name = "counts", tax_level = "Family",
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


dunnet_16S_LS_rareNR = ancom_Family_16S_LS_rareNR$res_dunn

write.csv(dunnet_16S_LS_rareNR,"dunnet_16S_LS_rareNR.csv")


df_fig_dunnet_16S_LS1 = dunnet_16S_LS_rareNR %>%
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

diff_16S_LS_3 = dunnet_16S_LS_rareNR %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_16S_LS2 = dunnet_16S_LS_rareNR %>%
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
  scale_color_identity(guide = FALSE)+
  labs(x=NULL,y=NULL, title = "Log fold changes as compared to individual animals")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

###High random
View(sample_data(randomhigh_rareNR))

sample_data(randomhigh_rareNR)$HRPool

sample_data(randomhigh_rareNR)$Number.in.Pool <- factor(sample_data(randomhigh_rareNR)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(randomhigh_rareNR)


ancom_Family_16S_HR_rareNR = ancombc2(data = randomhigh_rareNR, assay_name = "counts", tax_level = "Family",
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


dunnet_16S_HR_rareNR = ancom_Family_16S_HR_rareNR$res_dunn


write.csv(dunnet_16S_HR_rareNR,"dunnet_16S_HR_rareNR.csv")
write.csv(BCabundancetable_16S_HR, "BCLogAbundance_16S_HR_rareNR.csv")


df_fig_dunnet_16S_HR1 = dunnet_16S_HR_rareNR %>%
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

diff_16S_HR_3 = dunnet_16S_HR_rareNR %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_16S_HR2 = dunnet_16S_HR_rareNR %>%
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


####Low random
sample_data(randomlow_rareNR)$LRPool


sample_data(randomlow_rareNR)$Number.in.Pool <- factor(sample_data(randomlow_rareNR)$Number.in.Pool, levels=c("1","3","6","12"))

rank_names(randomlow_rareNR)


ancom_Family_16S_LR_rareNR = ancombc2(data = randomlow_rareNR, assay_name = "counts", tax_level = "Family",
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


dunnet_16S_LR_rareNR = ancom_Family_16S_LR_rareNR$res_dunn

write.csv(dunnet_16S_LR_rareNR,"dunnet_16S_LR_rareNR.csv")
write.csv(BCabundancetable_16S_LR, "BCLogAbundance_16S_LR_rareNR.csv")


df_fig_dunnet_16S_LR1 = dunnet_16S_LR_rareNR %>%
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

diff_16S_LR_3 = dunnet_16S_LR_rareNR %>%
  dplyr::filter(diff_Number.in.Pool3)


df_fig_dunnet_16S_LR2 = dunnet_16S_LR_rareNR %>%
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
