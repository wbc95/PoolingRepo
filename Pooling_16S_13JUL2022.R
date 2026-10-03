###Script from from LivAbs_16S (Lee)--indicated by ### before

### #load packages
library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr);
library(ggplot2); library(btools); library(randomcoloR); library(cowplot);
library(pairwiseAdonis); library(picante); library(gridExtra); library(grid); library(wrapr);
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(rstatix); library(viridis); library(ggprism); library(ggpubr);
library(UpSetR); library(vegan); library(BiodiversityR); library(FSA);


### #set wd
setwd("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling")

### # source some stuff
source("g_unifrac.R")
source("w_unifrac.R")
source("uw_unifrac.R")
source("changeGGTaxaNames.R") ###Scripts from Lee

### # import data
qiimedata <- import_biom("table-with-taxonomy.biom", "tree.nwk", "dna-sequences.fasta")
map_file <- import_qiime_sample_data("16Smetadata_complete.txt") # need to convert the date to date type


# combining sample data with the rest
data <- merge_phyloseq(qiimedata,map_file)

sample_names(data) <- sample_data(data)$MEG_ID
sample_names(data)

##Data exploration
data #163 samples with 65612 ASVS
sum(taxa_sums(data)==0)  # 2031 taxa with no counts, from samples not included in metadata
data <- prune_taxa(taxa_sums(data) > 0, data)
#poolqiimedata <- prune_taxa(taxa_sums(poolqiimedata) > 0, poolqiimedata)--old don't think I need to do this step but going ahead
sum(sample_sums(data)==0)

#sum(sample_sums(poolqiimedata)==0) # still 0
data



### # check the names of our ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) # beauty, now they are named properly
head(tax_table(data))


### # changing the GG style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- changeGGtaxa(tax.data) # this gets rid of the GG format


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


## # lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(data))
readplot <- ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill="skyblue", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) #  Talk to Lee about distribution
readplot



# some QC checks
min(sample_sums(data)) #357
max(sample_sums(data)) # 7207405
mean(sample_sums(data)) # 422817.5 good I think
median(sample_sums(data)) # 275248
sort(sample_sums(data)) # 357-> 7207405

View(sample_data(data))
##Subsetting
swab <- subset_samples(data, Type=="Individual")
any(sample_sums()==0) #no
sum(taxa_sums(swab)==0) # 28366
swab <- prune_taxa(taxa_sums(swab)>0,swab)
sum(taxa_sums(swab)==0) #none




#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

######### overall
## data preparation - calculation and adding metadata
alpha_div1 <- estimate_richness(data, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2 <- estimate_pd(data) # calculating Faith's PD

alpha_div <- cbind(alpha_div1, alpha_div2) # combining Faith's and other metrics
alpha_div # looks good, but have some duplicates so lets trim it a bit
alpha_div <- alpha_div[,c(1:5)]
alpha_div # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div.df <- as(sample_data(data), "data.frame")


alpha_div_meta <- cbind(alpha_div1, alpha_div.df)

sample_data(data)$Number.in.Pool <- as.factor(sample_data(data)$Number.in.Pool)

poolnumber_palette <- distinctColorPalette(k=4)

poolsizecomparisons <- list(c("1","3"),c("1","6"),c("1","12"),c("3","6"),c("3","12"),c("6","12"))

ggplot(alpha_div_meta, aes(x= Number.in.Pool, y= Observed, fill= Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y = 17500)

ggplot(alpha_div_meta, aes(x= Number.in.Pool, y= Shannon, fill= Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Shannon's Diversity Index", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y = 10)

##Trying by pool
Pool12_1<- subset_samples(data, PrevH=="Y"&MEG_ID!="H_12-02"&MEG_ID!="H_12-03"&MEG_ID!="H_12-04"&MEG_ID!="H_12-05"&MEG_ID!="H_12-06"&Number.in.Pool!="3"&Number.in.Pool!="6")
View(sample_data(Pool12_1.css))
Pool12_1
sum(taxa_sums(Pool12_1)==0)#45687


alpha_div1_Pool12_1 <- estimate_richness(Pool12_1, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2_Pool12_1 <- estimate_pd(Pool12_1) # calculating Faith's PD

alpha_div_Pool12_1 <- cbind(alpha_div1_Pool12_1, alpha_div2_Pool12_1) # combining Faith's and other metrics
alpha_div_Pool12_1 # looks good, but have some duplicates so lets trim it a bit
alpha_div_Pool12_1 <- alpha_div_Pool12_1[,c(1:5)]
alpha_div_Pool12_1 # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div_Pool12_1.df <- as(sample_data(Pool12_1), "data.frame")

alpha_div_meta_Pool12_1 <- cbind(alpha_div_Pool12_1, alpha_div_Pool12_1.df)

ggplot(alpha_div_meta_Pool12_1, aes(x= MEG_ID, y= Observed)) + theme_bw() +
  geom_bar(stat = "summary") + geom_errorbar(stat = "summary") 



###Pooling strategy
##Subset
subsetpooling <- subset_samples(data,Pooling.Type!="Random" )
View(sample_data(subsetpooling))
any(sample_sums(subsetpooling)==0)#no
sum(taxa_sums(subsetpooling)==0)#yes, 5843
subsetpooling<- prune_taxa(taxa_sums(subsetpooling) > 0, subsetpooling)
sum(taxa_sums(subsetpooling)==0) # no more!


alpha_div1_subset <- estimate_richness(subsetpooling, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2_subset <- estimate_pd(subsetpooling) # calculating Faith's PD

alpha_div_subset <- cbind(alpha_div1_subset, alpha_div2_subset) # combining Faith's and other metrics
alpha_div_subset # looks good, but have some duplicates so lets trim it a bit
alpha_div_subset <- alpha_div_subset[,c(1:5)]
alpha_div_subset # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div_subset.df <- as(sample_data(subsetpooling), "data.frame")

alpha_div_meta_subset <- cbind(alpha_div1_subset, alpha_div_subset.df)

ggplot(alpha_div_meta_subset, aes(x= Number.in.Pool, y= Observed, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~Number.in.Pool)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=17500)

#all different from individual except pool size 6

ggplot(alpha_div_meta_subset, aes(x= Number.in.Pool, y= Shannon, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Shannon's Diversity", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~Number.in.Pool)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=10)
#All different from individual except 6. 12 different from 6


##Random
randompooling <- subset_samples(data,Pooling.Type!="Subset" )
any(sample_sums(randompooling)==0) #no
sum(taxa_sums(randompooling)==0) #yes, 6297
randompooling <- prune_taxa(taxa_sums(randompooling) > 0, randompooling)


alpha_div1_random <- estimate_richness(randompooling, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2_random <- estimate_pd(randompooling) # calculating Faith's PD

alpha_div_random <- cbind(alpha_div1_random, alpha_div2_random) # combining Faith's and other metrics
alpha_div_random # looks good, but have some duplicates so lets trim it a bit
alpha_div_random <- alpha_div_random[,c(1:5)]
alpha_div_random # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div_random.df <- as(sample_data(randompooling), "data.frame")

alpha_div_meta_random <- cbind(alpha_div1_random, alpha_div_random.df)

ggplot(alpha_div_meta_random, aes(x= Number.in.Pool, y= Observed, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~Prevalence)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=17500)
#6 and 12 different from individual, 12 different from 3

ggplot(alpha_div_meta_random, aes(x= Number.in.Pool, y= Shannon, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Shannon's Diversity", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~PrevH)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=10)
#3 and 6 not different from individual, 12 different from all

###Subsetting by prevalence
data_highprev_subset <- subset_samples(subsetpooling, PrevH=="Y")
View(sample_data(data_highprev_subset))
any(sample_sums(data_highprev_subset)==0)
sum(taxa_sums(data_highprev_subset)==0)
data_highprev_subset <- prune_taxa(taxa_sums(data_highprev_subset)>0, data_highprev_subset)


data_lowprev_subset <- subset_samples(subsetpooling, PrevL=="Y")
any(sample_sums(data_lowprev_subset)==0)
sum(taxa_sums(data_lowprev_subset)==0)
data_lowprev_subset <- prune_taxa(taxa_sums(data_lowprev_subset)>0, data_lowprev_subset)



##high subset

alpha_div1_hsubset <- estimate_richness(data_highprev_subset, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2_hsubset <- estimate_pd(data_highprev_subset) # calculating Faith's PD

alpha_div_hsubset <- cbind(alpha_div1_hsubset, alpha_div2_hsubset) # combining Faith's and other metrics
alpha_div_hsubset # looks good, but have some duplicates so lets trim it a bit
alpha_div_hsubset <- alpha_div_hsubset[,c(1:5)]
alpha_div_hsubset # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div_hsubset.df <- as(sample_data(data_highprev_subset), "data.frame")

alpha_div_meta_hsubset <- cbind(alpha_div1_hsubset, alpha_div_hsubset.df)

ggplot(alpha_div_meta_hsubset, aes(x= Number.in.Pool, y= Observed, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~Prevalence)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=7500)
#Individual different from 3 and 12, others same

ggplot(alpha_div_meta_hsubset, aes(x= Number.in.Pool, y= Shannon, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Shannon's Diversity", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~PrevH)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=10)
#Same as Richness


##low subset
alpha_div1_lsubset <- estimate_richness(data_lowprev_subset, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2_lsubset <- estimate_pd(data_lowprev_subset) # calculating Faith's PD

alpha_div_lsubset <- cbind(alpha_div1_lsubset, alpha_div2_lsubset) # combining Faith's and other metrics
alpha_div_lsubset # looks good, but have some duplicates so lets trim it a bit
alpha_div_lsubset <- alpha_div_lsubset[,c(1:5)]
alpha_div_lsubset # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div_lsubset.df <- as(sample_data(data_lowprev_subset), "data.frame")

alpha_div_meta_lsubset <- cbind(alpha_div1_lsubset, alpha_div_lsubset.df)

ggplot(alpha_div_meta_lsubset, aes(x= Number.in.Pool, y= Observed, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~Prevalence)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=17500)

#Individual different from all, 12 different from 6

ggplot(alpha_div_meta_lsubset, aes(x= Number.in.Pool, y= Shannon, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Shannon's Diversity", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  #facet_wrap(~PrevH)+
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=10)
#Individual different from 1 and 3, all others statistically the same


###PrevalenceSubset Random
data_highprev_random <- subset_samples(randompooling, PrevH=="Y")
any(sample_sums(data_highprev_random)==0)#no
sum(taxa_sums(data_highprev_random)==0)#Yes, 34985
data_highprev_random <- prune_taxa(taxa_sums(data_highprev_random)>0,data_highprev_random)


data_lowprev_random <- subset_samples(randompooling, PrevL=="Y")
any(sample_sums(data_lowprev_random)==0)#no
sum(taxa_sums(data_lowprev_random)==0)#yes, 7694
data_lowprev_random <- prune_taxa(taxa_sums(data_lowprev_random)>0, data_lowprev_random)


##highprev
alpha_div1_hrandom <- estimate_richness(data_highprev_random, measures = c("Observed","Shannon","Simpson","InvSimpson"))
alpha_div2_hrandom <- estimate_pd(data_highprev_random)
alpha_div_hrandom <- cbind(alpha_div1_hrandom, alpha_div2_hrandom)
alpha_div_hrandom
alpha_div_hrandom <- alpha_div_hrandom[,c(1:5)]
alpha_div_hrandom
alpha_div_hrandom.df <- as(sample_data(data_highprev_random), "data.frame")
alpha_div_meta_hrandom <- cbind(alpha_div1_hrandom, alpha_div_hrandom.df)

ggplot(alpha_div_meta_hrandom, aes(x= Number.in.Pool, y= Observed, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  #stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=4000)
#No difference

ggplot(alpha_div_meta_hrandom, aes(x= Number.in.Pool, y= Shannon, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  #stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=6, label.x.npc = 0.5)
#No difference

##lowprev
alpha_div1_lrandom <- estimate_richness(data_lowprev_random, measures = c("Observed","Shannon","Simpson","InvSimpson"))
alpha_div2_lrandom <- estimate_pd(data_lowprev_random)
alpha_div_lrandom <- cbind(alpha_div1_lrandom, alpha_div2_lrandom)
alpha_div_lrandom
alpha_div_lrandom <- alpha_div_lrandom[,c(1:5)]
alpha_div_lrandom
alpha_div_lrandom.df <- as(sample_data(data_lowprev_random), "data.frame")
alpha_div_meta_lrandom <- cbind(alpha_div1_lrandom, alpha_div_lrandom.df)

ggplot(alpha_div_meta_lrandom, aes(x= Number.in.Pool, y= Observed, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=17500)

#Individual different from 6 and 12, 3 different, 12 different from all

ggplot(alpha_div_meta_lrandom, aes(x= Number.in.Pool, y= Shannon, fill=Number.in.Pool, color= Number.in.Pool)) + 
  theme_bw() +
  labs(title = "", y= "Richness", x= "Number in Pool") +
  #scale_y_continuous(expand = c(0.05,0,0.15,0)) +
  scale_fill_manual(values = poolnumber_palette) +
  scale_colour_manual(values = poolnumber_palette) +
  geom_boxplot(alpha = 0.6, size =1) +
  geom_point() +
  #scale_x_discrete(limits=c(1,3,6,12))+
  #geom_bar(stat = "summary", alpha= 0.5, size= 0.75) +
  #geom_errorbar(stat = "summary", width= 0.4, size= 0.75) +
  theme(plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 40),
        axis.title.y = element_text(size = 40, vjust = 2.5),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        #axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_compare_means(comparisons = poolsizecomparisons)+
  stat_compare_means(label.y=10)
#Same as richness, except 6 and 12 same

### from Lee
####################################   CSS TRANSFORM    ###################################
any(taxa_sums(data)==0) # QUADRUPLE CHECKING - no Taxa with 0 counts. Good
data.css <- phyloseq_transform_css(data, log = F) #should we separate pools and individuals before CSS?

any(taxa_sums(subsetpooling)==0)
subsetpooling.css <- phyloseq_transform_css(subsetpooling, log = F)
any(taxa_sums(randompooling)==0)
randompooling.css <- phyloseq_transform_css(randompooling, log = F)
any(taxa_sums(data_highprev_subset)==0)
data_highprev_subset.css <- phyloseq_transform_css(data_highprev_subset, log=F)
data_highprev_random.css <- phyloseq_transform_css(data_highprev_random, log=F)
data_lowprev_subset.css <- phyloseq_transform_css(data_lowprev_subset, log = F)
data_lowprev_random.css <- phyloseq_transform_css(data_lowprev_random, log = F)

plot(sort(sample_sums(data.css), T), type = "h", ylab = "reads", xlab= "samples")
plot(sort(sample_sums(subsetpooling.css), T), type = "h", ylab = "reads", xlab= "samples")
plot(sort(sample_sums(randompooling.css), T), type = "h", ylab = "reads", xlab= "samples")
plot(sort(sample_sums(data_highprev_random.css), T), type = "h", ylab = "reads", xlab= "samples")
plot(sort(sample_sums(data_highprev_subset.css), T), type = "h", ylab = "reads", xlab= "samples")
plot(sort(sample_sums(data_lowprev_subset.css), T), type = "h", ylab = "reads", xlab= "samples")
plot(sort(sample_sums(data_lowprev_random.css), T), type = "h", ylab = "reads", xlab= "samples")

data.css.df <- as(sample_data(data.css), "data.frame")

subsetpooling.css.df <- as(sample_data(subsetpooling.css), "data.frame")
randompooling.css.df <- as(sample_data(randompooling.css), "data.frame")
data_highprev_subset.css.df <- as(sample_data(data_highprev_subset.css), "data.frame")
data_highprev_random.css.df <- as(sample_data(data_highprev_random.css), "data.frame")
data_lowprev_random.css.df <- as(sample_data(data_lowprev_random.css), "data.frame")
data_lowprev_subset.css.df <- as(sample_data(data_lowprev_subset.css), "data.frame")

##All

gunifrac_all.dist <- gunifrac(data.css)
gunifrac_all.ord <- ordinate(data.css, method = "NMDS", distance = gunifrac_all.dist)
##Plot generalized
plot_ordination(data.css, gunifrac_all.ord, type = "samples", color = "Number.in.Pool") + theme_bw() +
  labs(x= "NMDS1",
       y= "NMDS2") +
  geom_point(size = 3) +
  stat_ellipse(geom= "polygon", lty = 2, alpha = 0.3, aes(fill= Number.in.Pool), size = 1, level = 0.99) +
  scale_colour_manual(values = Numberinpool_palette) +
  scale_fill_manual(values = Numberinpool_palette) +
  theme(plot.margin = unit(c(1,1,1,1),"lines"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 32),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 16, colour = "black"),
        axis.ticks = element_line(size = 0.9, colour = "black"))


gunifrac_subset.dist <- gunifrac(subsetpooling.css)
gunifrac_random.dist <- gunifrac(randompooling.css)
gunifrac_hr.dist <- gunifrac(data_highprev_random.css)
gunifrac_hs.dist <- gunifrac(data_highprev_subset.css)
gunifrac_ls.dist <- gunifrac(data_lowprev_subset.css)
gunifrac_lr.dist <- gunifrac(data_lowprev_random.css)






#############################################################################################
##############################         RELTV ABUNDANCE         ##############################
#############################################################################################
#############################################################################################

rel_abund_all <- transform_sample_counts(data.css, function(x) {x/sum(x)} * 100)
plot(sort(sample_sums(rel_abund_all), TRUE), type = "h", ylab = "reads", xlab= "samples")


ra_phylum_all <- tax_glom(rel_abund_all, taxrank = "Phylum", NArm = F) 
ra_phylum_all #57 phyla
ra_phylum_palette <- distinctColorPalette(k=57)
ra_phylum_palette


ra_phylum_filt <- filter_taxa(ra_phylum, function(x) mean(x)>1,T)
ra_phylum_filt #6 phyla > 1%

ra_class <- tax_glom(rel_abund, taxrank = "Class", NArm = F) 
ra_class # 175 classes
ra_class_palette <- distinctColorPalette(k=175)

ra_order <- tax_glom(rel_abund, taxrank = "Order", NArm = F) 
ra_order # 346 orders
ra_order_palette <- distinctColorPalette(k=346)

ra_family <- tax_glom(rel_abund, taxrank = "Family", NArm = F) 
ra_family # 586 families
ra_family_palette <- distinctColorPalette(k=586)

ra_genus <- tax_glom(rel_abund, taxrank = "Genus", NArm = F) 
ra_genus # 1278 genera
ra_genus_filt <- filter_taxa(ra_genus, function(x) mean(x) > 1, T)
ra_genus_filt # 10 over 1%
ra_genus_filt_palette <- distinctColorPalette(k=10)


write.csv(otu_table(ra_phylum),"ra_phylum_otus.csv")
write.csv(tax_table(ra_phylum),"ra_phylum_taxa.csv")
write.csv(otu_table(ra_class),"ra_class_otus.csv")
write.csv(tax_table(ra_class),"ra_class_taxa.csv")
write.csv(otu_table(ra_order),"ra_order_otus.csv")
write.csv(tax_table(ra_order),"ra_order_taxa.csv")
write.csv(otu_table(ra_family),"ra_family_otus.csv")
write.csv(tax_table(ra_family),"ra_family_taxa.csv")
write.csv(otu_table(ra_genus),"ra_genus_otus.csv")
write.csv(tax_table(ra_genus),"ra_genus_taxa.csv")

ra_phylum_melt <- psmelt(ra_phylum)
ra_phylum_filt_melt <- psmelt(ra_phylum_filt)


ra_class_melt <- psmelt(ra_class)

ra_order_melt <- psmelt(ra_order)

ra_genus_melt <- psmelt(ra_genus)

ra_genus_filt_melt <- psmelt(ra_genus_filt)

rabundgenusplot <- ggplot(ra_genus_filt_melt, aes(x= Sample, y= Abundance, fill = Genus)) +
  theme_bw() +
  facet_wrap(~Type,scales = "free_x")+
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = ra_genus_palette) +
  theme(axis.text.x=element_text(angle=90))
rabundgenusplot

rabundphylumplot <- ggplot(ra_phylum_melt, aes(x= Sample, y= Abundance, fill = Phylum)) +
  theme_bw() +
  facet_wrap(~Type,scales = "free_x")+
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = ra_phylum_palette) +
  theme(axis.text.x=element_text(angle=90))
rabundphylumplot



##At genus level pools for each pen look similar. Ropes and pools fairly similar, increases in mycoplasma with pools.
##Higher lactobacillus with ropes...oral samples vs nps
#Why is there [Clostridium] and Clostridium? Why is there 5-7N15?-[Clostridium] low abundance taxa?
#No Pasteurella/Mannheimia?


#rabundfamilyplot <- ggplot(ra_family_filt_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  #theme_bw() +
  #geom_bar(stat = "summary", colour = "black") +
  #scale_fill_manual(values = ra_family_palette) +
  #theme(axis.text.x=element_text(angle=90))#legend.position = "none")
#rabundfamilyplot

#rabundorderplot <- ggplot(ra_order_filt_melt, aes(x= Sample, y= Abundance, fill = Order)) +
  #theme_bw() +
  #geom_bar(stat = "summary", colour = "black") +
  #scale_fill_manual(values = ra_order_palette) +
  #theme(axis.text.x=element_text(angle=90))#legend.position = "none")
#rabundorderplot

#rabundclassplot <- ggplot(ra_class_filt_melt, aes(x= Sample, y= Abundance, fill = Class)) +
  #theme_bw() +
  #geom_bar(stat = "summary", colour = "black") +
  #scale_fill_manual(values = ra_class_palette) +
  #theme(axis.text.x=element_text(angle=90))#legend.position = "none")
#rabundclassplot


#############################################################################################
##############################         ALPHA DIVERSITY         ##############################
#############################################################################################
#############################################################################################

## data preparation - calculation and adding metadata
alpha_div1 <- estimate_richness(data, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2 <- estimate_pd(data) # calculating Faith's PD

alpha_div <- cbind(alpha_div1, alpha_div2) # combining Faith's and other metrics
alpha_div # looks good, but have some duplicates so lets trim it a bit
alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div.df <- as(sample_data(data), "data.frame") # making into DF for metadata
str(alpha_div.df)
alpha_div_meta <- cbind(alpha_div, alpha_div.df)
alpha_div_meta # great now we've got alpha_div values and metadata



ggplot(alpha_div_meta, aes(x= Type, y= PD, colour = Type, fill=Type)) +
  theme_bw() +
  geom_errorbar(stat = "summary", color="black", )+
  geom_bar(stat = "summary") +
  stat_compare_means(label.y = 550)



pool_alpha_div_meta <- alpha_div_meta[which(alpha_div_meta$Type=="Pool"),]

poolprevalenceplot <- ggplot(pool_alpha_div_meta, aes(x= Prevalence, y= Observed, fill = Prevalence, colour = Prevalence)) +
  theme_bw() +
  labs(y= "Observed ASVs") +
  geom_boxplot(alpha = 0.5, size = 0.75) +
  geom_point(size = 2) +
  scale_color_viridis_d(option="H") +
  scale_fill_viridis_d(option="H") +
  scale_y_continuous(expand = c(0,0), limits = c(-2,2800), breaks = c(0,900,1800,2700)) +
  theme(legend.position = "right",
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
        panel.border = element_rect(colour = "black", size = 1.25),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(colour = "black", size = 14))+
  stat_compare_means()
poolprevalenceplot

wilcox.test(pool_alpha_div_meta$Observed~pool_alpha_div_meta$Prevalence, p.adjust.method = "BH")
#P=.000143-different, more ASVs in Low Prevalence

poolprevalencestat <- poolprevalenceplot+stat_compare_means(label.y=2700)
poolprevalencestat

##Looking at Pooling ASVs

# split up swabs and pools - never going to be compared together
swabs <- subset_samples(data, Type =="Individual")
swabs


pools <- subset_samples(data,Type == "Pool")
pools

## DATA EXPLORATION

# pools
pools # we have 60 samples
sum(taxa_sums(pools)==0) # bunch (17872)
pools <- prune_taxa(taxa_sums(pools) > 0, pools) #that should get rid of them
sum(taxa_sums(pools)==0)# now no samples without any ASVs, good!
pools # 60 samples, 45709 taxa

# swabs
swabs # 103 samples
sum(taxa_sums(swabs)==0) #  bunch, 28366
swabs <- prune_taxa(taxa_sums(swabs) > 0, swabs) # they gone
sum(sample_sums(swabs)==0) # no samples without any ASVs, good!
swabs # 103 samples, 35215 taxa



#### some QC checks

# pools
min(sample_sums(pools)) #4561
max(sample_sums(pools)) # 7207405
mean(sample_sums(pools)) # ~720K
median(sample_sums(pools)) # ~595K
sort(sample_sums(pools)) # 4561, 5227, etc - should I trim any?
#pools <- subset_samples(pools, sample_sums(pools) > 2000) 
#pools # what should be trimmed if anythinig? doing nothing right now

# swabs
min(sample_sums(swabs)) # 357
max(sample_sums(swabs)) # 773555
mean(sample_sums(swabs)) # 245850
median(sample_sums(swabs)) # 176964
sort(sample_sums(swabs)) # 357; 106770.... 357 seams low. trimming.
swabs <- subset_samples(swabs, sample_sums(swabs)>2000)
swabs # 102 samples...lost one

##
###alpha diversity-subset
##

#pools
pools_alpha_div1 <- estimate_richness(pools, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
pools_alpha_div2 <- estimate_pd(pools) # calculating Faith's PD

pools_alpha_div <- cbind(pools_alpha_div1, pools_alpha_div2) # combining Faith's and other metrics
pools_alpha_div # looks good, but have some duplicates so lets trim it a bit
pools_alpha_div <- pools_alpha_div[,c(1:5)]
pools_alpha_div
pools_alpha_div.df <- as(sample_data(pools), "data.frame") # making into DF for metadata
str(pools_alpha_div.df)
pools_alpha_div_meta <- cbind(pools_alpha_div, pools_alpha_div.df)
pools_alpha_div_meta # great now we've got pools_alpha_div values and metadata

ggplot(pools_alpha_div_meta, aes(x= Prevalence, y= PD, colour = Prevalence, fill=Prevalence)) +
  theme_bw() +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary")+
  scale_fill_viridis_d(option="H")+
  scale_color_viridis_d(option="H")

#swabs
swabs_alpha_div1 <- estimate_richness(swabs, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
swabs_alpha_div2 <- estimate_pd(swabs) # calculating Faith's PD

swabs_alpha_div <- cbind(swabs_alpha_div1, swabs_alpha_div2) # combining Faith's and other metrics
swabs_alpha_div # looks good, but have some duplicates so lets trim it a bit
swabs_alpha_div <- swabs_alpha_div[,c(1:5)]
swabs_alpha_div
swabs_alpha_div.df <- as(sample_data(swabs), "data.frame") # making into DF for metadata
str(swabs_alpha_div.df)
swabs_alpha_div_meta <- cbind(swabs_alpha_div, swabs_alpha_div.df)
swabs_alpha_div_meta # great now we've got swabs_alpha_div values and metadata

ggplot(swabs_alpha_div_meta, aes(x= Mh.Cult, y= PD, colour = Mh.Cult, fill=Mh.Cult)) +
  theme_bw() +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary")+
  scale_fill_viridis_d(option="H")+
  scale_color_viridis_d(option="H")



#############################################################################################
##############################        POOLING QUESTION         ##############################
#############################################################################################
#############################################################################################

#TEST POOL SETUP

PoolH_12.01 <- subset_samples(data, H_12.01=="Y" )
PoolH_12.01
alpha_div1_PoolH_12.01 <- estimate_richness(PoolH_12.01, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
alpha_div2_PoolH_12.01 <- estimate_pd(PoolH12_1) # calculating Faith's PD
alpha_div_PoolH_12.01 <- cbind(alpha_div1_PoolH12_1, alpha_div2_PoolH12_1) # combining Faith's and other metrics

alpha_div_PoolH_12.01 # looks good, but have some duplicates so lets trim it a bit
alpha_div_PoolH_12.01 <- alpha_div_PoolH_12.01[,c(1:5)]
alpha_div_PoolH_12.01 # great, got rid of the duplicate, now have : ASVs (richness), Shannon, Simpson, InvSimpson, FaithsPD
alpha_div_PoolH_12.01.df <- as(sample_data(PoolH_12.01), "data.frame")

alpha_div_meta_PoolH_12.01 <- cbind(alpha_div_PoolH_12.01, alpha_div_PoolH_12.01.df)

#Richness
ggplot(alpha_div_meta_PoolH_12.01, aes(x= Type, y= Observed)) + theme_bw() +
  geom_bar(stat = "summary") + geom_errorbar(stat = "summary") 

pairwise.wilcox.test(alpha_div_meta_PoolH12_1$Observed, alpha_div_meta_PoolH_12.01$Type, p.adjust.method = "BH")

#Play with error bar delete and colors

#Shannon/Diversity 
ggplot(alpha_div_meta_PoolH_12.01, aes(x= Type, y= Shannon)) + theme_bw() +
  geom_bar(stat = "summary") + geom_errorbar(stat = "summary")


pairwise.wilcox.test(alpha_div_meta_PoolH_12.01$Shannon, alpha_div_meta_PoolH_12.01$Type, p.adjust.method = "BH")

#Set up rest of pools

###Alpha Diversity

###RA Plots
data.css <- phyloseq_transform_css(data, log = F)
PoolH_12.01.css <- subset_samples(data.css, H_12.01=="Y" )
PoolH_12.02.css <- subset_samples(data.css, H_12.02=="Y" )
PoolH_12.03.css <- subset_samples(data.css, H_12.03=="Y" )
PoolH_12.04.css <- subset_samples(data.css, H_12.04=="Y" )
PoolH_12.05.css <- subset_samples(data.css, H_12.05=="Y" )
PoolH_12.06.css <- subset_samples(data.css, H_12.06=="Y" )
PoolHP_06.01.css <- subset_samples(data.css, HP_06.01=="Y" )
PoolHP_06.02.css <- subset_samples(data.css, HP_06.02=="Y" )
PoolHP_06.03.css <- subset_samples(data.css, HP_06.03=="Y" )
PoolHP_06.04.css <- subset_samples(data.css, HP_06.04=="Y" )
PoolHP_06.05.css <- subset_samples(data.css, HP_06.05=="Y" )
PoolHP_06.06.css <- subset_samples(data.css, HP_06.06=="Y" )
PoolHP_03.01.css <- subset_samples(data.css, HP_03.01=="Y" )
PoolHP_03.02.css <- subset_samples(data.css, HP_03.02=="Y" )
PoolHP_03.03.css <- subset_samples(data.css, HP_03.03=="Y" )
PoolHP_03.04.css <- subset_samples(data.css, HP_03.04=="Y" )
PoolHP_03.05.css <- subset_samples(data.css, HP_03.05=="Y" )
PoolHP_03.06.css <- subset_samples(data.css, HP_03.06=="Y" )
PoolHR_06.01.css <- subset_samples(data.css, HR_06.01=="Y" )
PoolHR_06.02.css <- subset_samples(data.css, HR_06.02=="Y" )
PoolHR_06.03.css <- subset_samples(data.css, HR_06.03=="Y" )
PoolHR_06.04.css <- subset_samples(data.css, HR_06.04=="Y" )
PoolHR_06.05.css <- subset_samples(data.css, HR_06.05=="Y" )
PoolHR_06.06.css <- subset_samples(data.css, HR_06.06=="Y" )
PoolHR_03.01.css <- subset_samples(data.css, HR_03.01=="Y" )
PoolHR_03.02.css <- subset_samples(data.css, HR_03.02=="Y" )
PoolHR_03.03.css <- subset_samples(data.css, HR_03.03=="Y" )
PoolHR_03.04.css <- subset_samples(data.css, HR_03.04=="Y" )
PoolHR_03.05.css <- subset_samples(data.css, HR_03.05=="Y" )
PoolHR_03.06.css <- subset_samples(data.css, HR_03.06=="Y" )

PoolL_12.01.css <- subset_samples(data.css, L_12.01=="Y" )
PoolL_12.02.css <- subset_samples(data.css, L_12.02=="Y" )
PoolL_12.03.css <- subset_samples(data.css, L_12.03=="Y" )
PoolL_12.04.css <- subset_samples(data.css, L_12.04=="Y" )
PoolL_12.05.css <- subset_samples(data.css, L_12.05=="Y" )
PoolL_12.06.css <- subset_samples(data.css, L_12.06=="Y" )
PoolLP_06.01.css <- subset_samples(data.css, LP_06.01=="Y" )
PoolLP_06.02.css <- subset_samples(data.css, LP_06.02=="Y" )
PoolLP_06.03.css <- subset_samples(data.css, LP_06.03=="Y" )
PoolLP_06.04.css <- subset_samples(data.css, LP_06.04=="Y" )
PoolLP_06.05.css <- subset_samples(data.css, LP_06.05=="Y" )
PoolLP_06.06.css <- subset_samples(data.css, LP_06.06=="Y" )
PoolLP_03.01.css <- subset_samples(data.css, LP_03.01=="Y" )
PoolLP_03.02.css <- subset_samples(data.css, LP_03.02=="Y" )
PoolLP_03.03.css <- subset_samples(data.css, LP_03.03=="Y" )
PoolLP_03.04.css <- subset_samples(data.css, LP_03.04=="Y" )
PoolLP_03.05.css <- subset_samples(data.css, LP_03.05=="Y" )
PoolLP_03.06.css <- subset_samples(data.css, LP_03.06=="Y" )
PoolLR_06.01.css <- subset_samples(data.css, LR_06.01=="Y" )
PoolLR_06.02.css <- subset_samples(data.css, LR_06.02=="Y" )
PoolLR_06.03.css <- subset_samples(data.css, LR_06.03=="Y" )
PoolLR_06.04.css <- subset_samples(data.css, LR_06.04=="Y" )
PoolLR_06.05.css <- subset_samples(data.css, LR_06.05=="Y" )
PoolLR_06.06.css <- subset_samples(data.css, LR_06.06=="Y" )
PoolLR_03.01.css <- subset_samples(data.css, LR_03.01=="Y" )
PoolLR_03.02.css <- subset_samples(data.css, LR_03.02=="Y" )
PoolLR_03.03.css <- subset_samples(data.css, LR_03.03=="Y" )
PoolLR_03.04.css <- subset_samples(data.css, LR_03.04=="Y" )
PoolLR_03.05.css <- subset_samples(data.css, LR_03.05=="Y" )
PoolLR_03.06.css <- subset_samples(data.css, LR_03.06=="Y" )

highprev.css <- subset_samples(data.css, PrevH=="Y")
highprev.css

lowprev.css <- subset_samples(data.css, PrevL=="Y")
lowprev.css

###Relative Abundance
rel_abun_PoolH_12.01 <- transform_sample_counts(PoolH_12.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolH_12.02 <- transform_sample_counts(PoolH_12.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolH_12.03 <- transform_sample_counts(PoolH_12.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolH_12.04 <- transform_sample_counts(PoolH_12.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolH_12.05 <- transform_sample_counts(PoolH_12.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolH_12.06 <- transform_sample_counts(PoolH_12.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_06.01 <- transform_sample_counts(PoolHP_06.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_06.02 <- transform_sample_counts(PoolHP_06.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_06.03 <- transform_sample_counts(PoolHP_06.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_06.04 <- transform_sample_counts(PoolHP_06.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_06.05 <- transform_sample_counts(PoolHP_06.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_06.06 <- transform_sample_counts(PoolHP_06.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_03.01 <- transform_sample_counts(PoolHP_03.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_03.02 <- transform_sample_counts(PoolHP_03.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_03.03 <- transform_sample_counts(PoolHP_03.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_03.04 <- transform_sample_counts(PoolHP_03.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_03.05 <- transform_sample_counts(PoolHP_03.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHP_03.06 <- transform_sample_counts(PoolHP_03.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_06.01 <- transform_sample_counts(PoolHR_06.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_06.02 <- transform_sample_counts(PoolHR_06.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_06.03 <- transform_sample_counts(PoolHR_06.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_06.04 <- transform_sample_counts(PoolHR_06.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_06.05 <- transform_sample_counts(PoolHR_06.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_06.06 <- transform_sample_counts(PoolHR_06.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_03.01 <- transform_sample_counts(PoolHR_03.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_03.02 <- transform_sample_counts(PoolHR_03.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_03.03 <- transform_sample_counts(PoolHR_03.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_03.04 <- transform_sample_counts(PoolHR_03.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_03.05 <- transform_sample_counts(PoolHR_03.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolHR_03.06 <- transform_sample_counts(PoolHR_03.06.css,function(x) {x/sum(x)}*100)

rel_abun_PoolL_12.01 <- transform_sample_counts(PoolL_12.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolL_12.02 <- transform_sample_counts(PoolL_12.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolL_12.03 <- transform_sample_counts(PoolL_12.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolL_12.04 <- transform_sample_counts(PoolL_12.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolL_12.05 <- transform_sample_counts(PoolL_12.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolL_12.06 <- transform_sample_counts(PoolL_12.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_06.01 <- transform_sample_counts(PoolLP_06.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_06.02 <- transform_sample_counts(PoolLP_06.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_06.03 <- transform_sample_counts(PoolLP_06.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_06.04 <- transform_sample_counts(PoolLP_06.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_06.05 <- transform_sample_counts(PoolLP_06.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_06.06 <- transform_sample_counts(PoolLP_06.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_03.01 <- transform_sample_counts(PoolLP_03.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_03.02 <- transform_sample_counts(PoolLP_03.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_03.03 <- transform_sample_counts(PoolLP_03.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_03.04 <- transform_sample_counts(PoolLP_03.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_03.05 <- transform_sample_counts(PoolLP_03.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLP_03.06 <- transform_sample_counts(PoolLP_03.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_06.01 <- transform_sample_counts(PoolLR_06.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_06.02 <- transform_sample_counts(PoolLR_06.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_06.03 <- transform_sample_counts(PoolLR_06.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_06.04 <- transform_sample_counts(PoolLR_06.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_06.05 <- transform_sample_counts(PoolLR_06.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_06.06 <- transform_sample_counts(PoolLR_06.06.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_03.01 <- transform_sample_counts(PoolLR_03.01.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_03.02 <- transform_sample_counts(PoolLR_03.02.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_03.03 <- transform_sample_counts(PoolLR_03.03.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_03.04 <- transform_sample_counts(PoolLR_03.04.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_03.05 <- transform_sample_counts(PoolLR_03.05.css,function(x) {x/sum(x)}*100)
rel_abun_PoolLR_03.06 <- transform_sample_counts(PoolLR_03.06.css,function(x) {x/sum(x)}*100)

rel_abun_high <- transform_sample_counts(highprev.css,function(x) {x/sum(x)}*100)
rel_abun_low <- transform_sample_counts(lowprev.css,function(x) {x/sum(x)}*100)


#Phyla
RAPoolH_12.01_Phyla <- tax_glom(rel_abun_PoolH_12.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolH_12.02_Phyla <- tax_glom(rel_abun_PoolH_12.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolH_12.03_Phyla <- tax_glom(rel_abun_PoolH_12.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolH_12.04_Phyla <- tax_glom(rel_abun_PoolH_12.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolH_12.05_Phyla <- tax_glom(rel_abun_PoolH_12.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolH_12.06_Phyla <- tax_glom(rel_abun_PoolH_12.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_06.01_Phyla <- tax_glom(rel_abun_PoolHP_06.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_06.02_Phyla <- tax_glom(rel_abun_PoolHP_06.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_06.03_Phyla <- tax_glom(rel_abun_PoolHP_06.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_06.04_Phyla <- tax_glom(rel_abun_PoolHP_06.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_06.05_Phyla <- tax_glom(rel_abun_PoolHP_06.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_06.06_Phyla <- tax_glom(rel_abun_PoolHP_06.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_03.01_Phyla <- tax_glom(rel_abun_PoolHP_03.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_03.02_Phyla <- tax_glom(rel_abun_PoolHP_03.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_03.03_Phyla <- tax_glom(rel_abun_PoolHP_03.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_03.04_Phyla <- tax_glom(rel_abun_PoolHP_03.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_03.05_Phyla <- tax_glom(rel_abun_PoolHP_03.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHP_03.06_Phyla <- tax_glom(rel_abun_PoolHP_03.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_06.01_Phyla <- tax_glom(rel_abun_PoolHR_06.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_06.02_Phyla <- tax_glom(rel_abun_PoolHR_06.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_06.03_Phyla <- tax_glom(rel_abun_PoolHR_06.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_06.04_Phyla <- tax_glom(rel_abun_PoolHR_06.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_06.05_Phyla <- tax_glom(rel_abun_PoolHR_06.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_06.06_Phyla <- tax_glom(rel_abun_PoolHR_06.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_03.01_Phyla <- tax_glom(rel_abun_PoolHR_03.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_03.02_Phyla <- tax_glom(rel_abun_PoolHR_03.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_03.03_Phyla <- tax_glom(rel_abun_PoolHR_03.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_03.04_Phyla <- tax_glom(rel_abun_PoolHR_03.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_03.05_Phyla <- tax_glom(rel_abun_PoolHR_03.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolHR_03.06_Phyla <- tax_glom(rel_abun_PoolHR_03.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()

RAPoolL_12.01_Phyla <- tax_glom(rel_abun_PoolL_12.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolL_12.02_Phyla <- tax_glom(rel_abun_PoolL_12.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolL_12.03_Phyla <- tax_glom(rel_abun_PoolL_12.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolL_12.04_Phyla <- tax_glom(rel_abun_PoolL_12.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolL_12.05_Phyla <- tax_glom(rel_abun_PoolL_12.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolL_12.06_Phyla <- tax_glom(rel_abun_PoolL_12.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_06.01_Phyla <- tax_glom(rel_abun_PoolLP_06.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_06.02_Phyla <- tax_glom(rel_abun_PoolLP_06.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_06.03_Phyla <- tax_glom(rel_abun_PoolLP_06.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_06.04_Phyla <- tax_glom(rel_abun_PoolLP_06.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_06.05_Phyla <- tax_glom(rel_abun_PoolLP_06.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_06.06_Phyla <- tax_glom(rel_abun_PoolLP_06.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_03.01_Phyla <- tax_glom(rel_abun_PoolLP_03.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_03.02_Phyla <- tax_glom(rel_abun_PoolLP_03.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_03.03_Phyla <- tax_glom(rel_abun_PoolLP_03.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_03.04_Phyla <- tax_glom(rel_abun_PoolLP_03.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_03.05_Phyla <- tax_glom(rel_abun_PoolLP_03.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLP_03.06_Phyla <- tax_glom(rel_abun_PoolLP_03.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_06.01_Phyla <- tax_glom(rel_abun_PoolLR_06.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_06.02_Phyla <- tax_glom(rel_abun_PoolLR_06.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_06.03_Phyla <- tax_glom(rel_abun_PoolLR_06.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_06.04_Phyla <- tax_glom(rel_abun_PoolLR_06.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_06.05_Phyla <- tax_glom(rel_abun_PoolLR_06.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_06.06_Phyla <- tax_glom(rel_abun_PoolLR_06.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_03.01_Phyla <- tax_glom(rel_abun_PoolLR_03.01, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_03.02_Phyla <- tax_glom(rel_abun_PoolLR_03.02, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_03.03_Phyla <- tax_glom(rel_abun_PoolLR_03.03, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_03.04_Phyla <- tax_glom(rel_abun_PoolLR_03.04, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_03.05_Phyla <- tax_glom(rel_abun_PoolLR_03.05, taxrank = "Phylum", NArm = F) %>%
  psmelt()
RAPoolLR_03.06_Phyla <- tax_glom(rel_abun_PoolLR_03.06, taxrank = "Phylum", NArm = F) %>%
  psmelt()



high_Phyla <- tax_glom(rel_abun_high, taxrank = "Phylum", NArm = F) %>%
  psmelt()
low_Phyla <- tax_glom(rel_abun_low, taxrank = "Phylum", NArm = F) %>%
  psmelt()

Phylacolors <- distinctColorPalette(57)
#Don't run again if you like it!

raplotH_12.01 <- ggplot(RAPoolH_12.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.01

raplotH_12.02 <- ggplot(RAPoolH_12.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.02

raplotH_12.03 <- ggplot(RAPoolH_12.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.03

raplotH_12.04 <- ggplot(RAPoolH_12.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.04

raplotH_12.05 <- ggplot(RAPoolH_12.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.05

raplotH_12.06 <- ggplot(RAPoolH_12.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.06

raplotHP_06.01 <- ggplot(RAPoolHP_06.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.01

raplotHP_06.02 <- ggplot(RAPoolHP_06.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.02

raplotHP_06.03 <- ggplot(RAPoolHP_06.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.03

raplotHP_06.04 <- ggplot(RAPoolHP_06.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.04

raplotHP_06.05 <- ggplot(RAPoolHP_06.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.05

raplotHP_06.06 <- ggplot(RAPoolHP_06.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.06

raplotHP_03.01 <- ggplot(RAPoolHP_03.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.01

raplotHP_03.02 <- ggplot(RAPoolHP_03.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.02

raplotHP_03.03 <- ggplot(RAPoolHP_03.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.03

raplotHP_03.04 <- ggplot(RAPoolHP_03.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.04

raplotHP_03.05 <- ggplot(RAPoolHP_03.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.05

raplotHP_03.06 <- ggplot(RAPoolHP_03.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.06

raplotHR_06.01 <- ggplot(RAPoolHR_06.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.01

raplotHR_06.02 <- ggplot(RAPoolHR_06.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.02

raplotHR_06.03 <- ggplot(RAPoolHR_06.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.03

raplotHR_06.04 <- ggplot(RAPoolHR_06.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.04

raplotHR_06.05 <- ggplot(RAPoolHR_06.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.05

raplotHR_06.06 <- ggplot(RAPoolHR_06.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.06

raplotHR_03.01 <- ggplot(RAPoolHR_03.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.01

raplotHR_03.02 <- ggplot(RAPoolHR_03.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.02

raplotHR_03.03 <- ggplot(RAPoolHR_03.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.03

raplotHR_03.04 <- ggplot(RAPoolHR_03.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.04

raplotHR_03.05 <- ggplot(RAPoolHR_03.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.05

raplotHR_03.06 <- ggplot(RAPoolHR_03.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.06

raplotL_12.01 <- ggplot(RAPoolL_12.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.01

raplotL_12.02 <- ggplot(RAPoolL_12.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.02

raplotL_12.03 <- ggplot(RAPoolL_12.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.03

raplotL_12.04 <- ggplot(RAPoolL_12.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.04

raplotL_12.05 <- ggplot(RAPoolL_12.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.05

raplotL_12.06 <- ggplot(RAPoolL_12.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.06

raplotLP_06.01 <- ggplot(RAPoolLP_06.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.01

raplotLP_06.02 <- ggplot(RAPoolLP_06.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.02

raplotLP_06.03 <- ggplot(RAPoolLP_06.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.03

raplotLP_06.04 <- ggplot(RAPoolLP_06.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.04

raplotLP_06.05 <- ggplot(RAPoolLP_06.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.05

raplotLP_06.06 <- ggplot(RAPoolLP_06.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.06

raplotLP_03.01 <- ggplot(RAPoolLP_03.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.01

raplotLP_03.02 <- ggplot(RAPoolLP_03.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.02

raplotLP_03.03 <- ggplot(RAPoolLP_03.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.03

raplotLP_03.04 <- ggplot(RAPoolLP_03.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.04

raplotLP_03.05 <- ggplot(RAPoolLP_03.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.05

raplotLP_03.06 <- ggplot(RAPoolLP_03.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.06

raplotLR_06.01 <- ggplot(RAPoolLR_06.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.01

raplotLR_06.02 <- ggplot(RAPoolLR_06.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.02

raplotLR_06.03 <- ggplot(RAPoolLR_06.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.03

raplotLR_06.04 <- ggplot(RAPoolLR_06.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.04

raplotLR_06.05 <- ggplot(RAPoolLR_06.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.05

raplotLR_06.06 <- ggplot(RAPoolLR_06.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.06

raplotLR_03.01 <- ggplot(RAPoolLR_03.01_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.01", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.01

raplotLR_03.02 <- ggplot(RAPoolLR_03.02_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.02", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.02

raplotLR_03.03 <- ggplot(RAPoolLR_03.03_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.03", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.03

raplotLR_03.04 <- ggplot(RAPoolLR_03.04_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.04", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.04

raplotLR_03.05 <- ggplot(RAPoolLR_03.05_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.05", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.05

raplotLR_03.06 <- ggplot(RAPoolLR_03.06_Phyla, aes(x = Type, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.06", y="Relative Abundance (%)") +
  scale_x_discrete(limits=c("Individual","Pool"),
                   labels = c("Individual","Pool")) +
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.06

##Separate by individual samples

raplotH_12.01_individual <- ggplot(RAPoolH_12.01_Phyla, aes(x = MEG_ID, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.01_individual

raplotH_12.02_individual <- ggplot(RAPoolH_12.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.02_individual

raplotH_12.03_individual <- ggplot(RAPoolH_12.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.03_individual

raplotH_12.04_individual <- ggplot(RAPoolH_12.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.04_individual

raplotH_12.05_individual <- ggplot(RAPoolH_12.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.05_individual

raplotH_12.06_individual <- ggplot(RAPoolH_12.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL H_12.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotH_12.06_individual

raplotHP_06.01_individual <- ggplot(RAPoolHP_06.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.01_individual

raplotHP_06.02_individual <- ggplot(RAPoolHP_06.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.02_individual

raplotHP_06.03_individual <- ggplot(RAPoolHP_06.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.03_individual

raplotHP_06.04_individual <- ggplot(RAPoolHP_06.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.04_individual

raplotHP_06.05_individual <- ggplot(RAPoolHP_06.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.05_individual

raplotHP_06.06_individual <- ggplot(RAPoolHP_06.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_06.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_06.06_individual

raplotHP_03.01_individual <- ggplot(RAPoolHP_03.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.01_individual

raplotHP_03.02_individual <- ggplot(RAPoolHP_03.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.02_individual

raplotHP_03.03_individual <- ggplot(RAPoolHP_03.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.03_individual

raplotHP_03.04_individual <- ggplot(RAPoolHP_03.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.04_individual

raplotHP_03.05_individual <- ggplot(RAPoolHP_03.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.05_individual

raplotHP_03.06_individual <- ggplot(RAPoolHP_03.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HP_03.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHP_03.06_individual

raplotHR_06.01_individual <- ggplot(RAPoolHR_06.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.01_individual

raplotHR_06.02_individual <- ggplot(RAPoolHR_06.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.02_individual

raplotHR_06.03_individual <- ggplot(RAPoolHR_06.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.03_individual

raplotHR_06.04_individual <- ggplot(RAPoolHR_06.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.04_individual

raplotHR_06.05_individual <- ggplot(RAPoolHR_06.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.05_individual

raplotHR_06.06_individual <- ggplot(RAPoolHR_06.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_06.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_06.06_individual

raplotHR_03.01_individual <- ggplot(RAPoolHR_03.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.01_individual

raplotHR_03.02_individual <- ggplot(RAPoolHR_03.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.02_individual

raplotHR_03.03_individual <- ggplot(RAPoolHR_03.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.03_individual

raplotHR_03.04_individual <- ggplot(RAPoolHR_03.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.04_individual

raplotHR_03.05_individual <- ggplot(RAPoolHR_03.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.05_individual

raplotHR_03.06_individual <- ggplot(RAPoolHR_03.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL HR_03.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotHR_03.06_individual

raplotL_12.01_individual <- ggplot(RAPoolL_12.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.01_individual

raplotL_12.02_individual <- ggplot(RAPoolL_12.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.02_individuald

raplotL_12.03_individual <- ggplot(RAPoolL_12.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.03_individual

raplotL_12.04_individual <- ggplot(RAPoolL_12.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.04_individual

raplotL_12.05_individual <- ggplot(RAPoolL_12.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.05_individual

raplotL_12.06_individual <- ggplot(RAPoolL_12.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL L_12.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotL_12.06_individual

raplotLP_06.01_individual <- ggplot(RAPoolLP_06.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.01_individual

raplotLP_06.02_individual <- ggplot(RAPoolLP_06.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.02_individual

raplotLP_06.03_individual <- ggplot(RAPoolLP_06.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.03_individual

raplotLP_06.04_individual <- ggplot(RAPoolLP_06.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.04_individual

raplotLP_06.05_individual <- ggplot(RAPoolLP_06.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.05_individual

raplotLP_06.06_individual <- ggplot(RAPoolLP_06.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_06.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_06.06_individual

raplotLP_03.01_individual <- ggplot(RAPoolLP_03.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.01_individual

raplotLP_03.02_individual <- ggplot(RAPoolLP_03.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.02_individual

raplotLP_03.03_individual <- ggplot(RAPoolLP_03.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.03_individual

raplotLP_03.04_individual <- ggplot(RAPoolLP_03.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.04_individual

raplotLP_03.05_individual <- ggplot(RAPoolLP_03.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.05_individual

raplotLP_03.06_individual <- ggplot(RAPoolLP_03.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LP_03.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLP_03.06_individual

raplotLR_06.01_individual <- ggplot(RAPoolLR_06.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.01_individual

raplotLR_06.02_individual <- ggplot(RAPoolLR_06.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.02_individual

raplotLR_06.03_individual <- ggplot(RAPoolLR_06.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.03_individual

raplotLR_06.04_individual <- ggplot(RAPoolLR_06.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.04_individual

raplotLR_06.05_individual <- ggplot(RAPoolLR_06.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.05_individual

raplotLR_06.06_individual <- ggplot(RAPoolLR_06.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_06.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_06.06_individual

raplotLR_03.01_individual <- ggplot(RAPoolLR_03.01_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.01", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.01_individual

raplotLR_03.02_individual <- ggplot(RAPoolLR_03.02_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.02", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.02_individual

raplotLR_03.03_individual <- ggplot(RAPoolLR_03.03_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.03", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.03_individual

raplotLR_03.04_individual <- ggplot(RAPoolLR_03.04_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.04", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.04_individual

raplotLR_03.05_individual <- ggplot(RAPoolLR_03.05_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.05", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.05_individual

raplotLR_03.06_individual <- ggplot(RAPoolLR_03.06_Phyla, aes(x = sample.id, y= Abundance, fill= Phylum)) +
  theme_bw() +
  labs(title = "POOL LR_03.06", y="Relative Abundance (%)") +
  
                   
  geom_bar(stat = "summary", color= "black")+
  scale_fill_manual(values = Phylacolors)+
  theme(#legend.position = "none",
    plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 48),
    axis.title.y = element_text(size = 40, vjust = 2.5),
    axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = 0.95, vjust = 0.95),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank())
raplotLR_03.06_individual

### MANNHEIMIA
samples_ra_genus <- tax_glom(samples.css.ra, taxrank = "Genus", NArm = F) # 1172 genera
samples_ra_genus_filt <- merge_low_abundance(samples_ra_genus, threshold = 0.05) # 118
samples_ra_genus_filt_melt <- psmelt(samples_ra_genus_filt)

mannheimia <- subset_taxa(samples_ra_genus, Genus=="Mannheimia")
mannheimia_melt <- psmelt(mannheimia)

mannheimia_plot_palette <- c("dodgerblue3","orange3","dodgerblue1","mediumorchid")

ggplot(mannheimia_melt, aes(x= matrix, y= Abundance, fill = matrix, colour = matrix)) +
  theme_bw() + facet_wrap(~pen, nrow = 1) +
  labs(y= "Relative Abundance (%)", title = "Mannheimia") +
  geom_boxplot(alpha= 0.1) +
  #geom_bar(stat = "summary", alpha = 0.1) +
  #geom_errorbar(stat = "summary") +
  geom_point(shape = 18, size = 4) +
  scale_fill_manual(values = mannheimia_plot_palette) +
  scale_colour_manual(values = mannheimia_plot_palette) +
  scale_x_discrete(limits = c("swab","pooled swab","rope","water bowl")) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 22, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

mannheimia_pen8 <- mannheimia_melt[which(mannheimia_melt$pen=="eight"),]
dunnTest(mannheimia_pen8$Abundance, mannheimia_pen8$matrix, method = "bh") #NS

mannheimia_pen11 <- mannheimia_melt[which(mannheimia_melt$pen=="eleven"),]
dunnTest(mannheimia_pen11$Abundance, mannheimia_pen11$matrix, method = "bh") #NS

mannheimia_pen5 <- mannheimia_melt[which(mannheimia_melt$pen=="five"),]
dunnTest(mannheimia_pen5$Abundance, mannheimia_pen5$matrix, method = "bh") #NS

mannheimia_pen4 <- mannheimia_melt[which(mannheimia_melt$pen=="four"),]
dunnTest(mannheimia_pen4$Abundance, mannheimia_pen4$matrix, method = "bh") #NS

mannheimia_pen9 <- mannheimia_melt[which(mannheimia_melt$pen=="nine"),]
dunnTest(mannheimia_pen9$Abundance, mannheimia_pen9$matrix, method = "bh") #NS

mannheimia_pen1 <- mannheimia_melt[which(mannheimia_melt$pen=="one"),]
dunnTest(mannheimia_pen1$Abundance, mannheimia_pen1$matrix, method = "bh") # all but p.swab v i.swab

mannheimia_pen7 <- mannheimia_melt[which(mannheimia_melt$pen=="seven"),]
dunnTest(mannheimia_pen7$Abundance, mannheimia_pen7$matrix, method = "bh") #NS

mannheimia_pen6 <- mannheimia_melt[which(mannheimia_melt$pen=="six"),]
dunnTest(mannheimia_pen6$Abundance, mannheimia_pen6$matrix, method = "bh") #NS

mannheimia_pen3 <- mannheimia_melt[which(mannheimia_melt$pen=="three"),]
dunnTest(mannheimia_pen3$Abundance, mannheimia_pen3$matrix, method = "bh") #NS

mannheimia_pen2 <- mannheimia_melt[which(mannheimia_melt$pen=="two"),]
dunnTest(mannheimia_pen2$Abundance, mannheimia_pen2$matrix, method = "bh") #NS

### Mh by pool -> ind. vs pool
mannheimia_swab_only <- mannheimia_melt[which(mannheimia_melt$matrix=="swab" | mannheimia_melt$matrix =="pooled swab"),]

pool_palette <- c("dodgerblue3","dodgerblue1")

ggplot(mannheimia_swab_only, aes(x= matrix, y= Abundance, fill = matrix, colour = matrix)) +
  theme_bw() + facet_wrap(~pool, nrow = 1) +
  labs(y= "Relative Abundance (%)", title = "16S rRNA GENE SEQUENCING") +
  geom_boxplot(alpha= 0.1) +
  #geom_bar(stat = "summary", alpha = 0.25) +
  #geom_errorbar(stat = "summary") +
  geom_point(shape = 18, size = 4) +
  scale_fill_manual(values = pool_palette) +
  scale_colour_manual(values = pool_palette) +
  scale_x_discrete(limits = c("swab","pooled swab")) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour = "black", size = 22, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_blank())

