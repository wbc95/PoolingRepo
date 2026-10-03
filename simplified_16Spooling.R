library(phyloseq);library(metagenomeSeq);library(dplyr);library(scales);
library(pairwiseAdonis); library(vegan); library(metagMisc); library(stringr);
library(ggplot2); library(btools); library(randomcoloR); library(cowplot);
library(pairwiseAdonis); library(picante); library(gridExtra); library(grid); library(wrapr);
library(ggalt); library(ggforce); library(concaveman); library(ggdendro);
library(microbiome); library(rstatix); library(viridis); library(ggprism); library(ggpubr);


setwd("~/Documents/Research/Research Projects/USDA NIFA AMR/Pooling Pilot/16S Pooling")

### # source some stuff
source("g_unifrac.R")
source("w_unifrac.R")
source("uw_unifrac.R")
source("changeGGTaxaNames.R") ###Scripts from Lee

qiimedata <- import_biom("table-with-taxonomy.biom", "tree.nwk", "dna-sequences.fasta")
map_file <- import_qiime_sample_data("16Smetadata_complete.txt") # need to convert the date to date type


# combining sample data with the rest
data <- merge_phyloseq(qiimedata,map_file)

##Data exploration
data #163 samples with 65612 ASVS
sum(taxa_sums(data)==0) # 2031 taxa with no counts
data <- prune_taxa(taxa_sums(data) > 0, data)
sum(sample_sums(data)==0)

### # check the names of our ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) # beauty, now they are named properly

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
tail(tax_table(data), 20) # "unclassified Unassigned" is that what we want?
#heading names changed

sample_names(data) <- sample_data(data)$MEG_ID
sample_names(data)
sample_data(data)$sample.id <- sample_data(data)$MEG_ID

## # lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(data))
readplot <- ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 5000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) #  Talk to Lee about distribution
readplot #a few pretty low. One very high

# some QC checks
min(sample_sums(data)) #357
max(sample_sums(data)) # 7207405
mean(sample_sums(data)) # 422817.5 good I think
median(sample_sums(data)) # 275248
sort(sample_sums(data)) # 357-> 7207405, do we need to cut some more out?

####################################   CSS TRANSFORM    ###################################
any(taxa_sums(data)==0) # QUADRUPLE CHECKING - no Taxa with 0 counts. Good
data.css <- phyloseq_transform_css(data, log = F)
data.css.ra <- transform_sample_counts(data.css, function(x) {x/sum(x)} * 100)

data.css.df <- as(sample_data(data.css), "data.frame")


###############################################################################
############################ CHANGES IN THE RAs of ASVs #######################

ASVs <- as(phyloseq_to_df(data.css.ra),"data.frame")
colnames(ASVs)
rownames(ASVs) <- ASVs$OTU
ASVs <- ASVs[,-1]

#write.csv(tax_table(ASVs_top100), "top100ASVs_taxa.csv")
#write.csv(otu_table(ASVs_top100), "top100ASVs_otus.csv")

length(unique(ASVs$Phylum)) # 57 phyla
ASV_phyla_palette <- distinctColorPalette(57)

length(unique(ASVs$Class)) # 175 classes
ASV_class_palette <- distinctColorPalette(175)

length(unique(ASVs$Order)) # 344 orders
ASV_order_palette <- distinctColorPalette(344)

length(unique(ASVs$Family)) # 584 familes
ASV_family_palette <- distinctColorPalette(584)

length(unique(ASVs$Genus)) # 1263 genera
ASV_genus_palette <- distinctColorPalette(1263)

###### ASV correlation plots
#### High Prevalence
### Pools of 12

#12.01
corplotASV_PoolH_12.01 <- ggplot(ASVs) + theme_bw() +
  labs(y= "ASV RA (individual)", x= "ASV RA (pool)") +
  geom_point(aes(x = RespSwabPool_H_12_1, y= SwabCompGNP110, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y= SwabCompGNP124, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP144, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP147, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP210, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP219, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP220, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP231, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP242, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP245, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP247, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_1, y = SwabCompGNP256, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_12.01




cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP110) #0.714
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP124) #0.894
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP144) #0.794
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP147) #0.474
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP210) #0.833
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP219) #0.363
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP220) #0.593
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP231) #0.857
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP242) #0.733
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP245) #0.667
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP247) #0.807
cor.test(ASVs$RespSwabPool_H_12_1,ASVs$SwabCompGNP256) #0.091

(0.714 + 0.894 + 0.794 + 0.474 + 0.833 + 0.363 + 0.593 + 0.857 + 0.733 + 0.667 + 0.807 + 0.091)/12 # average 0.652

#12.02
corplotASV_PoolH_12.02 <- ggplot(ASVs) + theme_bw() +
  labs(y= "ASV RA (individual)", x= "ASV RA (pool)") +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP102, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP115, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP120, colour = Phylum), size = 2, shape = 19) + #NO SC_140
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP145, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP152, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP154, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP202, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP206, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP221, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP236, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP246, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_2, y = SwabCompGNP258, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_12.02

cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP102) #0.244
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP115) #0.335
#cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP120) # SC_140 not sequenced?
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP145) #0.355
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP152) #0.422
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP154) #0.367
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP202) #0.359
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP206) #0.235
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP221) #0.868
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP236) #0.476
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP246) #0.910
cor.test(ASVs$RespSwabPool_H_12_2,ASVs$SwabCompGNP258) #0.174

(0.244 + 0.335 + 0.355 + 0.422 + 0.367 + 0.359 + 0.235 + 0.868 + 0.476 + 0.910 + 0.174 )/11 # average 0.431

#12.03
corplotASV_PoolH_12.03 <- ggplot(ASVs) + theme_bw() +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP122, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP131, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP143, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP151, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP207, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP217, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP222, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP229, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP238, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP250, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP253, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_3, y = SwabCompGNP254, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_12.03

cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP122) #0.118
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP131) #0.959
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP143) #0.969
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP151) #0.968
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP207) #0.378
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP217) #0.474
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP222) #0.836
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP229) #0.904
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP238) #0.887
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP250) #0.114
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP253) #0.924
cor.test(ASVs$RespSwabPool_H_12_3,ASVs$SwabCompGNP254) #0.970

(0.118 + 0.959 + 0.969 + 0.968 + 0.378 + 0.474 + 0.836 + 0.904 + 0.887 + 0.114 + 0.924 + 0.970)/12 # average 0.708

#12.04
corplotASV_PoolH_12.04 <- ggplot(ASVs) + theme_bw() +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP127, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP128, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP142, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP160, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP201, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP203, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP213, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP214, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP227, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP235, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP240, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_12_4, y = SwabCompGNP257, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_12.04

cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP127) #0.800
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP128) #0.471
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP142) #0.810
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP160) #0.779
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP201) #0.191
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP203) #0.784
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP213) #0.609
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP214) #0.565
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP227) #0.817
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP235) #0.714
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP240) #0.773
cor.test(ASVs$RespSwabPool_H_12_4,ASVs$SwabCompGNP257) #0.062

#12.05
corplotASV_PoolH_12.05 <- ggplot(ASVs) + theme_bw() +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP107, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP111, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP118, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP119, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP129, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP204, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP209, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP218, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP228, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP243, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP244, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_5, y = SwabCompGNP255, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_12.05

cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP107) #0.306
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP111) #0.316
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP118) #0.383
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP119) #0.362
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP129) #0.196
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP204) #0.924
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP209) #0.689
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP218) #0.323
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP228) #0.377
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP243) #0.401
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP244) #0.278
cor.test(ASVs$RespSwabPool_H_12_5,ASVs$SwabCompGNP255) #0.423


#mean correlation
(0.306+0.316+0.383+0.362+0.196+0.924+0.689+0.323+0.377+0.401+0.278+0.423)/12 #0.415


#12.06
corplotASV_PoolH_12.06 <- 
  ggplot(ASVs) + theme_bw() +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP121, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP123, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP125, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP135, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP137, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP138, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP150, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP153, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP157, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP159, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP226, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_12_6, y = SwabCompGNP241, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_12.06

cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP121) #0.775
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP123) #0.842
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP125) #0.595
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP135) #0.772
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP137) #0.566
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP138) #0.759
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP150) #0.413
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP153) #0.059
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP157) #0.101
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP159) #0.800
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP226) #0.430
cor.test(ASVs$RespSwabPool_H_12_6,ASVs$SwabCompGNP241) #0.596

(0.775+0.842+0.595+0.772+0.566+0.759+0.413+0.059+0.101+0.800+0.430+0.596)/12 #0.559


###Pools of 6
##Subset

#06.01
#corplotASV_PoolHP_06.01 <- 
  ggplot(ASVs) + theme_bw() +
  labs(y= "ASV RA (individual)", x= "ASV RA (pool)") +
  geom_point(aes(x = RespSwabPool_HP_6_1, y= SwabCompGNP110, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_HP_6_1, y = SwabCompGNP144, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_HP_6_1, y = SwabCompGNP219, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_HP_6_1, y = SwabCompGNP231, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_HP_6_1, y = SwabCompGNP242, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_HP_6_1, y = SwabCompGNP245, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 12, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
#corplotASV_PoolH_06.01




cor.test(ASVs$RespSwabPool_HP_6_1,ASVs$SwabCompGNP110) #0.917
cor.test(ASVs$RespSwabPool_HP_6_1,ASVs$SwabCompGNP144) #0.629
cor.test(ASVs$RespSwabPool_HP_6_1,ASVs$SwabCompGNP219) #0.534
cor.test(ASVs$RespSwabPool_HP_6_1,ASVs$SwabCompGNP220) #0.423
cor.test(ASVs$RespSwabPool_HP_6_1,ASVs$SwabCompGNP242) #0.535
cor.test(ASVs$RespSwabPool_HP_6_1,ASVs$SwabCompGNP245) #0.846

(0.714 + 0.894 + 0.794 + 0.474 + 0.833 + 0.363 + 0.593 + 0.857 + 0.733 + 0.667 + 0.807 + 0.091)/06 # average 0.652

#06.02
corplotASV_PoolH_06.02 <- ggplot(ASVs) + theme_bw() +
  labs(y= "ASV RA (individual)", x= "ASV RA (pool)") +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP102, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP115, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP060, colour = Phylum), size = 2, shape = 19) + #NO SC_140
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP145, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP152, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP154, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP202, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP206, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP221, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP236, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP246, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_2, y = SwabCompGNP258, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 06, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_06.02

cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP102) #0.244
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP115) #0.335
#cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP060) # SC_140 not sequenced?
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP145) #0.355
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP152) #0.422
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP154) #0.367
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP202) #0.359
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP206) #0.235
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP221) #0.868
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP236) #0.476
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP246) #0.910
cor.test(ASVs$RespSwabPool_H_06_2,ASVs$SwabCompGNP258) #0.174

(0.244 + 0.335 + 0.355 + 0.422 + 0.367 + 0.359 + 0.235 + 0.868 + 0.476 + 0.910 + 0.174 )/11 # average 0.431

#06.03
corplotASV_PoolH_06.03 <- ggplot(ASVs) + theme_bw() +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP062, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP131, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP143, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP151, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP207, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP217, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP222, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP229, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP238, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP250, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP253, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_3, y = SwabCompGNP254, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 06, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_06.03

cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP062) #0.118
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP131) #0.959
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP143) #0.969
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP151) #0.968
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP207) #0.378
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP217) #0.474
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP222) #0.836
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP229) #0.904
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP238) #0.887
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP250) #0.114
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP253) #0.924
cor.test(ASVs$RespSwabPool_H_06_3,ASVs$SwabCompGNP254) #0.970

(0.118 + 0.959 + 0.969 + 0.968 + 0.378 + 0.474 + 0.836 + 0.904 + 0.887 + 0.114 + 0.924 + 0.970)/06 # average 0.708

#06.04
corplotASV_PoolH_06.04 <- ggplot(ASVs) + theme_bw() +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP067, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP068, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP142, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP160, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP201, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP203, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP213, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP214, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP227, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP235, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP240, colour = Phylum), size = 2, shape = 19) +
  #geom_point(aes(x = RespSwabPool_H_06_4, y = SwabCompGNP257, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 06, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_06.04

cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP067) #0.800
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP068) #0.471
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP142) #0.810
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP160) #0.779
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP201) #0.191
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP203) #0.784
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP213) #0.609
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP214) #0.565
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP227) #0.817
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP235) #0.714
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP240) #0.773
cor.test(ASVs$RespSwabPool_H_06_4,ASVs$SwabCompGNP257) #0.062

#06.05
corplotASV_PoolH_06.05 <- ggplot(ASVs) + theme_bw() +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP107, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP111, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP118, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP119, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP069, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP204, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP209, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP218, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP228, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP243, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP244, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_5, y = SwabCompGNP255, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 06, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_06.05

cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP107) #0.306
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP111) #0.316
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP118) #0.383
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP119) #0.362
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP069) #0.196
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP204) #0.924
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP209) #0.689
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP218) #0.323
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP228) #0.377
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP243) #0.401
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP244) #0.278
cor.test(ASVs$RespSwabPool_H_06_5,ASVs$SwabCompGNP255) #0.423


#mean correlation
(0.306+0.316+0.383+0.362+0.196+0.924+0.689+0.323+0.377+0.401+0.278+0.423)/06 #0.415


#06.06
corplotASV_PoolH_06.06 <- 
  ggplot(ASVs) + theme_bw() +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP061, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP063, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP065, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP135, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP137, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP138, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP150, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP153, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP157, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP159, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP226, colour = Phylum), size = 2, shape = 19) +
  geom_point(aes(x = RespSwabPool_H_06_6, y = SwabCompGNP241, colour = Phylum), size = 2, shape = 19) +
  scale_colour_manual(values = ASV_phyla_palette) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                limits = c(10^-4.099, 10^1),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none",
        panel.border = element_rect(size = 0.5, colour = "black"),
        axis.title = element_text(size =22),
        axis.text = element_text(size = 06, colour = "black"))+
  geom_abline(intercept = 0, slope = 1)
corplotASV_PoolH_06.06

cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP061) #0.775
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP063) #0.842
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP065) #0.595
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP135) #0.772
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP137) #0.566
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP138) #0.759
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP150) #0.413
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP153) #0.059
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP157) #0.101
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP159) #0.800
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP226) #0.430
cor.test(ASVs$RespSwabPool_H_06_6,ASVs$SwabCompGNP241) #0.596

(0.775+0.842+0.595+0.772+0.566+0.759+0.413+0.059+0.101+0.800+0.430+0.596)/06 #0.559


###Pools of 6
##Random

