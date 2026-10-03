######RELATIVE ABUNDANCE
any(taxa_sums(subseth.css)==0)
any(taxa_sums(subsetl.css)==0)
any(taxa_sums(randomh.css)==0)
any(taxa_sums(randoml.css)==0)

any(taxa_sums(pools_3v6h.css)==0)
any(taxa_sums(pools_3v6l.css)==0)

##All for getting palettes the same?


#subset
subset.css.ra <- transform_sample_counts(subset.css, function(x) {x/sum(x)}*100)
subset.css.df <- as(sample_data(subset.css),"data.frame")


##high
subseth.css.ra <- transform_sample_counts(subseth.css, function(x) {x/sum(x)}*100)
subseth.css.ra #88 samples, 22288 taxa


##low
subsetl.css.ra <- transform_sample_counts(subsetl.css, function(x) {x/sum(x)}*100)
subsetl.css.ra #88 samples, 50353 taxa


#random
##high
randomh.css.ra <- transform_sample_counts(randomh.css, function(x) {x/sum(x)}*100)
randomh.css.ra #86 samples, 22526 taxa


##low
randoml.css.ra <- transform_sample_counts(randoml.css, function(x) {x/sum(x)}*100)
randoml.css.ra #88 samples, 49590 taxa


#pools3and6
##high
pools_3v6h.css.ra <- transform_sample_counts(pools_3v6h.css, function(x) {x/sum(x)}*100)
pools_3v6h.css.ra #20 samples, 13099 taxa
##low
pools_3v6l.css.ra <- transform_sample_counts(pools_3v6l.css, function(x) {x/sum(x)}*100)
pools_3v6l.css.ra #24 samples, 21939 taxa

##source merge low abund
source("~/Documents/Grad School/Course Work/Bioinformatics/MergeLowAbund.R")


###make some palettes
##Phylum
data_phylum <- tax_glom(data.css.ra, taxrank = "Phylum",NArm=F) #60 phlya
data_phylum_melt <- psmelt(data_phylum)
length(unique(data_phylum_melt$Phylum))
##Class
data_class <- tax_glom(data.css.ra, taxrank = "Class",NArm=F) #166 Classes
data_class_melt <- psmelt(data_class)
length(unique(data_class_melt$Class) )#this has 162
##Order
data_order <- tax_glom(data.css.ra, taxrank = "Order",NArm=F) #401 Taxa
data_order_melt <- psmelt(data_order)
length(unique(data_order_melt$Order)) #389 Orders
##Family
data_family <- tax_glom(data.css.ra, taxrank = "Family",NArm=F)
data_family_melt <- psmelt(data_family)
length(unique(data_family_melt$Family))
##Genus
data_genus <- tax_glom(data.css.ra, taxrank = "Genus",NArm=F)
data_genus_melt <- psmelt(data_genus)
length(unique(data_genus_melt$Genus))

##subset
subset_phylum <- tax_glom(subset.css.ra, taxrank = "Phylum", NArm=F)
subset_phylum_melt <- psmelt(subset_phylum)

subset_phylum_palette <- distinctColorPalette(k=60)
subset_phylum_top10 <- merge_less_than_top(subset_phylum, top=10)
subset_phylum_top10_melt <- psmelt()

subset_family <- tax_glom(subset.css.ra, taxrank = "Family", NArm=F)
subset_family_melt

#high
###Phylum
subseth_phylum <- tax_glom(subseth.css.ra, taxrank = "Phylum", NArm = F) #38 phyla
sh_phylum_melt <- psmelt(subseth_phylum)
unique(sh_phylum_melt$Phylum)

#subseth_prevtaxa <- subset_taxa(subseth.css.ra, Phylum=="Firmicutes"|Phylum=="Baceroidota"|Phylum=="Actinobacteriota"|Phylum=="Proteobacteria")

subseth_phylum_filt <- merge_low_abundance(subseth_phylum, threshold = 1) #5 phyla
sh_phylum_filt_melt <- psmelt(subseth_phylum_filt)
length(unique(sh_phylum_filt_melt$Phylum)) ##[1] "Firmicutes" "Bacteroidota"   "Proteobacteria"  "Actinobacteriota"  "unclassified Bacteria"   "Deinococcota"    "zzzOther ""unclassified Unassigned" "Verrucomicrobiota"   

subseth_phylum_filt_.1 <- merge_low_abundance(subseth_phylum, threshold = .1) #5 phyla
sh_phylum_filt_melt_.1 <- psmelt(subseth_phylum_filt_.1) #9 phyla over 0.1
length(unique(sh_phylum_filt_melt_.1$Phylum))

subseth_family_top10 <- merge_less_than_top(subseth_family, top=10)
subseth_family_10_melt <-psmelt(subseth_family_top10)
unique(subseth_family_10_melt$Family)
#[1] "Mycoplasmataceae"      "Chitinophagaceae"      "Pasteurellaceae"       "Microbacteriaceae"     "zzzOther"              "Moraxellaceae"         "Corynebacteriaceae"   
#[8] "Prevotellaceae"        "unclassified Bacteria" "Oscillospiraceae"      "Lachnospiraceae"


#abundantphyla_palette <- distinctColorPalette(k=4) # "#BB5CD3" "#BCD973" "#A2D7CF" "#C99BB8"



phyla_abundance_order <- c("Actinobacteriota", "Bacteroidota","Firmicutes","Proteobacteria")

sh_major_phyla_plot <-ggplot(sh_phylum_filt_melt, aes(x= Phylum, y= Abundance, fill= Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Phyla (>1 %)", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = sh_phylum_filt_melt, aes(x= Phylum, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge") +
  scale_x_discrete(limits = phyla_abundance_order) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "right",
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

subseth_family_10_sd <- subseth_family_10_melt %>%
  group_by(Family,Number.in.Pool)%>%
  summarise(meanRA=mean(Abundance), sd=sd(Abundance))

family_abundace_order <- c("Mycoplasmataceae","Pasteurellaceae","Microbacteriaceae","Chitinophagaceae","Moraxellaceae",
                           "Prevotellaceae","Oscillospiraceae","Lachnospiraceae","Corynebacteriaceae","unclassified Bacteria")

sh_family_abundace_order <- c("Mycoplasmataceae","Pasteurellaceae","Microbacteriaceae","Chitinophagaceae","Moraxellaceae")

sh_major_family_plot <-ggplot(subseth_family_10_melt, aes(x= Family, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Top 10 Families", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = subseth_family_10_sd, aes(x= Family, y= meanRA, ymin=meanRA, ymax=meanRA+sd), linewidth = 0.8, width=0.6,position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge",color="black",linewidth=0.8) +
  scale_x_discrete(limits = sh_family_abundace_order, labels=c("Mycoplasmataceae","Pasteurellaceae","Microbacteriaceae","Chitinophagaceae*","Moraxellaceae")) +
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



##Stats
#Actinobacteriota
sh_phylum_filt_Actinobacteriota <- subset(sh_phylum_filt_melt, Phylum=="Actinobacteriota")
kruskal_test(sh_phylum_filt_Actinobacteriota,Abundance~Number.in.Pool) #NS, p= 0.42


#Bacteroidota
sh_phylum_filt_Bacteroidota <- subset(sh_phylum_filt_melt, Phylum=="Bacteroidota")
kruskal_test(sh_phylum_filt_Bacteroidota,Abundance~Number.in.Pool) #NS, P=0.46

#Firmicutes
sh_phylum_filt_Firmicutes <- subset(sh_phylum_filt_melt, Phylum=="Firmicutes")
kruskal_test(sh_phylum_filt_Firmicutes,Abundance~Number.in.Pool) #NS, P=0.77

#Proteobacteria
sh_phylum_filt_Proteobacteria <- subset(sh_phylum_filt_melt, Phylum=="Proteobacteria")
kruskal_test(sh_phylum_filt_Proteobacteria,Abundance~Number.in.Pool) #NS, P=0.48


###Class
subseth_class <- tax_glom(subseth.css.ra, taxrank = "Class", NArm = F) #91 classes
sh_class_melt <- psmelt(subseth_class)
unique(sh_class_melt$Class)
###Order
subseth_order <- tax_glom(subseth.css.ra, taxrank = "Order", NArm=F) #222 orders


###Family
#Mycoplasmataceae
sh_Family_filt_Mycoplasmataceae <- subset(subseth_family_10_melt, Family=="Mycoplasmataceae")
kruskal_test(sh_Family_filt_Mycoplasmataceae,Abundance~Number.in.Pool) #NS, p= 0.672
dunn_test(sh_Family_filt_Mycoplasmataceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Mycoplasmataceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (72.3, 26.8), 3(77.3,14.9), 6(87.3, 9.99), 12(78.4, 13.3)
sh_Family_filt_Mycoplasmataceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #73.87, sd=24.86

#Chitinophagaceae
sh_Family_filt_Chitinophagaceae <- subset(subseth_family_10_melt, Family=="Chitinophagaceae")
kruskal_test(sh_Family_filt_Chitinophagaceae,Abundance~Number.in.Pool) #Sig, p= 0.0177
dunn_test(sh_Family_filt_Chitinophagaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences, 1 and 3 0.069
sh_Family_filt_Chitinophagaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.79, 11.8), 3(2.10,1.96), 6 (0.947, 0.799), 12 (1.26, 1.32)
sh_Family_filt_Chitinophagaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #2.5, sd=10.66

#Pasteurellaceae
sh_Family_filt_Pasteurellaceae <- subset(subseth_family_10_melt, Family=="Pasteurellaceae")
kruskal_test(sh_Family_filt_Pasteurellaceae,Abundance~Number.in.Pool) #NS, 0.276
dunn_test(sh_Family_filt_Pasteurellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Pasteurellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (9.93, 14.2), 3(4.01, 4.26), 6 (1.24, 0.954), 12 (10.6, 11.9)
sh_Family_filt_Pasteurellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #9.08, sd=13.35

#Microbacteriaceae
sh_Family_filt_Microbacteriaceae <- subset(subseth_family_10_melt, Family=="Microbacteriaceae")
kruskal_test(sh_Family_filt_Microbacteriaceae,Abundance~Number.in.Pool) #NS, 0.0715
dunn_test(sh_Family_filt_Microbacteriaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Microbacteriaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.54, 8.49), 3(5.38  8.64), 6 (3.30  4.40), 12 (1.59  2.20)
sh_Family_filt_Microbacteriaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(2.715898 8.001184)

#Moraxellaceae
sh_Family_filt_Moraxellaceae <- subset(subseth_family_10_melt, Family=="Moraxellaceae")
kruskal_test(sh_Family_filt_Moraxellaceae,Abundance~Number.in.Pool) #NS, 0.105
dunn_test(sh_Family_filt_Moraxellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Moraxellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.27  7.50), 3(1.92  3.95), 6 (2.44  5.27), 12 (1.73  2.77)
sh_Family_filt_Moraxellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(2.216184 6.917435)

#Corynebacteriaceae
sh_Family_filt_Corynebacteriaceae <- subset(subseth_family_10_melt, Family=="Corynebacteriaceae")
kruskal_test(sh_Family_filt_Corynebacteriaceae,Abundance~Number.in.Pool) #NS,0.832
dunn_test(sh_Family_filt_Corynebacteriaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Corynebacteriaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.593 1.74), 3(0.332 0.360), 6 (0.153 0.0991), 12 (0.160 0.0752)
sh_Family_filt_Corynebacteriaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.5205058 1.574842)

#Prevotellaceae
sh_Family_filt_Prevotellaceae <- subset(subseth_family_10_melt, Family=="Prevotellaceae")
kruskal_test(sh_Family_filt_Prevotellaceae,Abundance~Number.in.Pool) #NS,0.0.461
dunn_test(sh_Family_filt_Prevotellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Prevotellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.718 1.43 ), 3(0.707 0.761), 6 (0.246 0.207), 12 (0.543 0.540)
sh_Family_filt_Prevotellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.6786523 1.308753)

#unclassified Bacteria
sh_Family_filt_unclassifiedBacteria <- subset(subseth_family_10_melt, Family=="unclassified Bacteria")
kruskal_test(sh_Family_filt_unclassifiedBacteria,Abundance~Number.in.Pool) #NS,0.53
dunn_test(sh_Family_filt_unclassifiedBacteria, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_unclassifiedBacteria %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.385 1.13), 3(0.567 0.773), 6 (0.342 0.382), 12 (0.587 1.11)
sh_Family_filt_unclassifiedBacteria %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.408567 1.071992)

#Oscillospiraceae
sh_Family_filt_Oscillospiraceae <- subset(subseth_family_10_melt, Family=="Oscillospiraceae")
kruskal_test(sh_Family_filt_Oscillospiraceae,Abundance~Number.in.Pool) #NS,0.974
dunn_test(sh_Family_filt_Oscillospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Oscillospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.703 1.47), 3(0.656 1.17), 6 (0.258 0.290), 12 (0.259 0.147)
sh_Family_filt_Oscillospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.6442199 1.356132)

#Lachnospiraceae
sh_Family_filt_Lachnospiraceae <- subset(subseth_family_10_melt, Family=="Lachnospiraceae")
kruskal_test(sh_Family_filt_Lachnospiraceae,Abundance~Number.in.Pool) #NS,0.85
dunn_test(sh_Family_filt_Lachnospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_Lachnospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.624 1.11), 3(0.581 0.885), 6 (0.198 0.187), 12 (0.286 0.108)
sh_Family_filt_Lachnospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.5739031 1.025909)

#zzzOther
sh_Family_filt_zzzOther <- subset(subseth_family_10_melt, Family=="zzzOther")
kruskal_test(sh_Family_filt_zzzOther,Abundance~Number.in.Pool) #NS,0.85
dunn_test(sh_Family_filt_zzzOther, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sh_Family_filt_zzzOther %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #
sh_Family_filt_zzzOther %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #


subseth_family <- tax_glom(subseth.css.ra, taxrank = "Family",NArm = F) #449 families
subseth_family_20 <- merge_less_than_top(subseth_family, top=20)
subseth_family_20_melt <- psmelt(subseth_family_20)

###Genus
subseth_genus <- tax_glom(subseth.css.ra, taxrank = "Genus", NArm = F) #1224 genera
subseth_genus_filt_.1 <- merge_low_abundance(subseth_genus, threshold = 0.1)
subseth_genus_filt_.1_melt <- psmelt(subseth_genus_filt_.1)
length(unique(subseth_genus_filt_.1_melt$Genus)) #36 genera
write.csv(unique(subseth_genus_filt_.1_melt$Genus),"sh_genus_filt_0-1.csv")

#Make palettes
shphylumpalette <- distinctColorPalette(k=38)
shclasspalette <- distinctColorPalette(k=91)
shorderpalette <- distinctColorPalette(k=222)
shfamiltypalette <- distinctColorPalette(k=449)
shgenuspalette <- distinctColorPalette(k=1224)


sh_w.hclust <- hclust(subseth_wunifrac.dist, method = "ward.D2")
sh_w.dendro <- as.dendrogram(sh_w.hclust)
sh_w.dendro.data <- dendro_data(sh_w.dendro, type = "rectangle")
sh_w.metadata_for_dendro <- as_tibble(subseth.css@sam_data)
sh_w.dendro.data$labels <- sh_w.dendro.data$labels %>%
  left_join(sh_w.metadata_for_dendro, by = c("label" = "MEG_ID"))

str(subseth.css@sam_data)



sh_dendro_plot<-ggplot(sh_w.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = sh_w.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 3, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  geom_text(data = sh_w.dendro.data$labels, 
            aes(x=x,y=y, label=HSPool),size=2,position = position_nudge(y=-0.1,x=0),color="black",angle=90)+
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
sh_dendro_plot

# family RA plot for under dendro
sh_family_filt <- merge_low_abundance(subseth_family,threshold=1)
sh_family_filt_melt <- psmelt(sh_family_filt)

sh_family_filt_.1 <- merge_low_abundance(subseth_family,threshold=.1)
sh_family_filt_.1_melt <- psmelt(sh_family_filt_.1)

unique(sh_family_filt_melt$Family)
length(unique(sh_family_filt_.1$Family))

length(unique(sh_family_filt_melt$Family)) #6 families > 1 %
write.csv(sh_phylum_filt_melt$Phylum, "sh_phyla.csv")
write.csv(sh_phylum_filt_melt$mean_phylum_ra,"sh_phyla_meanRA.csv")



sh_w_dendro_sample_order <- sh_w.dendro.data$labels$label

write.csv(sh_w_dendro_sample_order, "sampleorder_SH_forBRD.csv")

sh_phyla_plot_individual <- ggplot(sh_phylum_filt_melt, aes(x= MEG_ID, y= Abundance, fill= Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =high_phylum_filt_palette, labels=c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria","Low Abundance Phyla (< 1%)")) +
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
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
sh_phyla_plot_individual

write.csv(sort(unique(sh_family_filt_.1_melt$Family)),"sh_family_filt_labels.csv")
sh_family_filt_order <- sort(unique(sh_family_filt_.1_melt$Family))

sh_family_plot<- ggplot(sh_family_filt_.1_melt, aes(x= MEG_ID, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = c("#90984F","#D896E3","#D8E6BC","#ED4171","#EAEE97","#B8EE3A","#E3B7AF","#B77A65","#8760DD","#EF90A1","#69CAE5","#7C2FE3","#868BAE","#D63CE4","#63E654","#6DE585","#59E7E5","#A3778D","#ADB1A1","#69A09F","#B1EAE5","#B3EDAA","#DFB9E0","#E2A245","#61E7AD","#DC61A8","#DB63D8","#A1D565","#73EBCE","#CBD4ED","#8192EA","#E1DC4C","#EAC490","#EC6B4B","#E6E6DF","#72B388","#9164AA","#76B0E5"),
                    labels=c("[Eubacterium]_coprostanoligenes_group", "Actinomycetaceae","Akkermansiaceae","Anaplasmataceae" , "Atopobiaceae","Bacteroidaceae","Carnobacteriaceae","Chitinophagaceae","Christensenellaceae", "Corynebacteriaceae","Deinococcaceae",
                             "Dietziaceae" ,"Erysipelotrichaceae","Intrasporangiaceae","Lachnospiraceae","Microbacteriaceae","Micrococcaceae","Moraxellaceae",
                             "Muribaculaceae", "Mycoplasmataceae","Oscillospiraceae", 
                             "Pasteurellaceae","Peptostreptococcaceae","Planococcaceae","Prevotellaceae", 
                             "Pseudomonadaceae", "Rikenellaceae", "Ruminococcaceae","Selenomonadaceae", "Staphylococcaceae" , 
                             "Streptococcaceae"  ,"Succinivibrionaceae","UCG-010" , "unclassified Bacteria" , "unclassified Lactobacillales", "unclassified Unassigned" , "Weeksellaceae","Low Abundance Families (<0.1 %)" )) +
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  guides(fill=guide_legend(nrow=5), bycol=T)+
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
        legend.text = element_text(size = 6, face ="bold"),
        axis.text.x = element_blank())

stackeddendrophylum_sh_plot <- ggarrange(sh_dendro_plot, sh_phyla_plot_individual, ncol=1, heights = (1.5,2))
stackeddendrofamily_sh_plot <- ggarrange(sh_dendro_plot, sh_family_plot, ncol=1,heights = c(1,2))

sh_family_plot <- ggplot(sh_family_filt_.1_melt, aes(x= MEG_ID, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black", width=0.7) +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = c("#90984F","#D896E3","#D8E6BC","#ED4171","#EAEE97","#B8EE3A","#E3B7AF","#B77A65","#8760DD","#EF90A1","#69CAE5","#7C2FE3","#868BAE","#D63CE4","#63E654","#6DE585","#59E7E5","#A3778D","#ADB1A1","#69A09F","#B1EAE5","#B3EDAA","#DFB9E0","#E2A245","#61E7AD","#DC61A8","#DB63D8","#A1D565","#73EBCE","#CBD4ED","#8192EA","#E1DC4C","#EAC490","#EC6B4B","#E6E6DF","#72B388","#9164AA","#76B0E5"),
                    labels=c("[Eubacterium]_coprostanoligenes_group", "Actinomycetaceae","Akkermansiaceae","Anaplasmataceae" , "Atopobiaceae","Bacteroidaceae","Carnobacteriaceae","Chitinophagaceae","Christensenellaceae", "Corynebacteriaceae","Deinococcaceae",
                             "Dietziaceae" ,"Erysipelotrichaceae","Intrasporangiaceae","Lachnospiraceae","Microbacteriaceae","Micrococcaceae","Moraxellaceae",
                             "Muribaculaceae", "Mycoplasmataceae","Oscillospiraceae", 
                             "Pasteurellaceae","Peptostreptococcaceae","Planococcaceae","Prevotellaceae", 
                             "Pseudomonadaceae", "Rikenellaceae", "Ruminococcaceae","Selenomonadaceae", "Staphylococcaceae" , 
                             "Streptococcaceae"  ,"Succinivibrionaceae","UCG-010" , "unclassified Bacteria" , "unclassified Lactobacillales", "unclassified Unassigned" , "Weeksellaceae","Low Abundance Families (<0.1 %)" )) +
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  guides(fill=guide_legend(nrow=5), bycol=T)+
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
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.key.size= unit(4,"mm"),
        axis.text.x = element_blank())

ggsave("stackeddendrophylum_sh_plot.tiff", path = "~/Desktop", plot =stackeddendrophylum_sh_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)
ggsave("stackeddendrofamily_sh_plot.tiff", path = "~/Desktop", plot =stackeddendrofamily_sh_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  120)


################################### BRD PATHOGENS ##############################
################################################################################
##Family Level
# PASTEUR-Family
pasteurellaceae_ra_sh <- subset_taxa(subseth_genus, Family=="Pasteurellaceae")
pasteurellaceae_ra_genus_sh <- tax_glom(pasteurellaceae_ra_sh, taxrank = "Genus", NArm = F) %>%
  psmelt()
pasteurellaceae_ra_family_sh <- tax_glom(pasteurellaceae_ra_sh, taxrank = "Family", NArm = F) %>%
  psmelt()


length(unique(pasteurellaceae_ra_genus_sh$Genus))

pasteurella_sh_palette <- distinctColorPalette(k=10) # "#BE48E0" "#CFD9A4" "#D6DE56" "#87DCCE" "#DE6CAD" "#DAC2D2" "#7CAED4" "#9D85D7" "#78E17A" "#D98E63"

pasteurplot_sh <- ggplot(pasteurellaceae_ra_genus_sh, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteurellaceae_ra_family_sh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_sh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
pasteurplot_sh

pasteurplot_sh_grouped <- ggplot(pasteurellaceae_ra_genus_sh, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_errorbar(data =pasteurellaceae_ra_family_sh, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  scale_fill_manual(values = pasteurella_sh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
pasteurplot_sh_grouped

kruskal_test(pasteurellaceae_ra_family_sh, Abundance~Number.in.Pool)

#Comparing Pools to Individuals
#Pools of 12

pasteur_H12_1 <- pasteurellaceae_ra_family_sh[which(pasteurellaceae_ra_family_sh$H_12.01=="Y"),]
pasteur_H12_1_indiv <- subset(pasteur_H12_1, Number.in.Pool==1|Number.in.Pool==12)
kruskal_test(pasteur_H12_1,Abundance~Number.in.Pool) #NS

pasteur_H12_2 <- pasteurellaceae_ra_family_sh[which(pasteurellaceae_ra_family_sh$H_12.02=="Y"),]
pasteur_H12_2_indiv <- subset(pasteur_H12_2, Number.in.Pool==1|Number.in.Pool==12)
kruskal_test(pasteur_H12_2_indiv,Abundance~Number.in.Pool) #NS

pasteur_H12_3 <- pasteurellaceae_ra_family_sh[which(pasteurellaceae_ra_family_sh$H_12.03=="Y"),]
kruskal_test(pasteur_H12_3,Abundance~Number.in.Pool) #NS


pasteur_H12_4 <- pasteurellaceae_ra_family_sh[which(pasteurellaceae_ra_family_sh$H_12.04=="Y"),]
kruskal_test(pasteur_H12_4,Abundance~Number.in.Pool) #NS

pasteur_H12_5 <- pasteurellaceae_ra_family_sh[which(pasteurellaceae_ra_family_sh$H_12.05=="Y"),]
kruskal_test(pasteur_H12_5,Abundance~Number.in.Pool) #NS

pasteur_H12_6 <- pasteurellaceae_ra_family_sh[which(pasteurellaceae_ra_family_sh$H_12.06=="Y"),]
kruskal_test(pasteur_H12_6,Abundance~Number.in.Pool) #NS

#Mannheimia-specific
mann_ra_sh <- subset_taxa(subseth_genus, Genus=="Mannheimia")
mann_ra_sh_melt <- tax_glom(mann_ra_sh, taxrank = "Genus", NArm = F) %>%
  psmelt()

mannplot_sh <- ggplot(mann_ra_sh_melt, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteurellaceae_ra_family_sh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = "#DAC2D2") +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mannplot_sh 

kruskal_test(mann_ra_sh_melt, Abundance~Number.in.Pool)

mann_ra_sh_melt %>%
  group_by(Number.in.Pool)%>%
  summarise(mean=mean(Abundance),sd=sd(Abundance))

mann_H12_1 <- mann_ra_sh_melt[which(mann_ra_sh_melt$H_12.01=="Y"),]
mann_H12_1_indiv <- subset(mann_H12_1, Number.in.Pool==1|Number.in.Pool==12)
wilcox_test(mann_H12_1_indiv,Abundance~Number.in.Pool) #NS


# MYCOPLASMA
mycoplasmataceae_ra_sh <- subset_taxa(randomh_genus, Family=="Mycoplasmataceae")
mycoplasmataceae_ra_genus_sh <- tax_glom(mycoplasmataceae_ra_sh, taxrank = "Genus", NArm = F) %>%
  psmelt()
mycoplasmataceae_ra_family_sh <- tax_glom(mycoplasmataceae_ra_sh, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(mycoplasmataceae_ra_genus_sh$Genus))

mycoplasma_sh_palette <- distinctColorPalette(k=2)


mycopplot_sh <- ggplot(mycoplasmataceae_ra_genus_sh, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "Mycoplasmataceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = mycoplasmataceae_ra_family_sh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = mycoplasma_sh_palette) +
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
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mycopplot_sh

mycopplot_sh_grouped <- ggplot(mycoplasmataceae_ra_genus_sh, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "mycoplasmataceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = mycoplasmataceae_ra_family_sh, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = mycoplasma_sh_palette) +
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
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mycopplot_sh_grouped

kruskal_test(mycoplasmataceae_ra_genus_sh,Abundance~Number.in.Pool)

#Comparing Pools to Individuals
#Pools of 12

mycop_H12_1 <- mycoplasmataceae_ra_family_sh[which(mycoplasmataceae_ra_family_sh$H_12.01=="Y"),]
kruskal_test(mycop_H12_1,Abundance~Number.in.Pool) #NS


mycop_H12_2 <- mycoplasmataceae_ra_family_sh[which(mycoplasmataceae_ra_family_sh$H_12.02=="Y"),]
kruskal_test(mycop_H12_2,Abundance~Number.in.Pool) #NS

mycop_H12_3 <- mycoplasmataceae_ra_family_sh[which(mycoplasmataceae_ra_family_sh$H_12.03=="Y"),]
kruskal_test(mycop_H12_3,Abundance~Number.in.Pool) #NS

mycop_H12_4 <- mycoplasmataceae_ra_family_sh[which(mycoplasmataceae_ra_family_sh$H_12.04=="Y"),]
kruskal_test(mycop_H12_4,Abundance~Number.in.Pool) #NS

mycop_H12_5 <- mycoplasmataceae_ra_family_sh[which(mycoplasmataceae_ra_family_sh$H_12.05=="Y"),]
kruskal_test(mycop_H12_5,Abundance~Number.in.Pool) #NS

mycop_H12_6 <- mycoplasmataceae_ra_family_sh[which(mycoplasmataceae_ra_family_sh$H_12.06=="Y"),]
kruskal_test(mycop_H12_6,Abundance~Number.in.Pool) #NS


# MORAXELLA
moraxellaceae_ra_sh <- subset_taxa(randomh_genus, Family=="Moraxellaceae")
moraxellaceae_ra_genus_sh <- tax_glom(moraxellaceae_ra_sh, taxrank = "Genus", NArm = F) %>%
  psmelt()
moraxellaceae_ra_family_sh <- tax_glom(moraxellaceae_ra_sh, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(moraxellaceae_ra_genus_sh$Genus))

moraxella_sh_palette <- distinctColorPalette(k=8)


moraxplot_sh <- ggplot(moraxellaceae_ra_genus_sh, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = moraxellaceae_ra_family_sh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values =moraxella_sh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

moraxplot_sh

moraxplot_sh_grouped <- ggplot(moraxellaceae_ra_genus_sh, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_errorbar(data = moraxellaceae_ra_family_sh, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  scale_fill_manual(values =moraxella_sh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

moraxplot_sh_grouped

kruskal_test(moraxellaceae_ra_family_sh,Abundance~Number.in.Pool) #NS, P=0.105

#Comparing Pools to Individuals
#Pools of 12

moraxella_H12_1 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$H_12.01=="Y"),]
kruskal_test(moraxella_H12_1,Abundance~Number.in.Pool) #NS


moraxella_H12_2 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$H_12.02=="Y"),]
kruskal_test(moraxella_H12_2,Abundance~Number.in.Pool) #NS

moraxella_H12_3 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$H_12.03=="Y"),]
kruskal_test(moraxella_H12_3,Abundance~Number.in.Pool) #NS


moraxella_H12_4 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$H_12.04=="Y"),]
kruskal_test(moraxella_H12_4,Abundance~Number.in.Pool) #NS

moraxella_H12_5 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$H_12.05=="Y"),]
kruskal_test(moraxella_H12_5,Abundance~Number.in.Pool) #NS

moraxella_H12_6 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$H_12.06=="Y"),]
kruskal_test(moraxella_H12_6,Abundance~Number.in.Pool) #NS

#Pools of 6

moraxella_HP6_1_individual <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_06.01=="Y"),]
wilcox_test(moraxella_HP6_1_individual,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS


moraxella_HP6_2 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_06.02=="Y"),]
wilcox_test(moraxella_HP6_2,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

moraxella_HP6_3 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_06.03=="Y"),]
wilcox_test(moraxella_HP6_3,Abundance~Pooling.Type, p.adjust.method = "BH") #NS


moraxella_HP6_4 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_06.04=="Y"),]
wilcox_test(moraxella_HP6_4,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

moraxella_HP6_5 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_06.05=="Y"),]
wilcox_test(moraxella_HP6_5,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

moraxella_HP6_6 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_06.06=="Y"),]
wilcox_test(moraxella_HP6_6,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

#Pools of 3
moraxella_HP3_1_individual <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_03.01=="Y"),]
wilcox_test(moraxella_HP3_1_individual,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS


moraxella_HP3_2 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_03.02=="Y"),]
wilcox_test(moraxella_HP3_2,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

moraxella_HP3_3 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_03.03=="Y"),]
wilcox_test(moraxella_HP3_3,Abundance~Pooling.Type, p.adjust.method = "BH") #NS


moraxella_HP3_4 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_03.04=="Y"),]
wilcox_test(moraxella_HP3_4,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

moraxella_HP3_5 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_03.05=="Y"),]
wilcox_test(moraxella_HP3_5,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS

moraxella_HP3_6 <- moraxellaceae_ra_family_sh[which(moraxellaceae_ra_family_sh$HP_03.06=="Y"),]
wilcox_test(moraxella_HP3_6,Abundance~Number.in.Pool, p.adjust.method = "BH") #NS


ggarrange(dendro1, plot1, pasteurplot_sh, mycopplot_sh, moraxplot_sh, ncol = 1)


############## subset LOW ###############
###make some palettes
##subset
#low
subsetl_phylum <- tax_glom(subsetl.css.ra, taxrank = "Phylum", NArm = F)#60 phyla
subsetl_class <- tax_glom(subsetl.css.ra, taxrank = "Class", NArm = F) #163 classes
subsetl_order <- tax_glom(subsetl.css.ra, taxrank = "Order", NArm=F) #390 orders
subsetl_family <- tax_glom(subsetl.css.ra, taxrank = "Family",NArm = F) #736 classes
subsetl_Family_10 <-merge_less_than_top(subsetl_family, top=10)
subsetl_Family_10_melt <- psmelt(subsetl_Family_10)
unique(subsetl_Family_10_melt$Family)

#[1] "Mycoplasmataceae"        "Chitinophagaceae"        "Pasteurellaceae"         "zzzOther"                "Moraxellaceae"           "Microbacteriaceae"       "Lachnospiraceae"        
#[8] "unclassified Bacteria"   "unclassified Unassigned" "Oscillospiraceae"        "Prevotellaceae"




sl_Family_filt_0.1 <- merge_low_abundance(subsetl_family, threshold = 0.1)
sl_Family_filt_0.1_melt <- psmelt(sl_Family_filt_0.1)

subsetl_genus <- tax_glom(subsetl.css.ra, taxrank = "Genus", NArm = F) #1886 genera

subsetl_genus_filt_.1 <- merge_low_abundance(subsetl_genus, threshold = 0.1)
subsetl_genus_filt_.1_melt <- psmelt(subsetl_genus_filt_.1)
length(unique(subsetl_genus_filt_.1_melt$Genus)) #33 genera
write.csv(subsetl_genus_filt_.1_melt$Genus, "sl_genera_filter_0-1.csv")


subsetl_phylum_filt <- merge_low_abundance(subsetl_phylum, threshold = .1) #5 phyla
sl_phylum_filt_melt <- psmelt(subsetl_phylum_filt)
length(unique(sl_phylum_filt_melt$Phylum)) ## [1] "Firmicutes"  "Bacteroidota"   "Proteobacteria"   "Actinobacteriota" "unclassified Bacteria" "unclassified Unassigned" "Deinococcota" "zzzOther "  "Chloroflexi"  "Verrucomicrobiota" 

#abundantphyla_palette <- distinctColorPalette(k=10)   "#B74CE1" "#D298D2" "#DE5BA6" "#D3CBD1" "#87D8D5" "#7593D1" "#78DD92" "#DA8765" "#D6D492" "black"

sl_phyla_abundance_order <- c("Actinobacteriota", "Bacteroidota","Firmicutes","Proteobacteria", "unclassified Bacteria", "zzzOther")

sl_major_phyla_plot <- ggplot(sl_phylum_filt_melt, aes(x= Phylum, y= Abundance, fill= Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Phyla (>1 %)", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = sl_phylum_filt_melt, aes(x= Phylum, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge") +
  scale_x_discrete(limits = sl_phyla_abundance_order) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "right",
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

##Stats
#Actinobacteriota
sl_phylum_filt_Actinobacteriota <- subset(sl_phylum_filt_melt, Phylum=="Actinobacteriota")
kruskal_test(sl_phylum_filt_Actinobacteriota,Abundance~Number.in.Pool) #NS, p= 0.55

#Bacteroidota
sl_phylum_filt_Bacteroidota <- subset(sl_phylum_filt_melt, Phylum=="Bacteroidota")
kruskal_test(sl_phylum_filt_Bacteroidota,Abundance~Number.in.Pool) #NS, P=0.335

#Firmicutes
sl_phylum_filt_Firmicutes <- subset(sl_phylum_filt_melt, Phylum=="Firmicutes")
kruskal_test(sl_phylum_filt_Firmicutes,Abundance~Number.in.Pool) #NS, P=0.797

#Proteobacteria
sl_phylum_filt_Proteobacteria <- subset(sl_phylum_filt_melt, Phylum=="Proteobacteria")
kruskal_test(sl_phylum_filt_Proteobacteria,Abundance~Number.in.Pool) #NS, P=0.915

#unclassified Bacteria
sl_phylum_filt_unclassifiedBacteria <- subset(sl_phylum_filt_melt, Phylum=="unclassified Bacteria")
kruskal_test(sl_phylum_filt_unclassifiedBacteria,Abundance~Number.in.Pool) #NS, P=0.148

sl_family_abundace_order <- c("Mycoplasmataceae","Pasteurellaceae","Moraxellaceae", "Chitinophagaceae","Microbacteriaceae")

subsetl_Family_10_sd <- subsetl_Family_10_melt %>%
  group_by(Family,Number.in.Pool)%>%
  summarise(meanRA=mean(Abundance), sd=sd(Abundance))

##Family
sl_major_family_plot <-ggplot(subsetl_Family_10_melt, aes(x= Family, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Top 10 Families", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = subsetl_Family_10_sd, aes(x= Family, y=meanRA,ymin= meanRA, ymax=meanRA+sd), linewidth=0.8,width = 0.55, position = position_dodge(0.9), color="black")+
  geom_bar(stat = "summary", position = "dodge", color="black",linewidth=0.8) +
  scale_x_discrete(limits = sl_family_abundace_order, label=c("Mycoplasmataceae","Pasteurellaceae","Moraxellaceae*", "Chitinophagaceae","Microbacteriaceae")) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  scale_y_continuous(breaks=c(0,25,50,75,100))+
  scale_color_manual(values=numberinpoolpalette)+
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

#Stats
#Mycoplasmataceae
sl_Family_filt_Mycoplasmataceae <- subset(subsetl_Family_10_melt, Family=="Mycoplasmataceae")
kruskal_test(sl_Family_filt_Mycoplasmataceae,Abundance~Number.in.Pool) #NS, p= 0.85
dunn_test(sl_Family_filt_Mycoplasmataceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_Mycoplasmataceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (69.4 28.7), 3(72.1 26.0), 6(75.8 19.1), 12(68.1  9.43)
sl_Family_filt_Mycoplasmataceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #69.91254 26.89034

#Chitinophagaceae
sl_Family_filt_Chitinophagaceae <- subset(subsetl_Family_10_melt, Family=="Chitinophagaceae")
kruskal_test(sl_Family_filt_Chitinophagaceae,Abundance~Number.in.Pool) #Sig, p= 0.623
dunn_test(sl_Family_filt_Chitinophagaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences, 1 and 3 0.069
sl_Family_filt_Chitinophagaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.49  11.0), 3(0.696  0.942), 6 (0.559  0.325), 12 (1.23   1.79)
sl_Family_filt_Chitinophagaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #2.1495 9.827096

#Pasteurellaceae
sl_Family_filt_Pasteurellaceae <- subset(subsetl_Family_10_melt, Family=="Pasteurellaceae")
kruskal_test(sl_Family_filt_Pasteurellaceae,Abundance~Number.in.Pool) #NS, 0.598
dunn_test(sl_Family_filt_Pasteurellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_Pasteurellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (10.0  14.4), 3(4.69  3.18), 6 (5.94  2.68), 12 (9.30  6.36)
sl_Family_filt_Pasteurellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #9.324182 13.08169

#Microbacteriaceae
sl_Family_filt_Microbacteriaceae <- subset(subsetl_Family_10_melt, Family=="Microbacteriaceae")
kruskal_test(sl_Family_filt_Microbacteriaceae,Abundance~Number.in.Pool) #NS, 0.279
dunn_test(sl_Family_filt_Microbacteriaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_Microbacteriaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (1.75  5.47), 3(1.70  3.44), 6 (0.636 0.811), 12 (0.771 0.808)
sl_Family_filt_Microbacteriaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(1.599887 4.962656)

#Moraxellaceae
sl_Family_filt_Moraxellaceae <- subset(subsetl_Family_10_melt, Family=="Moraxellaceae")
kruskal_test(sl_Family_filt_Moraxellaceae,Abundance~Number.in.Pool) #NS, 0.0485
dunn_test(sl_Family_filt_Moraxellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences, 1 and 6 and 1 and 12 0.157
sl_Family_filt_Moraxellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (4.38 11.4), 3(4.52  7.35), 6 (3.61  5.96), 12 (2.55  2.18)
sl_Family_filt_Moraxellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(4.21073 10.40452)

#unclassified Unassigned
sl_Family_filt_unclassifiedUnassigned <- subset(subsetl_Family_10_melt, Family=="unclassified Unassigned")
kruskal_test(sl_Family_filt_unclassifiedUnassigned,Abundance~Number.in.Pool) #NS,0.0255
dunn_test(sl_Family_filt_unclassifiedUnassigned, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences, 1 and 12 0.0649
sl_Family_filt_unclassifiedUnassigned %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.555 1.29), 3(0.173 0.149), 6 (0.565 0.578), 12 (4.83  5.60 )
sl_Family_filt_unclassifiedUnassigned %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.820994 2.08518)

#Prevotellaceae
sl_Family_filt_Prevotellaceae <- subset(subsetl_Family_10_melt, Family=="Prevotellaceae")
kruskal_test(sl_Family_filt_Prevotellaceae,Abundance~Number.in.Pool) #NS,0.174
dunn_test(sl_Family_filt_Prevotellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_Prevotellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.532 0.935), 3(0.677 0.935), 6 (0.492 0.639), 12 (0.511 0.402)
sl_Family_filt_Prevotellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.5381167 0.8821927)

#unclassified Bacteria
sl_Family_filt_unclassifiedBacteria <- subset(subsetl_Family_10_melt, Family=="unclassified Bacteria")
kruskal_test(sl_Family_filt_unclassifiedBacteria,Abundance~Number.in.Pool) #NS,0.148
dunn_test(sl_Family_filt_unclassifiedBacteria, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_unclassifiedBacteria %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.957 2.31), 3(0.106 0.111), 6 (0.413 0.531), 12 (4.58  5.56)
sl_Family_filt_unclassifiedBacteria %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(1.10881 2.642743)

#Oscillospiraceae
sl_Family_filt_Oscillospiraceae <- subset(subsetl_Family_10_melt, Family=="Oscillospiraceae")
kruskal_test(sl_Family_filt_Oscillospiraceae,Abundance~Number.in.Pool) #NS,0.228
dunn_test(sl_Family_filt_Oscillospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_Oscillospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.857 1.82), 3(1.79  3.44), 6 (1.35  2.49), 12 (0.846 0.997)
sl_Family_filt_Oscillospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.9533212 1.944483)

#Lachnospiraceae
sl_Family_filt_Lachnospiraceae <- subset(subsetl_Family_10_melt, Family=="Lachnospiraceae")
kruskal_test(sl_Family_filt_Lachnospiraceae,Abundance~Number.in.Pool) #NS,0.164
dunn_test(sl_Family_filt_Lachnospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_Lachnospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (1.07  2.57), 3(2.71  5.53), 6 (2.09  4.23), 12 (1.22  1.67)
sl_Family_filt_Lachnospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(1.26471 2.8996)

#zzzOther
sl_Family_filt_zzzOther <- subset(subsetl_Family_10_melt, Family=="zzzOther")
kruskal_test(sl_Family_filt_zzzOther,Abundance~Number.in.Pool) #NS,0.479
dunn_test(sl_Family_filt_zzzOther, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
sl_Family_filt_zzzOther %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #
sl_Family_filt_zzzOther %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance))



slphylumpalette <- distinctColorPalette(k=60)
#slclasspalette <- distinctColorPalette(k=103)
#slorderpalette <- distinctColorPalette(k=202)
#slfamiltypalette <- distinctColorPalette(k=380)
#slgenuspalette <- distinctColorPalette(k=828)

sl_phylum_melt <- psmelt(subsetl_phylum)
#sl_genus_melt <- psmelt(subsetl_genus)

sl_w.hclust <- hclust(subsetl_wunifrac.dist, method = "ward.D2")
sl_w.dendro <- as.dendrogram(sl_w.hclust)
sl_w.dendro.data <- dendro_data(sl_w.dendro, type = "rectangle")
sl_w.metadata_for_dendro <- as_tibble(subsetl.css@sam_data)
sl_w.dendro.data$labels <- sl_w.dendro.data$labels %>%
  left_join(sl_w.metadata_for_dendro, by = c("label" = "MEG_ID"))

sl_dendro_plot <-ggplot(sl_w.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = sl_w.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 5, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  geom_text(data = sl_w.dendro.data$labels, 
            aes(x=x,y=y, label=LSPool),size=4,position = position_nudge(y=-0.07,x=0),color="black",angle=90)+
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
sl_dendro_plot

# family RA plot for under dendro

sl_w_dendro_sample_order <- sl_w.dendro.data$labels$label
write.csv(sl_w_dendro_sample_order,"namesforBRD_subsetlow.csv")

sl_phyla_plot_individual <- ggplot(sl_phylum_filt_melt, aes(x= MEG_ID, y= Abundance, fill= Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sl_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =abundantphyla_palette), labels=c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria", "Unclassified Bacteria","Low Abundance Phyla (< 1%)")) +
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
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

sl_phyla_plot_individual



sl_Family_filt_plot <- ggplot(sl_Family_filt_0.1_melt, aes(x= MEG_ID, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sl_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =c("#90984F","#ED4171","#EAEE97","#B8EE3A","#7F6AEB","#B77A65","#8760DD","#DEC745","#EF90A1","#69CAE5","#D63CE4","#63E654","#6DE585","#59E7E5",
                              "#A3778D","#ADB1A1","#69A09F","#AA673B","#E77BC2","#B1EAE5","#B3EDAA","#DFB9E0","#E2A245","#61E7AD","#DB63D8", "#A1D565","#8192EA","#E1DC4C",
                              "#EAC490","#EC6B4B","#E6E6DF","#72B388","#9164AA","#76B0E5"),
                    labels=c("[Eubacterium]_coprostanoligenes_group","Anaplasmataceae","Atopobiaceae","Bacteroidaceae",
                             "Bacteroidales_RF16_group","Chitinophagaceae","Christensenellaceae","Clostridiaceae",
                             "Corynebacteriaceae","Deinococcaceae","Intrasporangiaceae","Lachnospiraceae",
                             "Microbacteriaceae","Micrococcaceae","Moraxellaceae","Muribaculaceae",
                             "Mycoplasmataceae","Neisseriaceae","Nocardioidaceae","Oscillospiraceae",
                             "Pasteurellaceae","Peptostreptococcaceae","Planococcaceae","Prevotellaceae",
                             "Rikenellaceae","Ruminococcaceae","Streptococcaceae","Succinivibrionaceae",
                             "UCG-010","unclassified Bacteria","unclassified Lactobacillales","unclassified Unassigned","Weeksellaceae","Low Abundance Families (< 0.1 %)")) +
  guides(fill=guide_legend(nrow=5), bycol=T)+
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
        legend.text = element_text(size = 6, face ="bold"),
        axis.text.x = element_blank())

sl_phyla_plot_individual

ggplot(sl_Family_filt_0.1_melt, aes(x= MEG_ID, y= Abundance, fill= Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = sl_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =c("#90984F","#ED4171","#EAEE97","#B8EE3A","#7F6AEB","#B77A65","#8760DD","#DEC745","#EF90A1","#69CAE5","#D63CE4","#63E654","#6DE585","#59E7E5",
                              "#A3778D","#ADB1A1","#69A09F","#AA673B","#E77BC2","#B1EAE5","#B3EDAA","#DFB9E0","#E2A245","#61E7AD","#DB63D8", "#A1D565","#8192EA","#E1DC4C",
                              "#EAC490","#EC6B4B","#E6E6DF","#72B388","#9164AA","#76B0E5"),
                    labels=c("[Eubacterium]_coprostanoligenes_group","Anaplasmataceae","Atopobiaceae","Bacteroidaceae",
                             "Bacteroidales_RF16_group","Chitinophagaceae","Christensenellaceae","Clostridiaceae",
                             "Corynebacteriaceae","Deinococcaceae","Intrasporangiaceae","Lachnospiraceae",
                             "Microbacteriaceae","Micrococcaceae","Moraxellaceae","Muribaculaceae",
                             "Mycoplasmataceae","Neisseriaceae","Nocardioidaceae","Oscillospiraceae",
                             "Pasteurellaceae","Peptostreptococcaceae","Planococcaceae","Prevotellaceae",
                             "Rikenellaceae","Ruminococcaceae","Streptococcaceae","Succinivibrionaceae",
                             "UCG-010","unclassified Bacteria","unclassified Lactobacillales","unclassified Unassigned","Weeksellaceae","Low Abundance Families (< 0.1 %)"))+
  guides(fill=guide_legend(nrow=5), bycol=T)+
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
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.key.size= unit(4,"mm"),
        axis.text.x = element_blank())

stackeddendrophylum_sl_plot <- ggarrange(sl_dendro_plot, sl_phyla_plot_individual, ncol=1)
stackeddendrofamily_sl_plot <- ggarrange(sl_dendro_plot, sl_Family_filt_plot, ncol=1)

ggsave("stackeddendrophylum_sl_plot.tiff", path = "~/Desktop", plot =stackeddendrophylum_sl_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)


################################### BRD PATHOGENS ##############################
################################################################################

mann_ra_sl_melt <- subset_taxa(subsetl_genus, Genus=="Mannheimia") %>%
  psmelt()

kruskal_test(mann_ra_sl_melt, Abundance~Number.in.Pool)

mann_ra_sl_melt %>%
  group_by(Number.in.Pool)%>%
  summarise(mean=mean(Abundance),sd=sd(Abundance))

# PASTEUR
pasteurellaceae_ra_sl <- subset_taxa(subsetl_genus, Family=="Pasteurellaceae")
pasteurellaceae_ra_genus_sl <- tax_glom(pasteurellaceae_ra_sl, taxrank = "Genus", NArm = F) %>%
  psmelt()
pasteurellaceae_ra_family_sl <- tax_glom(pasteurellaceae_ra_sl, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(pasteurellaceae_ra_genus_sl$Genus))

pasteurella_sl_palette <- distinctColorPalette(k=10)

pasteurplot_sl <- ggplot(pasteurellaceae_ra_genus_sl, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteurellaceae_ra_family_sl, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_sl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
pasteurplot_sl

#Comparing Pools to Individuals
#Pools of 12

pasteur_L12_1 <- pasteurellaceae_ra_family_sl[which(pasteurellaceae_ra_family_sl$L_12.01=="Y"),]
kruskal_test(pasteur_L12_1,Abundance~Number.in.Pool) #NS

mann_sl <- subset_taxa(subsetl_genus,Genus=="Mannheimia")
mann_genus_sl <- tax_glom(mann_sl, taxrank = "Genus", NArm = F) %>%
  psmelt()
mann_L12_1 <- mann_genus_sl[which(mann_genus_sl$L_12.01=="Y"),]
kruskal_test(mann_L12_1,Abundance~Number.in.Pool)


pasteur_L12_2 <- pasteurellaceae_ra_family_sl[which(pasteurellaceae_ra_family_sl$L_12.02=="Y"),]
kruskal_test(pasteur_L12_2,Abundance~Number.in.Pool) #NS

pasteur_L12_3 <- pasteurellaceae_ra_family_sl[which(pasteurellaceae_ra_family_sl$L_12.03=="Y"),]
kruskal_test(pasteur_L12_3,Abundance~Number.in.Pool) #NS


pasteur_L12_4 <- pasteurellaceae_ra_family_sl[which(pasteurellaceae_ra_family_sl$L_12.04=="Y"),]
kruskal_test(pasteur_L12_4,Abundance~Number.in.Pool) #NS

pasteur_L12_5 <- pasteurellaceae_ra_family_sl[which(pasteurellaceae_ra_family_sl$L_12.05=="Y"),]
kruskal_test(pasteur_L12_5,Abundance~Number.in.Pool) #NS

pasteur_L12_6 <- pasteurellaceae_ra_family_sl[which(pasteurellaceae_ra_family_sl$L_12.06=="Y"),]
kruskal_test(pasteur_L12_6,Abundance~Number.in.Pool) #NS


# MYCOPLASMA
mycoplasmataceae_ra_sl <- subset_taxa(subsetl_genus, Family=="Mycoplasmataceae")
mycoplasmataceae_ra_genus_sl <- tax_glom(mycoplasmataceae_ra_sl, taxrank = "Genus", NArm = F) %>%
  psmelt()
mycoplasmataceae_ra_family_sl <- tax_glom(mycoplasmataceae_ra_sl, taxrank = "Family", NArm = F) %>%
  psmelt()

unique(mycoplasmataceae_ra_genus_sl$Genus)

mycoplasma_sl_palette <- distinctColorPalette(k=2)


mycopplot_sl <- ggplot(mycoplasmataceae_ra_genus_sl, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "mycoplasmataceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = mycoplasmataceae_ra_family_sl, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = mycoplasma_sl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mycopplot_sl 

#Comparing Pools to Individuals
#Pools of 12

mycop_L12_1 <- mycoplasmataceae_ra_family_sl[which(mycoplasmataceae_ra_family_sl$L_12.01=="Y"),]
kruskal_test(mycop_L12_1,Abundance~Number.in.Pool) #NS


mycop_L12_2 <- mycoplasmataceae_ra_family_sl[which(mycoplasmataceae_ra_family_sl$L_12.02=="Y"),]
kruskal_test(mycop_L12_2,Abundance~Number.in.Pool) #NS

mycop_L12_3 <- mycoplasmataceae_ra_family_sl[which(mycoplasmataceae_ra_family_sl$L_12.03=="Y"),]
kruskal_test(mycop_L12_3,Abundance~Number.in.Pool) #NS


mycop_L12_4 <- mycoplasmataceae_ra_family_sl[which(mycoplasmataceae_ra_family_sl$L_12.04=="Y"),]
kruskal_test(mycop_L12_4,Abundance~Number.in.Pool) #NS

mycop_L12_5 <- mycoplasmataceae_ra_family_sl[which(mycoplasmataceae_ra_family_sl$L_12.05=="Y"),]
kruskal_test(mycop_L12_5,Abundance~Number.in.Pool) #NS

mycop_L12_6 <- mycoplasmataceae_ra_family_sl[which(mycoplasmataceae_ra_family_sl$L_12.06=="Y"),]
kruskal_test(mycop_L12_6,Abundance~Number.in.Pool) #NS


# MORAXELLA
moraxellaceae_ra_sl <- subset_taxa(subsetl_genus, Family=="Moraxellaceae")
moraxellaceae_ra_genus_sl <- tax_glom(moraxellaceae_ra_sl, taxrank = "Genus", NArm = F) %>%
  psmelt()
moraxellaceae_ra_family_sl <- tax_glom(moraxellaceae_ra_sl, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(moraxellaceae_ra_genus_sl$Genus))

moraxella_sl_palette <- distinctColorPalette(k=11)


moraxplot_sl <- ggplot(moraxellaceae_ra_genus_sl, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = moraxellaceae_ra_family_sl, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values =moraxella_sl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
moraxplot_sl

#Comparing Pools to Individuals
#Pools of 12

moraxella_L12_1 <- moraxellaceae_ra_family_sl[which(moraxellaceae_ra_family_sl$L_12.01=="Y"),]
kruskal_test(moraxella_L12_1,Abundance~Number.in.Pool) #NS

ggplot(moraxella_L12_1, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_bar(aes(fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = moraxella_L12_1, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = moraxella_sl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = sh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))


moraxella_L12_2 <- moraxellaceae_ra_family_sl[which(moraxellaceae_ra_family_sl$L_12.02=="Y"),]
kruskal_test(moraxella_L12_2,Abundance~Number.in.Pool) #NS

moraxella_L12_3 <- moraxellaceae_ra_family_sl[which(moraxellaceae_ra_family_sl$L_12.03=="Y"),]
kruskal_test(moraxella_L12_3,Abundance~Number.in.Pool) #NS


moraxella_L12_4 <- moraxellaceae_ra_family_sl[which(moraxellaceae_ra_family_sl$L_12.04=="Y"),]
kruskal_test(moraxella_L12_4,Abundance~Number.in.Pool) #NS

moraxella_L12_5 <- moraxellaceae_ra_family_sl[which(moraxellaceae_ra_family_sl$L_12.05=="Y"),]
kruskal_test(moraxella_L12_5,Abundance~Number.in.Pool) #NS

moraxella_L12_6 <- moraxellaceae_ra_family_sl[which(moraxellaceae_ra_family_sl$L_12.06=="Y"),]
kruskal_test(moraxella_L12_6,Abundance~Number.in.Pool) #NS


##random
#high
###Phylum
randomh_phylum <- tax_glom(randomh.css.ra, taxrank = "Phylum", NArm = F) #37 phyla
rh_phylum_melt <- psmelt(randomh_phylum)
length(unique(rh_phylum_melt$Phylum))
###Class
randomh_class <- tax_glom(randomh.css.ra, taxrank = "Class", NArm = F) #93 classes
rh_class_melt <- psmelt(randomh_class)
length(unique(rh_class_melt$Class))
###Order
randomh_order <- tax_glom(randomh.css.ra, taxrank = "Order", NArm=F) #223 orders


###Family
randomh_family <- tax_glom(randomh.css.ra, taxrank = "Family",NArm = F) #444 families
randomh_Family_10 <- merge_less_than_top(randomh_family, top=10)
randomh_Family_10_melt <- psmelt(randomh_Family_10)

unique(rh_Family_top10_melt$Family)

###Genus
randomh_genus <- tax_glom(randomh.css.ra, taxrank = "Genus", NArm = F) #1226 genera


#Make palettes
rhphylumpalette <- distinctColorPalette(k=37)
rhclasspalette <- distinctColorPalette(k=93)
rhorderpalette <- distinctColorPalette(k=223)
rhfamiltypalette <- distinctColorPalette(k=444)
rhgenuspalette <- distinctColorPalette(k=1226)

randomh_phylum_filt <- merge_low_abundance(randomh_phylum, threshold = 1) #4 phyla
rh_phylum_filt_melt <- psmelt(randomh_phylum_filt)
unique(rh_phylum_filt_melt$Phylum) ##"Firmicutes"       "Bacteroidota"     "Proteobacteria"   "Actinobacteriota" "zzzOther " "unclassified Bacteria"

abundantphyla_palette <- distinctColorPalette(k=5)

rh_phyla_abundance_order <- c("Actinobacteriota", "Bacteroidota","Firmicutes","Proteobacteria")

rh_major_phyla_plot <- ggplot(rh_phylum_filt_melt, aes(x= Phylum, y= Abundance, fill= Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Phyla (>1 %)", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = rh_phylum_filt_melt, aes(x= Phylum, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge") +
  scale_x_discrete(limits = rh_phyla_abundance_order ) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "right",
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

##Stats
#Actinobacteriota
rh_phylum_filt_Actinobacteriota <- subset(rh_phylum_filt_melt, Phylum=="Actinobacteriota")
kruskal_test(rh_phylum_filt_Actinobacteriota,Abundance~Number.in.Pool) #NS, p= 0.275

#Bacteroidota
rh_phylum_filt_Bacteroidota <- subset(rh_phylum_filt_melt, Phylum=="Bacteroidota")
kruskal_test(rh_phylum_filt_Bacteroidota,Abundance~Number.in.Pool) #NS, P=0.371

#Firmicutes
rh_phylum_filt_Firmicutes <- subset(rh_phylum_filt_melt, Phylum=="Firmicutes")
kruskal_test(rh_phylum_filt_Firmicutes,Abundance~Number.in.Pool) #NS, P=0.975

#Proteobacteria
rh_phylum_filt_Proteobacteria <- subset(rh_phylum_filt_melt, Phylum=="Proteobacteria")
kruskal_test(rh_phylum_filt_Proteobacteria,Abundance~Number.in.Pool) #NS, P=0.708


#[1] "Mycoplasmataceae"      "Chitinophagaceae"      "Pasteurellaceae"       "Microbacteriaceae"     "zzzOther"              "Moraxellaceae"         "Corynebacteriaceae"   
#[8] "Prevotellaceae"        "unclassified Bacteria" "Oscillospiraceae"      "Lachnospiraceae"

rh_family_abundace_order <- c("Mycoplasmataceae","Pasteurellaceae","Microbacteriaceae","Chitinophagaceae","Moraxellaceae")

randomh_Family_10_sd <- randomh_Family_10_melt %>%
  group_by(Family,Number.in.Pool)%>%
  summarise(meanRA=mean(Abundance), sd=sd(Abundance))

##Family
rh_major_family_plot <-ggplot(randomh_Family_10_melt, aes(x= Family, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Top 10 Families", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = randomh_Family_10_sd, aes(x= Family, y= meanRA, ymin=meanRA, ymax=meanRA+sd), width = 0.55, position = position_dodge(0.9), color="black", linewidth=0.8)+
  geom_bar(stat = "summary", position = "dodge", color="black", linewidth=0.8) +
  scale_x_discrete(limits = rh_family_abundace_order) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
  scale_color_manual(values=numberinpoolpalette)+
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

#Stats
#Mycoplasmataceae
rh_Family_filt_Mycoplasmataceae <- subset(randomh_Family_10_melt, Family=="Mycoplasmataceae")
kruskal_test(rh_Family_filt_Mycoplasmataceae,Abundance~Number.in.Pool) #NS, p= 0.968
dunn_test(rh_Family_filt_Mycoplasmataceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Mycoplasmataceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (72.3, 26.8), 3(82.5  13.6), 6(79.3  16.1), 12(78.4, 13.3)
rh_Family_filt_Mycoplasmataceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #73.60891 25.03991

#Chitinophagaceae
rh_Family_filt_Chitinophagaceae <- subset(randomh_Family_10_melt, Family=="Chitinophagaceae")
kruskal_test(rh_Family_filt_Chitinophagaceae,Abundance~Number.in.Pool) #Sig, p= 0.204
dunn_test(rh_Family_filt_Chitinophagaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences, 1 and 3 0.069
rh_Family_filt_Chitinophagaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.79, 11.8), 3(0.229  0.216), 6 (0.762  1.02), 12 (1.26, 1.32)
rh_Family_filt_Chitinophagaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #2.436913 10.78406

#Pasteurellaceae
rh_Family_filt_Pasteurellaceae <- subset(randomh_Family_10_melt, Family=="Pasteurellaceae")
kruskal_test(rh_Family_filt_Pasteurellaceae,Abundance~Number.in.Pool) #NS, 0.596
dunn_test(rh_Family_filt_Pasteurellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Pasteurellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (9.93, 14.2), 3(9.04  9.26), 6 (10.9  15.5 ), 12 (10.6, 11.9)
rh_Family_filt_Pasteurellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #9.97097 13.70427

#Microbacteriaceae
rh_Family_filt_Microbacteriaceae <- subset(randomh_Family_10_melt, Family=="Microbacteriaceae")
kruskal_test(rh_Family_filt_Microbacteriaceae,Abundance~Number.in.Pool) #NS, 0.11
dunn_test(rh_Family_filt_Microbacteriaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Microbacteriaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.54, 8.49), 3(0.242 0.414), 6 (2.07  3.49 ), 12 (1.59  2.20)
rh_Family_filt_Microbacteriaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(2.715898 8.001184)

#Moraxellaceae
rh_Family_filt_Moraxellaceae <- subset(randomh_Family_10_melt, Family=="Moraxellaceae")
kruskal_test(rh_Family_filt_Moraxellaceae,Abundance~Number.in.Pool) #NS, 0.117
dunn_test(rh_Family_filt_Moraxellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Moraxellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.27  7.50 ), 3(5.09  6.90 ), 6 (0.908 0.907), 12 (1.73  2.77)
rh_Family_filt_Moraxellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(2.330398 7.040223)

#Corynebacteriaceae
rh_Family_filt_Corynebacteriaceae <- subset(randomh_Family_10_melt, Family=="Corynebacteriaceae")
kruskal_test(rh_Family_filt_Corynebacteriaceae,Abundance~Number.in.Pool) #NS,0.68
dunn_test(rh_Family_filt_Corynebacteriaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Corynebacteriaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.593 1.74), 3(0.113 0.138), 6 (0.201 0.289), 12 (0.160 0.0752)
rh_Family_filt_Corynebacteriaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.516472 1.593222)

#Prevotellaceae
rh_Family_filt_Prevotellaceae <- subset(randomh_Family_10_melt, Family=="Prevotellaceae")
kruskal_test(rh_Family_filt_Prevotellaceae,Abundance~Number.in.Pool) #NS,0.644
dunn_test(rh_Family_filt_Prevotellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Prevotellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.718 1.43 ), 3(0.305 0.369), 6 (0.637 0.782), 12 (0.543 0.540)
rh_Family_filt_Prevotellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.6781369 1.31999)

#unclassified Bacteria
rh_Family_filt_unclassifiedBacteria <- subset(randomh_Family_10_melt, Family=="unclassified Bacteria")
kruskal_test(rh_Family_filt_unclassifiedBacteria,Abundance~Number.in.Pool) #NS,0.658
dunn_test(rh_Family_filt_unclassifiedBacteria, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_unclassifiedBacteria %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.385 1.13), 3(0.0845 0.0318), 6 (0.0834 0.0511), 12 (0.587 1.11)
rh_Family_filt_unclassifiedBacteria %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.3674297 1.068518)

#Oscillospiraceae
rh_Family_filt_Oscillospiraceae <- subset(randomh_Family_10_melt, Family=="Oscillospiraceae")
kruskal_test(rh_Family_filt_Oscillospiraceae,Abundance~Number.in.Pool) #NS,0.75
dunn_test(rh_Family_filt_Oscillospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Oscillospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.703 1.47), 3(0.137 0.137), 6 (0.267 0.326), 12 (0.259 0.147)
rh_Family_filt_Oscillospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.618813 1.347414)

#Lachnospiraceae
rh_Family_filt_Lachnospiraceae <- subset(randomh_Family_10_melt, Family=="Lachnospiraceae")
kruskal_test(rh_Family_filt_Lachnospiraceae,Abundance~Number.in.Pool) #NS,0.728
dunn_test(rh_Family_filt_Lachnospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_Lachnospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.624 1.11), 3(0.140 0.116), 6 (0.341 0.310), 12 (0.286 0.108)
rh_Family_filt_Lachnospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.5591685 1.019127)

#zzzOther
rh_Family_filt_zzzOther <- subset(randomh_Family_10_melt, Family=="zzzOther")
kruskal_test(rh_Family_filt_zzzOther,Abundance~Number.in.Pool) #NS,0.623
dunn_test(rh_Family_filt_zzzOther, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rh_Family_filt_zzzOther %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #
rh_Family_filt_zzzOther %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance))


rh_w.hclust <- hclust(randomh_wunifrac.dist, method = "ward.D2")
rh_w.dendro <- as.dendrogram(rh_w.hclust)
rh_w.dendro.data <- dendro_data(rh_w.dendro, type = "rectangle")
rh_w.metadata_for_dendro <- as_tibble(randomh.css@sam_data)
rh_w.dendro.data$labels <- rh_w.dendro.data$labels %>%
  left_join(rh_w.metadata_for_dendro, by = c("label" = "MEG_ID"))



rh_w_dendro_sample_order <- rh_w.dendro.data$labels$label

rh_dendro_plot <-ggplot(rh_w.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = rh_w.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 5, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  geom_text(data = rh_w.dendro.data$labels, 
            aes(x=x,y=y, label=HRPool),size=3,position = position_nudge(y=-0.05,x=0),color="black",angle=90)+
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
rh_dendro_plot

# family RA plot for under dendro
rh_family_filt <- merge_low_abundance(randomh_family,threshold=1)
rh_family_filt_melt <- psmelt(rh_family_filt)

unique(rh_family_filt_melt$Family)

length(unique(sh_family_filt_melt$Family)) #6 families > 1 %
write.csv(sh_phylum_filt_melt$Phylum, "sh_phyla.csv")
write.csv(sh_phylum_filt_melt$mean_phylum_ra,"sh_phyla_meanRA.csv")


rh_phyla_plot_individual <- ggplot(rh_phylum_filt_melt, aes(x= MEG_ID, y= Abundance, fill= Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =high_phylum_filt_palette,labels=c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria","Low Abundance Phyla (< 1%)")) +
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
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

rh_phyla_plot_individual

stackeddendrophylum_rh_plot <- ggarrange(rh_dendro_plot, rh_phyla_plot_individual, ncol=1)

ggsave("stackeddendrophylum_sl_plot.tiff", path = "~/Desktop", plot =stackeddendrophylum_sl_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)


################################### BRD PATHOGENS ##############################
################################################################################

mann_ra_rh_melt <- subset_taxa(randomh_genus, Genus=="Mannheimia") %>%
  psmelt()

kruskal_test(mann_ra_rh_melt, Abundance~Number.in.Pool)

mann_ra_rh_melt %>%
  group_by(Number.in.Pool)%>%
  summarise(mean=mean(Abundance),sd=sd(Abundance))

# PASTEUR
pasteurellaceae_ra_rh <- subset_taxa(randomh_genus, Family=="Pasteurellaceae")
pasteurellaceae_ra_genus_rh <- tax_glom(pasteurellaceae_ra_rh, taxrank = "Genus", NArm = F) %>%
  psmelt()
pasteurellaceae_ra_family_rh <- tax_glom(pasteurellaceae_ra_rh, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(pasteurellaceae_ra_genus_rh$Genus))

pasteurella_rh_palette <- distinctColorPalette(k=8)

pasteurplot_rh <- ggplot(pasteurellaceae_ra_genus_rh, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteurellaceae_ra_family_rh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_rh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
pasteurplot_rh

pasteurplot_rh_grouped <- ggplot(pasteurellaceae_ra_genus_rh, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_errorbar(data =pasteurellaceae_ra_family_rh, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  scale_fill_manual(values = pasteurella_rh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
pasteurplot_rh_grouped

kruskal_test(pasteurellaceae_ra_family_rh, Abundance~Number.in.Pool)

#Comparing Pools to Individuals
#Pools of 12

pasteur_HR6_1 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_06.01=="Y"),]
wilcox_test(pasteur_HR6_1,Abundance~Type) #NS

ggplot(pasteur_HR6_1, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteur_HR6_1, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_rh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

pasteur_HR6_2 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_06.02=="Y"),]
wilcox_test(pasteur_HR6_2,Abundance~Type) #NS

ggplot(pasteur_HR6_2, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteur_HR6_2, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_rh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

pasteur_HR6_3 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_06.03=="Y"),]
wilcox_test(pasteur_HR6_3,Abundance~Type) #NS

pasteur_HR6_4 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_06.04=="Y"),]
wilcox_test(pasteur_HR6_4,Abundance~Type) #NS

pasteur_HR6_5 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_06.05=="Y"),]
wilcox_test(pasteur_HR6_5,Abundance~Type) #Pool removed with trimming

pasteur_HR6_6 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_06.06=="Y"),]
wilcox_test(pasteur_HR6_6,Abundance~Type) #NSPool removed with trimming

##Pools of 3
pasteur_HR03_1 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_03.01=="Y"),]
wilcox_test(pasteur_HR03_1,Abundance~Type) #NS

ggplot(pasteur_HR03_1, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteur_HR03_1, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_rh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

pasteur_HR03_2 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_03.02=="Y"),]
wilcox_test(pasteur_HR03_2,Abundance~Type) #NS

ggplot(pasteur_HR03_2, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteur_HR03_2, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_rh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

pasteur_HR03_3 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_03.03=="Y"),]
wilcox_test(pasteur_HR03_3,Abundance~Type) #NS

pasteur_HR03_4 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_03.04=="Y"),]
wilcox_test(pasteur_HR03_4,Abundance~Type) #NS

pasteur_HR03_5 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_03.05=="Y"),]
wilcox_test(pasteur_HR03_5,Abundance~Type) #Pool removed with trimming

pasteur_HR03_6 <- pasteurellaceae_ra_family_rh[which(pasteurellaceae_ra_family_rh$HR_03.06=="Y"),]
wilcox_test(pasteur_HR03_6,Abundance~Type) #Pool removed with trimming



# MYCOPLASMA
mycoplasmataceae_ra_rh <- subset_taxa(randomh_genus, Family=="Mycoplasmataceae")
mycoplasmataceae_ra_genus_rh <- tax_glom(mycoplasmataceae_ra_rh, taxrank = "Genus", NArm = F) %>%
  psmelt()
mycoplasmataceae_ra_family_rh <- tax_glom(mycoplasmataceae_ra_rh, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(mycoplasmataceae_ra_genus_rh$Genus))

mycoplasma_rh_palette <- distinctColorPalette(k=2)


mycopplot_rh <- ggplot(mycoplasmataceae_ra_genus_rh, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "Mycoplasmataceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = mycoplasmataceae_ra_family_rh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = mycoplasma_sh_palette) +
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
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mycopplot_rh

mycopplot_rh_grouped <- ggplot(mycoplasmataceae_ra_genus_rh, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "mycoplasmataceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = mycoplasmataceae_ra_family_rh, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = mycoplasma_sh_palette) +
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
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mycopplot_rh_grouped

kruskal_test(mycoplasmataceae_ra_genus_rh,Abundance~Number.in.Pool)

#Comparing Pools to Individuals
#Pools of 6
mycop_H06_1 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_06.01=="Y"),]
wilcox_test(mycop_H06_1,Abundance~Type) #NS

mycop_H06_2 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_06.02=="Y"),]
wilcox_test(mycop_H06_2,Abundance~Type) #NS

mycop_H06_3 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_06.03=="Y"),]
wilcox_test(mycop_H06_3,Abundance~Type) #NS

mycop_H06_4 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_06.04=="Y"),]
wilcox_test(mycop_H06_4,Abundance~Type) #NS

mycop_H06_5 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_06.05=="Y"),]
wilcox_test(mycop_H06_5,Abundance~Type) #NS...Pool trimmed out

mycop_H06_6 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_06.06=="Y"),]
wilcox_test(mycop_H06_6,Abundance~Type) #NS...Pool trimmed out

#Pool of 3
mycop_H03_1 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_03.01=="Y"),]
wilcox_test(mycop_H03_1,Abundance~Type) #NS

mycop_H03_2 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_03.02=="Y"),]
wilcox_test(mycop_H03_2,Abundance~Type) #NS

mycop_H03_3 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_03.03=="Y"),]
wilcox_test(mycop_H03_3,Abundance~Type) #NS

mycop_H03_4 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_03.04=="Y"),]
wilcox_test(mycop_H03_4,Abundance~Type) #NS

mycop_H03_5 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_03.05=="Y"),]
wilcox_test(mycop_H03_5,Abundance~Type) #NS

mycop_H03_6 <- mycoplasmataceae_ra_family_rh[which(mycoplasmataceae_ra_family_rh$HR_03.06=="Y"),]
wilcox_test(mycop_H03_6,Abundance~Type) #NS...Pool trimmed out


# MORAXELLA
moraxellaceae_ra_rh <- subset_taxa(randomh_genus, Family=="Moraxellaceae")
moraxellaceae_ra_genus_rh <- tax_glom(moraxellaceae_ra_rh, taxrank = "Genus", NArm = F) %>%
  psmelt()
moraxellaceae_ra_family_rh <- tax_glom(moraxellaceae_ra_rh, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(moraxellaceae_ra_genus_rh$Genus))

moraxella_rh_palette <- distinctColorPalette(k=8)


moraxplot_rh <- ggplot(moraxellaceae_ra_genus_rh, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = moraxellaceae_ra_family_rh, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values =moraxella_sh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

moraxplot_rh

moraxplot_rh_grouped <- ggplot(moraxellaceae_ra_genus_rh, aes(x= Number.in.Pool, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_errorbar(data = moraxellaceae_ra_family_rh, aes(x= Number.in.Pool, y= Abundance), stat = "summary", width = 0.55) +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  scale_fill_manual(values =moraxella_sh_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))

moraxplot_rh_grouped

kruskal_test(moraxellaceae_ra_family_rh,Abundance~Number.in.Pool) #NS, P=0.117

#Comparing Pools to Individuals

#Pools of 6

moraxella_HR6_1 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_06.01=="Y"),]
wilcox_test(moraxella_HR6_1,Abundance~Type) #NS

moraxella_HR6_2 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_06.02=="Y"),]
wilcox_test(moraxella_HR6_2,Abundance~Type) #NS

moraxella_HR6_3 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_06.03=="Y"),]
wilcox_test(moraxella_HP6_3,Abundance~Type) #NS


moraxella_HR6_4 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_06.04=="Y"),]
wilcox_test(moraxella_HR6_4,Abundance~Type) #NS

moraxella_HR6_5 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_06.05=="Y"),]
wilcox_test(moraxella_HR6_5,Abundance~Type) #NS, pool trimmed out

moraxella_HR6_6 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_06.06=="Y"),]
wilcox_test(moraxella_HR6_6,Abundance~Type) #NS, pool trimmed out

#Pools of 3
moraxella_HR3_1 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_03.01=="Y"),]
wilcox_test(moraxella_HR3_1,Abundance~Type) #NS


moraxella_HR3_2 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_03.02=="Y"),]
wilcox_test(moraxella_HR3_2,Abundance~Type) #NS

moraxella_HR3_3 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_03.03=="Y"),]
wilcox_test(moraxella_HR3_3,Abundance~Type) #NS


moraxella_HR3_4 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_03.04=="Y"),]
wilcox_test(moraxella_HR3_4,Abundance~Type) #NS

moraxella_HR3_5 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_03.05=="Y"),]
wilcox_test(moraxella_HR3_5,Abundance~Type) #NS

moraxella_HR3_6 <- moraxellaceae_ra_family_rh[which(moraxellaceae_ra_family_rh$HR_03.06=="Y"),]
wilcox_test(moraxella_HR3_6,Abundance~Type) #NS, pool trimmed out


ggarrange(dendro1, plot1, pasteurplot_rh, mycopplot_rh, moraxplot_rh, ncol = 1)


############## random LOW ###############
###make some palettes
##random
#low
randoml_phylum <- tax_glom(randoml.css.ra, taxrank = "Phylum", NArm = F)#58 phyla
randoml_class <- tax_glom(randoml.css.ra, taxrank = "Class", NArm = F) #158 classes
randoml_order <- tax_glom(randoml.css.ra, taxrank = "Order", NArm=F) #374 orders
randoml_family <- tax_glom(randoml.css.ra, taxrank = "Family",NArm = F) #715 families

randoml_Family_10 <- merge_less_than_top(randoml_family, top=10)
randoml_Family_10_melt <- psmelt(randoml_Family_10)


randoml_genus <- tax_glom(randoml.css.ra, taxrank = "Genus", NArm = F) #1817 genera


rlphylumpalette <- distinctColorPalette(k=60)
#slclasspalette <- distinctColorPalette(k=103)
#slorderpalette <- distinctColorPalette(k=202)
#slfamiltypalette <- distinctColorPalette(k=380)
#slgenuspalette <- distinctColorPalette(k=828)

rl_phylum_melt <- psmelt(randoml_phylum)
rl_genus_melt <- psmelt(randoml_genus)

randoml_phylum_filt <- merge_low_abundance(randoml_phylum, threshold = 1) #4 phyla
rl_phylum_filt_melt <- psmelt(randoml_phylum_filt)
unique(rl_phylum_filt_melt$Phylum) ##"Firmicutes"       "Bacteroidota"     "Proteobacteria"   "Actinobacteriota"  "unclassified Bacteria"

abundantphyla_palette <- distinctColorPalette(k=5)

rl_phyla_abundance_order <- c("Actinobacteriota", "Bacteroidota","Firmicutes","Proteobacteria", "unclassified Bacteria")

rl_major_phyla_plot <- ggplot(rl_phylum_filt_melt, aes(x= Phylum, y= Abundance, fill= Number.in.Pool, color=Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Phyla (>1 %)", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = rl_phylum_filt_melt, aes(x= Phylum, y= Abundance, color=Number.in.Pool), stat = "summary", width = 0.55, position = position_dodge(0.9))+
  geom_bar(stat = "summary", position = "dodge") +
  scale_x_discrete(limits = rl_phyla_abundance_order ) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))+
  scale_color_manual(values=numberinpoolpalette)+
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme_minimal()+
  theme(legend.position = "right",
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

##Stats
#Actinobacteriota
rl_phylum_filt_Actinobacteriota <- subset(rl_phylum_filt_melt, Phylum=="Actinobacteriota")
kruskal_test(rl_phylum_filt_Actinobacteriota,Abundance~Number.in.Pool) #NS, p= 0.254

#Bacteroidota
rl_phylum_filt_Bacteroidota <- subset(rl_phylum_filt_melt, Phylum=="Bacteroidota")
kruskal_test(rl_phylum_filt_Bacteroidota,Abundance~Number.in.Pool) #NS, P=0.209

#Firmicutes
rl_phylum_filt_Firmicutes <- subset(rl_phylum_filt_melt, Phylum=="Firmicutes")
kruskal_test(rl_phylum_filt_Firmicutes,Abundance~Number.in.Pool) #NS, P=0.676

#Proteobacteria
rl_phylum_filt_Proteobacteria <- subset(rl_phylum_filt_melt, Phylum=="Proteobacteria")
kruskal_test(rl_phylum_filt_Proteobacteria,Abundance~Number.in.Pool) #NS, P=0.755

#unclassified Bacteria
rl_phylum_filt_unclassified<- subset(rl_phylum_filt_melt, Phylum=="unclassified Bacteria")
kruskal_test(rl_phylum_filt_unclassified,Abundance~Number.in.Pool) #NS, P=0.021

rl_w.hclust <- hclust(randoml_wunifrac.dist, method = "ward.D2")
rl_w.dendro <- as.dendrogram(rl_w.hclust)
rl_w.dendro.data <- dendro_data(rl_w.dendro, type = "rectangle")
rl_w.metadata_for_dendro <- as_tibble(randoml.css@sam_data)
rl_w.dendro.data$labels <- rl_w.dendro.data$labels %>%
  left_join(rl_w.metadata_for_dendro, by = c("label" = "MEG_ID"))

rl_family_abundace_order <- c("Mycoplasmataceae","Pasteurellaceae","Moraxellaceae", "Chitinophagaceae","Microbacteriaceae")

randoml_Family_10_sd <- randoml_Family_10_melt %>%
  group_by(Family,Number.in.Pool)%>%
  summarise(meanRA=mean(Abundance), sd=sd(Abundance))

##Family
rl_major_family_plot <-ggplot(randoml_Family_10_melt, aes(x= Family, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
  labs(y= "Relative Abundance (%)", title = "Top 10 Families", color="Pool Size", fill="Pool Size", alpha="Pool Size") +
  geom_errorbar(data = randoml_Family_10_sd, aes(x= Family, y= meanRA, ymin=meanRA, ymax=meanRA+sd), width = 0.55, position = position_dodge(0.9), linewidth=0.8, color="black")+
  geom_bar(stat = "summary", position = "dodge", linewidth=0.8, color="black") +
  scale_x_discrete(limits = rl_family_abundace_order) +
  scale_fill_manual(values =numberinpoolpalette) +
  scale_alpha_manual(values=c(0.7,0.7,0.7,0.7))+
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

#Stats
#Mycoplasmataceae
rl_Family_filt_Mycoplasmataceae <- subset(randoml_Family_10_melt, Family=="Mycoplasmataceae")
kruskal_test(rl_Family_filt_Mycoplasmataceae,Abundance~Number.in.Pool) #NS, p= 0.601
dunn_test(rl_Family_filt_Mycoplasmataceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Mycoplasmataceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #1 (69.4 28.7), 3(80.2 17.9 ), 6(66.4 15.0), 12(68.1  9.43)
rl_Family_filt_Mycoplasmataceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #69.82553 26.45816

#Chitinophagaceae
rl_Family_filt_Chitinophagaceae <- subset(randoml_Family_10_melt, Family=="Chitinophagaceae")
kruskal_test(rl_Family_filt_Chitinophagaceae,Abundance~Number.in.Pool) #NS, p= 0.332
dunn_test(rl_Family_filt_Chitinophagaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Chitinophagaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (2.49  11.0), 3(0.261  0.237), 6 (2.89   5.51), 12 (1.23   1.79)
rl_Family_filt_Chitinophagaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #2.278461 9.910074

#Pasteurellaceae
rl_Family_filt_Pasteurellaceae <- subset(randoml_Family_10_melt, Family=="Pasteurellaceae")
kruskal_test(rl_Family_filt_Pasteurellaceae,Abundance~Number.in.Pool) #NS, 0.326
dunn_test(rl_Family_filt_Pasteurellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Pasteurellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (10.0  14.4), 3(9.15  9.12), 6 (13.3  10.2), 12 (9.30  6.36)
rl_Family_filt_Pasteurellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #10.13048 13.38088

#Microbacteriaceae
rl_Family_filt_Microbacteriaceae <- subset(randoml_Family_10_melt, Family=="Microbacteriaceae")
kruskal_test(rl_Family_filt_Microbacteriaceae,Abundance~Number.in.Pool) #NS, 0.164
dunn_test(rl_Family_filt_Microbacteriaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Microbacteriaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (1.75  5.47), 3(1.27  1.55 ), 6 (2.14  4.38), 12 ( 0.771 0.808)
rl_Family_filt_Microbacteriaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(1.673536 5.01092)

#Moraxellaceae
rl_Family_filt_Moraxellaceae <- subset(randoml_Family_10_melt, Family=="Moraxellaceae")
kruskal_test(rl_Family_filt_Moraxellaceae,Abundance~Number.in.Pool) #NS, 0.175
dunn_test(rl_Family_filt_Moraxellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences, 1 and 6 and 1 and 12 0.157
rl_Family_filt_Moraxellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (4.38 11.4), 3(3.71   7.34), 6 (0.760  0.721), 12 (2.55   2.18)
rl_Family_filt_Moraxellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(3.961652 10.34291)

#unclassified Unassigned
rl_Family_filt_unclassifiedUnassigned <- subset(randoml_Family_10_melt, Family=="unclassified Unassigned")
kruskal_test(rl_Family_filt_unclassifiedUnassigned,Abundance~Number.in.Pool) #Sig,0.000313
dunn_test(rl_Family_filt_unclassifiedUnassigned, Abundance~Number.in.Pool, p.adjust.method = "BH") #Pairwise differences, 1 different fromall (0.0192)
rl_Family_filt_unclassifiedUnassigned %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.555 1.29), 3(0.773 0.351), 6 (1.12  0.749), 12 (4.83  5.60)
rl_Family_filt_unclassifiedUnassigned %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.9000739 2.081703)

#Prevotellaceae
rl_Family_filt_Prevotellaceae <- subset(randoml_Family_10_melt, Family=="Prevotellaceae")
kruskal_test(rl_Family_filt_Prevotellaceae,Abundance~Number.in.Pool) #NS,0.103
dunn_test(rl_Family_filt_Prevotellaceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Prevotellaceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.532 0.935), 3(0.143 0.123), 6 (0.767 0.681), 12 (0.511 0.402)
rl_Family_filt_Prevotellaceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.5204521 0.8629437)

#unclassified Bacteria
rl_Family_filt_unclassifiedBacteria <- subset(randoml_Family_10_melt, Family=="unclassified Bacteria")
kruskal_test(rl_Family_filt_unclassifiedBacteria,Abundance~Number.in.Pool) #NS,0.081
dunn_test(rl_Family_filt_unclassifiedBacteria, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_unclassifiedBacteria %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.957 2.31), 3(0.430 0.183), 6 (0.799 0.613), 12 (4.58  5.56)
rl_Family_filt_unclassifiedBacteria %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(1.157225 2.631347)

#Oscillospiraceae
rl_Family_filt_Oscillospiraceae <- subset(randoml_Family_10_melt, Family=="Oscillospiraceae")
kruskal_test(rl_Family_filt_Oscillospiraceae,Abundance~Number.in.Pool) #NS,0.297
dunn_test(rl_Family_filt_Oscillospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Oscillospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (0.857 1.82), 3(0.253 0.275), 6 (1.38  2.58), 12 (0.846 0.997)
rl_Family_filt_Oscillospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(0.8505621 1.763119)

#Lachnospiraceae
rl_Family_filt_Lachnospiraceae <- subset(randoml_Family_10_melt, Family=="Lachnospiraceae")
kruskal_test(rl_Family_filt_Lachnospiraceae,Abundance~Number.in.Pool) #NS,0.183
dunn_test(rl_Family_filt_Lachnospiraceae, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_Lachnospiraceae %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #RA very close; 1 (1.07  2.57), 3(0.245 0.288), 6 (1.99  3.84), 12 (1.22  1.67)
rl_Family_filt_Lachnospiraceae %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #(1.089206 2.5213)

#zzzOther
rl_Family_filt_zzzOther <- subset(randoml_Family_10_melt, Family=="zzzOther")
kruskal_test(rl_Family_filt_zzzOther,Abundance~Number.in.Pool) #NS,0.573
dunn_test(rl_Family_filt_zzzOther, Abundance~Number.in.Pool, p.adjust.method = "BH") #No pairwise differences
rl_Family_filt_zzzOther %>%
  group_by(Number.in.Pool)%>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance)) #
rl_Family_filt_zzzOther %>%
  summarise(n=n(),MeanAbundance=mean(Abundance),SD=sd(Abundance))

dendro4<-ggplot(rl_w.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = rl_w.dendro.data$labels, 
             aes(x=x,y=y, colour = Number.in.Pool, fill= Number.in.Pool),
             size = 11, shape=22, stroke =1.5, position = position_nudge(y=-0.18)) +
  scale_y_continuous(limits = c(-.25,2)) +
  #scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
dendro4

# family RA plot for under dendro
rl_phylum_filt_melt <- rl_phylum_melt %>%
  group_by(Phylum) %>%
  mutate(mean_phylum_ra = mean(Abundance))

rl_phylum_filt_melt$Phylum[rl_phylum_filt_melt$mean_phylum_ra < 0.01] <- "zLow Abundance Phyla (<0.01 %)"
length(unique(rl_phylum_filt_melt$Phylum)) #20 phyla >0.01
write.csv(rl_phylum_filt_melt$Phylum, "rl_phyla.csv")
write.csv(rl_phylum_filt_melt$mean_phylum_ra,"rl_phyla_meanRA.csv")
filtphylapalette_rl <- distinctColorPalette(k=20)

rl_w_dendro_sample_order <- rl_w.dendro.data$labels$label

plot4<- ggplot(rl_phylum_filt_melt, aes(x= MEG_ID, y= Abundance, fill= Phylum)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = rl_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = filtphylapalette_rl) +
  #facet_wrap(~Number.in.Pool,scales = "free_x")+
  theme(#legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.y = element_line(size = 0.7, colour = "black"),
    axis.ticks.y = element_line(colour = "black", size = 0.75),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=-45))

plot4


ggarrange(dendro4, plot4, ncol=1)


################################### BRD PATHOGENS ##############################
################################################################################
mann_ra_rl_melt <- subset_taxa(randoml_genus, Genus=="Mannheimia") %>%
  psmelt()

kruskal_test(mann_ra_rl_melt, Abundance~Number.in.Pool) #P=0.07
dunn_test(mann_ra_rl_melt, Abundance~Number.in.Pool, p.adjust.method = "BH")

mann_ra_rl_melt %>%
  group_by(Number.in.Pool)%>%
  summarise(mean=mean(Abundance),sd=sd(Abundance))



# PASTEUR
pasteurellaceae_ra_rl <- subset_taxa(randoml_genus, Family=="Pasteurellaceae")

pasteurellaceae_ra_genus_rl <- tax_glom(pasteurellaceae_ra_rl, taxrank = "Genus", NArm = F) %>%
  psmelt()
pasteurellaceae_ra_family_rl <- tax_glom(pasteurellaceae_ra_rl, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(pasteurellaceae_ra_genus_rl$Genus))

pasteurella_rl_palette <- distinctColorPalette(k=10)

pasteurplot_rl <- ggplot(pasteurellaceae_ra_genus_rl, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "PASTEURELLACEAE") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = pasteurellaceae_ra_family_rl, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = pasteurella_rl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = rl_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
pasteurplot_rl

#Comparing Pools to Individuals
#Pools of 6

pasteur_LR06_1 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_06.01=="Y"),]
wilcox_test(pasteur_LR06_1,Abundance~Type) #NS


pasteur_LR06_2 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_06.02=="Y"),]
wilcox_test(pasteur_LR06_2,Abundance~Type) #NS

pasteur_LR06_3 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_06.03=="Y"),]
wilcox_test(pasteur_LR06_3,Abundance~Type) #NS


pasteur_LR06_4 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_06.04=="Y"),]
wilcox_test(pasteur_LR06_4,Abundance~Type) #NS

pasteur_LR06_5 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_06.05=="Y"),]
wilcox_test(pasteur_LR06_5,Abundance~Type) #NS

pasteur_LR06_6 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_06.06=="Y"),]
wilcox_test(pasteur_LR06_6,Abundance~Type) #NS

#Pools of 3

pasteur_LR03_1 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_03.01=="Y"),]
wilcox_test(pasteur_LR03_1,Abundance~Type) #NS


pasteur_LR03_2 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_03.02=="Y"),]
wilcox_test(pasteur_LR03_2,Abundance~Type) #NS

pasteur_LR03_3 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_03.03=="Y"),]
wilcox_test(pasteur_LR03_3,Abundance~Type) #NS


pasteur_LR03_4 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_03.04=="Y"),]
wilcox_test(pasteur_LR03_4,Abundance~Type) #NS

pasteur_LR03_5 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_03.05=="Y"),]
wilcox_test(pasteur_LR03_5,Abundance~Type) #NS

pasteur_LR03_6 <- pasteurellaceae_ra_family_rl[which(pasteurellaceae_ra_family_rl$LR_03.06=="Y"),]
wilcox_test(pasteur_LR03_6,Abundance~Type) #NS


# MYCOPLASMA
mycoplasmataceae_ra_rl <- subset_taxa(randoml_genus, Family=="Mycoplasmataceae")
mycoplasmataceae_ra_genus_rl <- tax_glom(mycoplasmataceae_ra_rl, taxrank = "Genus", NArm = F) %>%
  psmelt()
mycoplasmataceae_ra_family_rl <- tax_glom(mycoplasmataceae_ra_rl, taxrank = "Family", NArm = F) %>%
  psmelt()

unique(mycoplasmataceae_ra_genus_rl$Genus)

mycoplasma_rl_palette <- distinctColorPalette(k=2)


mycopplot_rl <- ggplot(mycoplasmataceae_ra_genus_rl, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "mycoplasmataceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = mycoplasmataceae_ra_family_rl, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values = mycoplasma_rl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  scale_x_discrete(limits = rl_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
mycopplot_rl 

#Comparing Pools to Individuals
#Pools of 6

mycop_LR06_1 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_06.01=="Y"),]
wilcox_test(mycop_LR06_1,Abundance~Type) #NS


mycop_LR06_2 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_06.02=="Y"),]
wilcox_test(mycop_LR06_2,Abundance~Type) #NS

mycop_LR06_3 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_06.03=="Y"),]
wilcox_test(mycop_LR06_3,Abundance~Type) #NS


mycop_LR06_4 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_06.04=="Y"),]
wilcox_test(mycop_LR06_4,Abundance~Type) #NS

mycop_LR06_5 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_06.05=="Y"),]
wilcox_test(mycop_LR06_5,Abundance~Type) #NS

mycop_LR06_6 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_06.06=="Y"),]
wilcox_test(mycop_LR06_6,Abundance~Type) #NS

#Pools of 3

mycop_LR03_1 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_03.01=="Y"),]
wilcox_test(mycop_LR03_1,Abundance~Type) #NS


mycop_LR03_2 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_03.02=="Y"),]
wilcox_test(mycop_LR03_2,Abundance~Type) #NS

mycop_LR03_3 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_03.03=="Y"),]
wilcox_test(mycop_LR03_3,Abundance~Type) #NS


mycop_LR03_4 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_03.04=="Y"),]
wilcox_test(mycop_LR03_4,Abundance~Type) #NS

mycop_LR03_5 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_03.05=="Y"),]
wilcox_test(mycop_LR03_5,Abundance~Type) #NS

mycop_LR03_6 <- mycoplasmataceae_ra_family_rl[which(mycoplasmataceae_ra_family_rl$LR_03.06=="Y"),]
wilcox_test(mycop_LR03_6,Abundance~Type) #NS


# MORAXELLA
moraxellaceae_ra_rl <- subset_taxa(randoml_genus, Family=="Moraxellaceae")
moraxellaceae_ra_genus_rl <- tax_glom(moraxellaceae_ra_rl, taxrank = "Genus", NArm = F) %>%
  psmelt()
moraxellaceae_ra_family_rl <- tax_glom(moraxellaceae_ra_rl, taxrank = "Family", NArm = F) %>%
  psmelt()

length(unique(moraxellaceae_ra_genus_rl$Genus))

moraxella_rl_palette <- distinctColorPalette(k=11)


moraxplot_rl <- ggplot(moraxellaceae_ra_genus_rl, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  geom_errorbar(data = moraxellaceae_ra_family_rl, aes(x= MEG_ID, y= Abundance), stat = "summary", width = 0.55) +
  scale_fill_manual(values =moraxella_rl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))
moraxplot_rl

#Comparing Pools to Individuals
#Pools of 12

moraxella_LR06_1 <- moraxellaceae_ra_family_rl[which(moraxellaceae_ra_family_rl$LR_06.01=="Y"),]
wilcox_test(moraxella_LR06_1,Abundance~Type) #NS

moraxella_LR06_1_figure <- moraxellaceae_ra_genus_rl[which(moraxellaceae_ra_genus_rl$LR_06.01=="Y"),]

ggplot(moraxella_LR06_1_figure, aes(x= MEG_ID, y= Abundance)) +
  theme_bw() + 
  #facet_wrap(~Number.in.Pool, nrow = 2, scales = "free_x") +
  labs(y= "Relative Abundance (%)", title = "moraxellaceae") +
  geom_bar(aes(fill = Genus), stat = "summary", colour = "black") +
  scale_fill_manual(values = moraxella_rl_palette) +
  scale_y_continuous(expand = c(0,0.001,0.1,0)) +
  #scale_x_discrete(limits = rh_w_dendro_sample_order, expand = c(0.03,0,0.03,0)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        strip.background = element_rect(fill = "gray44", colour = "black"),
        strip.text = element_text(face = "bold", colour = "white", size = 16),
        plot.title = element_text(size = 22),
        axis.title.y = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 0.95, vjust = 0.95),
        axis.ticks.x = element_line(colour = "black", size = 0.75))


moraxella_LR06_2 <- moraxellaceae_ra_family_rl[which(moraxellaceae_ra_family_rl$LR_06.02=="Y"),]
wilcox_test(moraxella_LR06_2,Abundance~Type) #NS

moraxella_LR06_3 <- moraxellaceae_ra_family_rl[which(moraxellaceae_ra_family_rl$LR_06.03=="Y"),]
wilcox_test(moraxella_LR06_3,Abundance~Type) #NS


moraxella_LR06_4 <- moraxellaceae_ra_family_rl[which(moraxellaceae_ra_family_rl$LR_06.04=="Y"),]
wilcox_test(moraxella_LR06_4,Abundance~Type) #NS

moraxella_LR06_5 <- moraxellaceae_ra_family_rl[which(moraxellaceae_ra_family_rl$LR_06.05=="Y"),]
wilcox_test(moraxella_LR06_5,Abundance~Type) #NS

moraxella_LR06_6 <- moraxellaceae_ra_family_rl[which(moraxellaceae_ra_family_rl$LR_06.06=="Y"),]
wilcox_test(moraxella_LR06_6,Abundance~Type) #NS


###Combine figures
###overall important phyla
majorphyla_combined <- ggarrange(sh_major_phyla_plot, sl_major_phyla_plot, rh_major_phyla_plot, sl_major_phyla_plot, common.legend = T, ncol = 2, nrow = 2, labels = "AUTO")

ggsave("phyla_RA_plot.tiff",plot= majorphyla_combined,device = "tiff", path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/ForPublication/Figures", units = "mm", height=200, width = 180, dpi=600)

majorfamilies_combined <- ggarrange(sh_major_family_plot, sl_major_family_plot, rh_major_family_plot, rl_major_family_plot, ncol = 2, nrow = 2, labels = "AUTO")

ggsave("families_RA_plot.tiff",plot=majorfamilies_combined,device = "tiff", path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/16S Pooling/ForPublication/Figures", units = "mm", height=200, width = 180, dpi=600)

