###BETA DIVERSITY####
ARG_high.css <- phyloseq_transform_css(ARG_high_trimmed, log=F)
ARG_low.css <- phyloseq_transform_css(ARG_low,log=F)

MH_high.css <- phyloseq_transform_css(MH_high, log = F)
MH_low.css <- phyloseq_transform_css(MH_low, log = F)

###Dista
ARG_high.dist <- vegdist(t(otu_table(ARG_high.css)),method = "bray")
ARG_low.dist <- vegdist(t(otu_table(ARG_low.css)),method = "bray")

MH_high.dist <- vegdist(t(otu_table(MH_high.css)),method = "jaccard")
MH_low.dist <- vegdist(t(otu_table(MH_low.css)), method = "jaccard")

high.ord <- vegan::metaMDS(comm = t(otu_table(ARG_high.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)

low.ord <- vegan::metaMDS(comm = t(otu_table(ARG_low.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)

MH_high.ord <-vegan::metaMDS(comm = t(otu_table(MH_high.css)),distance = "jaccard",try = 20, trymax = 50, autotransform = F)
MH_low.ord <-vegan::metaMDS(comm = t(otu_table(MH_low.css)),distance = "jaccard",try = 20, trymax = 50, autotransform = F)



plot_ordination(ARG_high.css,high.ord,  color = "Number.in.Pool") +
  scale_color_manual(values = numberinpoolpalette)+
  stat_ellipse()


plot_ordination(ARG_low.css, low.ord, color="Number.in.Pool")+
  scale_color_manual(values = poolonly_palette)+
  stat_ellipse()

plot_ordination(MH_high.css,MH_high.ord, color = "Number.in.Pool") +
  scale_color_manual(values=numberinpoolpalette)+
  stat_ellipse()

plot_ordination(MH_low.css, MH_low.ord,color = "Number.in.Pool") +
  scale_color_manual(values = poolonly_palette) +
  stat_ellipse()

##Adding centroids
#### findingcenters
#hs
high_plot1 <- ordiplot(high.ord$points)
high_siteslong <-sites.long(high_plot1,ARG_high_div.df)


high_centroids <- envfit(high.ord~ARG_high_div.df$Number.in.Pool)
high_centroids

high_col <- c("1","3","6","12")
#poolonly_col <- c("3","6","12")
high_NMDS_col1 <- c(0.0147,0.1275,0.0968,-0.2667)
high_NMDS_col2 <- c(0.0256,0.0297,-0.1463,-0.0102)
high_centroids.df <-data.frame(high_col,high_NMDS_col1,high_NMDS_col2)


#ls
low_plot1 <- ordiplot(low.ord$points)
low_siteslong <- sites.long(low_plot1,ARG_low_div.df)
low_centroids <- envfit(low.ord~ARG_low_div.df$Number.in.Pool)
low_centroids

low_col <- c("3","6","12")
low_NMDS_col1 <- c(-0.0804,0.1743,-0.0939)
low_NMDS_col2 <- c(-0.1857,0.1876,-0.0019)
low_centroids.df <-data.frame(low_col,low_NMDS_col1,low_NMDS_col2)


## high
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, high Prevalence") +
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
high.adonis <- adonis2(ARG_high.dist~Number.in.Pool, ARG_high_div.df, perm = 9999)
high.adonis #NS, 0.5
write.csv(high.adonis, "TE_highprev_adonisSNV.csv")

high.disper <- betadisper(ARG_high.dist, ARG_high_div.df$Number.in.Pool)
plot(high.disper)
high.permdisp <- permutest(high.disper, permutations = 9999, pairwise = F) #SIG, 0.0001
write.csv(high.permdisp[["pairwise"]][["permuted"]],"highprev_TE_bray_permdisp_SNV.csv")

## high_trimmed
high_ARG_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, high Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = high_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool, fill=Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(18,19,19,19))+
  scale_size_manual(values = c(2,2,2,2), guide = "none") +
  stat_ellipse(data = high_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = high_centroids.df, aes(x=high_NMDS_col1, y=high_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = high_centroids.df, aes(x=high_NMDS_col1, y=high_NMDS_col2, label = high_col), colour = "white", size = 2, fontface = "bold") +
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

ggsave("hs_arg_nmds.tiff", plot = high_ARG_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/RProjects", width = 90, height = 50, units = "mm")




## low_trimmed
low_ARG_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, low_trimmed Prevalence") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous() +
  scale_y_continuous(breaks = c(-1,0,1)) +
  scale_shape_manual(values = c(19,19,19)) +
  geom_point(data = low_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_size_manual(values = c(2,2,2), guide = "none") +
  stat_ellipse(data = low_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = low_centroids.df, aes(x=low_NMDS_col1, y=low_NMDS_col2), fill = poolonly_palette, colour = poolonly_palette, size = 3, shape = c(19,19,19)) +
  geom_text(data = low_centroids.df, aes(x=low_NMDS_col1, y=low_NMDS_col2, label = low_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = poolonly_palette) +
  scale_fill_manual(values = poolonly_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"),
        plot.title=element_blank())

ggsave("ls_arg_nmds.tiff", plot = low_ARG_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/RProjects", width = 90, height = 50, units = "mm")

# stats
low.adonis <- adonis2(ARG_low.dist~Number.in.Pool, ARG_low_div.df,perm = 9999) #NS
low.adonis #NS, P=0.501
write.csv(low.adonis, "TE_lowprev_adonis_SNV.csv")

low.disper <- betadisper(ARG_low.dist, ARG_low_div.df$Number.in.Pool)
plot(low.disper)
low.permdisp <- permutest(low.disper, permutations = 9999, pairwise = F) #No dif, 0.4674
write.csv(low.permdisp[["pairwise"]][["permuted"]],"lowprev_TE_bray_permdisp_SNV.csv")

#### CLUSTERING ON BRAY-CURTIS
high.hclust <- hclust(ARG_high.dist, method = "ward.D2")
plot(high.hclust)
high.dendro <- as.dendrogram(high.hclust)
high.dendro.data <- dendro_data(high.dendro, type = "rectangle")
high.dendro.data

high_ARG_dendro_metadata <- as_tibble(ARG_high.css@sam_data)
high.dendro.data$labels <- high.dendro.data$labels %>%
  left_join(high_ARG_dendro_metadata, by = c("label" = "samplename"))

high_amr_dendro_plot<-ggplot(high.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = high.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 5, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  geom_text(data = high.dendro.data$labels, 
            aes(x=x,y=y, label=HSPool),size=3,position = position_nudge(y=-0.3,x=0),color="black",angle=90)+
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


high_amr_order <- high.dendro.data$labels$label

low.hclust <- hclust(ARG_low.dist, method = "ward.D2")
plot(low.hclust)
low.dendro <- as.dendrogram(low.hclust)
low.dendro.data <- dendro_data(low.dendro, type = "rectangle")
low.dendro.data

low_dendro_metadata <- as_tibble(ARG_low.css@sam_data)
low.dendro.data$labels <- low.dendro.data$labels %>%
  left_join(low_dendro_metadata, by = c("label" = "MEG_ID"))

low_amr_dendro_plot<-ggplot(low.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance", fill="Pool Size") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = low.dendro.data$labels, 
             aes(x=x,y=y, fill= Number.in.Pool),
             size = 5, shape=22, stroke =0.5, position = position_nudge(y=-0.02,x=0), color="black") +
  geom_text(data = low.dendro.data$labels, 
            aes(x=x,y=y, label=HSPool),size=3,position = position_nudge(y=-0.3,x=0),color="black",angle=90)+
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


low_amr_order <- low.dendro.data$labels$label


#MH betadiversity
MH_high.adonis <- adonis2(MH_high.dist~Number.in.Pool, MHhigh.df, perm = 9999)
MH_high.adonis #NS, 0.6336
write.csv(high.adonis, "TE_highprev_adonisSNV.csv")

MH_high.disper <- betadisper(MH_high.dist, MHhigh.df$Number.in.Pool)
plot(MH_high.disper)
MH_high.permdisp <- permutest(MH_high.disper, permutations = 9999, pairwise = F) #NS, 0.6013

MH_low.adonis <- adonis2(MH_low.dist~Number.in.Pool, MHlow.df, perm = 9999)
MH_low.adonis #NS, 0.5537
write.csv(low.adonis, "TE_lowprev_adonisSNV.csv")

MH_low.disper <- betadisper(MH_low.dist, MHlow.df$Number.in.Pool)
plot(MH_low.disper)
MH_low.permdisp <- permutest(MH_low.disper, permutations = 9999, pairwise = F) #NS, 0.8533

####RELATIVE ABUNDANCE



####high####
data_ra <- transform_sample_counts(data.css, function(x) {x/sum(x)}*100)

data_ra_group <- tax_glom(data_ra, taxrank ="Group")
ra_group_melt <- psmelt(data_ra_group)

unique(ra_group_melt$Group)

class_order <- c("Aminocoumarins","Aminoglycosides","Bacitracin","betalactams",
                 "Cationic_antimicrobial_peptides","Glycopeptides","MLS","Multi-drug_resistance",
                 "Mupirocin","Nucleosides","Phenicol","Rifampin","Sulfonamides","Tetracyclines",
                 "Biocide_and_metal_resistance","Drug_and_biocide_and_metal_resistance",
                 "Drug_and_biocide_resistance","Acetate_resistance","Acid_resistance","Multi-biocide_resistance",
                 "Peroxide_resistance","Phenolic_compound_resistance","Quaternary_Ammonium_Compounds_(QACs)_resistance",
                 "Aluminum_resistance","Arsenic_resistance","Copper_resistance","Iron_resistance",
                 "Mercury_resistance","Multi-metal_resistance","Nickel_resistance","Sodium_resistance",
                 "Tungsten_Resistance","Zinc_resistance")
class_label <- c("Aminocoumarins","Aminoglycosides","Bacitracin","Beta-Lactams",
                 "Cationic antimicrobial peptides","Glycopeptides","MLS","Multi-drug resistance",
                 "Mupirocin","Nucleosides","Phenicol","Rifampin","Sulfonamides","Tetracyclines",
                 "Biocide and metal resistance","Drug and biocide and metal resistance",
                 "Drug and biocide resistance","Acetate resistance","Acid resistance","Multi-biocide resistance",
                 "Peroxide resistance","Phenolic compound resistance","Quaternary Ammonium Compounds (QACs) resistance",
                 "Aluminum resistance","Arsenic resistance","Copper resistance","Iron resistance",
                 "Mercury resistance","Multi-metal resistance","Nickel resistance","Sodium resistance",
                 "Tungsten Resistance","Zinc resistance")

rel_abund_high <- transform_sample_counts(ARG_high.css, function(x) {x/sum(x)}*100) #251 "taxa", 29 samples
high_ra_class <- tax_glom(rel_abund_high, taxrank = "Class")
high_ra_class_melt <- psmelt(high_ra_class)
write.csv(high_ra_group_melt,"ra_ARG_high.csv")

length(unique(high_ra_group_melt$Class))
sort(unique(high_ra_group_melt$Class))

palette_33 <- randomColor(count = 33)

#[1] "#27b4c6" "#c5a7ef" "#ea68ed" "#f7ccff" "#98f495" "#f9b6d2" "#17b71d" "#8075d8" "#37cc95" "#68c1c4" "#f2b5d6" "#c96c44" "#70ff6d"
#[14] "#8820bc" "#f3ffaf" "#fff766" "#c0f3f7" "#e2b68a" "#bef49c" "#db4376" "#bded76" "#e8c609" "#17b58d" "#baebf4" "#7f31a5" "#0108c1"
#[27] "#b075f4" "#9e85e2" "#48c988" "#e08170" "#f6fc58" "#44c486" "#f75980"


ra_class_high_plot <- ggplot(high_ra_class_melt, aes(x= samplename, y= Abundance, fill= Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = high_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(limits=class_order, values = palette_33, labels=class_label) +
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
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8, face ="bold"),
        axis.text.x = element_blank())

HP_63 <- high_ra_group_melt[which(high_ra_group_melt$HP_06.03=="Y"),]

ggplot(HP_63, aes(x=samplename, y=Abundance, fill=Class))+
  theme_minimal()+
  geom_bar(stat="identity",color="gray")+
  scale_fill_manual(values=randomColor(40))

stackeddendro_class_high_plot <- ggarrange(high_amr_dendro_plot, ra_class_high_plot, ncol=1, heights = c(1.5,2))
##Tetracyclines
high_tetracylines <- subset_taxa(high_ra_group, Class=="Tetracyclines")
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
  #scale_fill_manual(values=c("#8048a5","#2c0845","#c076db","#de9be3","#792b7c","#df73d1","#681552"), limits=c("TET16S","TET40","TETH","TETQ","TETR","TETW","TETZ")) +
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


ggsave("stackeddendroclass_highplot.tiff", path = "../DATA", plot =stackeddendro_class_high_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)
ggsave("grouphigh_plot.tiff",path="~/Desktop",plot=ra_group_high_plot, device="tiff",dpi=600)

ra_high_class <- tax_glom(rel_abund_high, taxrank = "Class") #9 taxa
ra_high_class_melt <- psmelt(ra_high_class)

unique(ra_high_class_melt$Class)

class_order <- c("Tetracyclines","Sulfonamides","Aminoglycosides","MLS","Phenicol","betalactams","Elfamycins","Drug_and_biocide_resistance","Mercury_resistance")

high_mean_sd <- ra_high_class_melt %>%
  group_by(Number.in.Pool,Class)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))


high_class_plot <-ggplot(high_ra_class_melt, aes(x= Class, y= Abundance, fill= Number.in.Pool, alpha=Number.in.Pool)) +
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

####low####

rel_abund_low <- transform_sample_counts(ARG_low.css, function(x) {x/sum(x)}*100) #251 "taxa", 29 samples
#rel_abund_low_trimmed <- transform_sample_counts(low_trimmed.css, function(x) {x/sum(x)}*100) #265 taxa, 17 samples

# agglomerate at different levels
rel_abund_group_low <- tax_glom(rel_abund_low, taxrank = "Group") # 63 groups
ra_group_melt_low <- psmelt(rel_abund_group_low)
ra_class_low <- tax_glom(rel_abund_low,taxrank = "Class")
ra_class_melt_low <- psmelt(ra_class_low)

#ra_low_trimmed_group <- tax_glom(rel_abund_low_trimmed, taxrank = "Group") #57 groups
#ra_low_trimmed_group_melt <- psmelt(ra_low_trimmed_group)

length(unique(ra_class_melt_low$Class))
sort(unique(ra_class_melt_low$Class))
ra_group_palette_low <- distinctColorPalette(57)#--done, can redo if don't like
#ra_group_50_low <- prune_taxa(names(sort(taxa_sums(rel_abund_group_low),T)[1:50]), rel_abund_group_low)

write.csv(ra_group_melt_low,"ra_ARG_low.csv")

write.csv(ra_low_trimmed_group_melt,"lowgroups.csv")



write.csv(otu_table(rel_abund_class_low),"low_class_otus.csv")
write.csv(tax_table(rel_abund_class_low),"low_class_taxa.csv")

abundance_group_low_df <- ra_low_trimmed_group_melt %>%
  group_by(Class,Group)%>%
  reframe(mean=mean(Abundance))

View(abundance_group_low_df)



low_group_order <- c("A16S","ANT3-DPRIME","ANT6","ANT9","APH2-DPRIME","APH3-DPRIME","APH3-PRIME","APH6","RRS","RRSA","RRSC","RRSH",
                     "ROB","CMX","FLOR",
                     "CFR","ERMA","ERMB","ERMF","ERMX","LNUA","LNUC","MEFA","MEFE","MLS23S","MPHE","MSRE","MYRA","VATE",
                     "FOLP","SULI","SULII",
                     "TET32","TET33","TET40","TET44","TETH","TETM","TETO","TETQ","TETR","TETW","TETX","TETY",
                     "ASR",
                     "TCRA",
                     "TCRY",
                     "TCRZ",
                     "QACE",
                     "QACEDELTA1",
                     "TUFAB",
                     "BLE",
                     "BRP",
                     "CFRA",
                     "EMRE",
                     "OPTRA",
                     "O23S")

all_group_palette

"#cb3f54" "#df6f70" "#8c282f" "#d94d4d" "#90281f" "#d0735a" "#b84828" "#e08053" "#cf6c2f" "#7d3f13" "#d38725" "#ca9457" "#e16fbc" "#b53c8f" "#bb689c"
[16] "#ace268" "#468027" "#a1dd8a" "#65b652" "#296021" "#569c56" "#5bd581" "#399f69" "#6fdfa2" "#48f2a7" "#37a886" "#55e5c1" "#3be6ea" "#5999e0" "#4663ab"
[31] "#a27521" "#dcab56" "#ceaa34" "#8048a5" "#633c7c" "#420e66" "#2c0845" "#3e1957" "#c076db" "#a068b1" "#52115d" "#de9be3" "#792b7c" "#df73d1" "#a73a9a"
[46] "#8a3577" "#681552" "#c5c575" "#b6c142" "#828f29" "#c1d168" "#4663ab" "#5b7fee" "#8587de" "#373e89" "#182465" "#5356bb" "#251d77" "#9778e2" "#1c0b53"
[61] "#311c66" "#240568"

low_group_palette <- c("#cb3f54","#df6f70","#8c282f","#d94d4d","#90281f","#d0735a","#b84828","#e08053","#cf6c2f","#7d3f13","#d38725","#ca9457",
                       "#e16fbc","#b53c8f","#bb689c",
                       "#ace268","#a1dd8a","#65b652","#296021","#569c56","#5bd581","#399f69","#6fdfa2","#48f2a7","#37a886","#55e5c1","#3be6ea","#5999e0","#4663ab",
                       "#a27521","#dcab56","#ceaa34",
                       "#633c7c","#420e66","#2c0845","#3e1957","#c076db","#a068b1","#52115d","#de9be3","#792b7c","#df73d1","#a73a9a","#8a3577",
                       "#c5c575","#b6c142","#828f29","#c1d168","#4663ab","#5b7fee","#8587de","#373e89","#182465","#9778e2","#1c0b53","#311c66","#240568")

ra_low_trimmed_group_melt$Group <-factor(ra_low_trimmed_group_melt$Group, levels=low_group_order)

low_arg_abundance <- ra_low_trimmed_group_melt %>%
  group_by(Group, Number.in.Pool)%>%
  reframe(mean=mean(Abundance),sd=sd(Abundance))

write.csv(low_arg_abundance, "low_arg_abundance.csv")

ra_group_low_plot <- ggplot(ra_class_melt_low, aes(x= samplename, y= Abundance, fill= Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = low_amr_order, expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(limits=class_order, values = palette_33, labels=class_label) +
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
        legend.title = element_text(size=10),
        legend.text = element_text(size = 8, face ="bold"),
        axis.text.x = element_blank())

stackeddendro_group_low_trimmed_plot <- ggarrange(low_amr_dendro_plot, ra_group_low_plot, ncol=1,heights = c(30,65))


combined_dendros_ra <- ggarrange(stackeddendro_class_high_plot,stackeddendro_group_low_trimmed_plot, ncol=1)

ggsave("stackeddendrogroup_lowplot.tiff", path = "../DATA", plot =stackeddendro_group_low_trimmed_plot,device = "tiff", dpi =600, units = "mm", width = 180, height =  100)
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

######MH#####
any(taxa_sums(MH_high)==0)
any(taxa_sums(MH_low)==0)
any(taxa_sums(MH_indiv)==0)

MH_high.css <- phyloseq_transform_css(MH_high, log=F)
MH_low.css <- phyloseq_transform_css(MH_low,log=F)
MH_indiv.css <- phyloseq_transform_css(MH_indiv, log=F)

MH_high_rel_abun <- transform_sample_counts(MH_high.css, function(x) {x/sum(x)}*100)
MH_low_rel_abun <- transform_sample_counts(MH_low.css, function(x) {x/sum(x)}*100)
MH_indiv_rel_abun <- transform_sample_counts(MH_indiv.css, function(x) {x/sum(x)}*100)

MH_high_melt <- psmelt(MH_high_rel_abun)
MH_low_melt <- psmelt(MH_low_rel_abun)
MH_indiv_melt <- psmelt(MH_indiv_rel_abun)


####BETA DIVERSITY####
###Distance
MH_high.dist <- vegdist(t(otu_table(MH_high.css)),method = "bray")
MH_low.dist <- vegdist(t(otu_table(MH_low.css)),method = "bray")
MH_indiv.dist <- vegdist(t(otu_table(MH_indiv.css)),method = "bray")

MH_high.ord <- vegan::metaMDS(comm = t(otu_table(MH_high.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)
MH_low.ord <- vegan::metaMDS(comm = t(otu_table(MH_low.css)), distance = "bray", try = 20, trymax = 50, autotransform = F)
MH_indiv.ord <- vegan::metaMDS(comm=t(otu_table(MH_indiv.css)),distance = "bray", try = 20, trymax = 50, autotransform = F)

plot_ordination(MH_high.css,MH_high.ord,  color = "Number.in.Pool") +
  scale_color_manual(values = numberinpoolpalette)+
  stat_ellipse()
plot_ordination(MH_low.css,MH_low.ord,  color = "Number.in.Pool") +
  scale_color_manual(values = poolonly_palette)+
  stat_ellipse()
plot_ordination(MH_indiv.css,MH_indiv.ord,  color = "Mh.Cult") +
  stat_ellipse()

##Adding centroids
#### findingcenters
MH_high_plot1 <- ordiplot(MH_high.ord$points)
MH_high_siteslong <-sites.long(MH_high_plot1,MHhigh.df)

MH_low_plot1 <- ordiplot(MH_low.ord$points)
MH_low_siteslong <-sites.long(MH_low_plot1,MHlow.df)

MH_indiv_plot1 <- ordiplot(MH_indiv.ord$points)
MH_indiv_siteslong <-sites.long(MH_indiv_plot1,MHindiv.df)

MH_high_centroids <- envfit(MH_high.ord~MHhigh.df$Number.in.Pool)
MH_high_centroids

high_col <- c("1","3","6","12")
#poolonly_col <- c("3","6","12")
MH_high_NMDS_col1 <- c(0.0020,-0.1802,-0.1965,0.3361)
MH_high_NMDS_col2 <- c(-0.0172,0.0826,0.0108,-0.0229)
MH_high_centroids.df <-data.frame(high_col,MH_high_NMDS_col1,MH_high_NMDS_col2)


MH_low_centroids <- envfit(MH_low.ord~MHlow.df$Number.in.Pool)
MH_low_centroids

low_col <- c("3","6","12")
MH_low_NMDS_col1 <- c(-0.1998,-0.2533,0.4531)
MH_low_NMDS_col2 <- c(-0.0597,0.0515,0.0081)
MH_low_centroids.df <-data.frame(low_col,MH_low_NMDS_col1,MH_low_NMDS_col2)

MH_indiv_centroids <- envfit(MH_indiv.ord~MHindiv.df$Mh.Cult)
MH_indiv_centroids

indiv_col <- c("N","Y")
#poolonly_col <- c("3","6","12")
MH_indiv_NMDS_col1 <- c(0.0537,-0.0455)
MH_indiv_NMDS_col2 <- c(-0.0380,0.0321)
MH_indiv_centroids.df <-data.frame(indiv_col,MH_indiv_NMDS_col1,MH_indiv_NMDS_col2)

## high
MH_high_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "MH Beta Diversity, High") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = MH_high_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool, fill=Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(18,19,19,19))+
  scale_size_manual(values = c(5,5,5,5), guide = "none") +
  stat_ellipse(data = MH_high_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = MH_high_centroids.df, aes(x=MH_high_NMDS_col1, y=MH_high_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(18,19,19,19)) +
  geom_text(data = MH_high_centroids.df, aes(x=MH_high_NMDS_col1, y=MH_high_NMDS_col2, label = high_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", size = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, colour = "black"))


# stats
MH_high.adonis <- adonis2(MH_high.dist~Number.in.Pool, MHhigh.df, perm = 9999)
MH_high.adonis #NS, 0.177
write.csv(MH_high.adonis, "PSV_highprev_adonis.csv")

MH_high.disper <- betadisper(MH_high.dist, MHhigh.df$Number.in.Pool)
plot(MH_high.disper)
MH_high.permdisp <- permutest(MH_high.disper, permutations = 9999, pairwise = F) #SIG, 0.0088
write.csv(MH_high.permdisp[["pairwise"]][["permuted"]],"highprev_PSV_bray_permdisp.csv")

ggsave("MH_high_nmds.tiff", plot = MH_high_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/MH", width = 90, height = 50, units = "mm")


## low_trimmed
MH_low_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "MH Beta Diversity, low") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = MH_low_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool, fill=Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(19,19,19))+
  scale_size_manual(values = c(5,5,5), guide = "none") +
  stat_ellipse(data = MH_low_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = MH_low_centroids.df, aes(x=MH_low_NMDS_col1, y=MH_low_NMDS_col2), fill = poolonly_palette, colour = poolonly_palette, size = 10, shape = 19) +
  geom_text(data = MH_low_centroids.df, aes(x=MH_low_NMDS_col1, y=MH_low_NMDS_col2, label = low_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = poolonly_palette) +
  scale_fill_manual(values = poolonly_palette) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", size = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, colour = "black"))

ggsave("MH_low_nmds.tiff", plot = MH_low_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/MH", width = 90, height = 50, units = "mm")

# stats
MH_low.adonis <- adonis2(MH_low.dist~Number.in.Pool, MHlow.df, perm = 9999)
MH_low.adonis #NS, 0.0784
write.csv(MH_low.adonis, "PSV_lowprev_adonis.csv")

MH_low.disper <- betadisper(MH_low.dist, MHlow.df$Number.in.Pool)
plot(MH_low.disper)
MH_low.permdisp <- permutest(MH_low.disper, permutations = 9999, pairwise = F) #NS, 0.146
write.csv(MH_low.permdisp[["pairwise"]][["permuted"]],"lowprev_PSV_bray_permdisp.csv")

## indiv
MH_indiv_nmds_plot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "MH Beta Diversity, indiv") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = MH_indiv_siteslong, aes(x=axis1,y=axis2, colour= Mh.Cult, shape = Mh.Cult, alpha = Mh.Cult, size = Mh.Cult, fill=Mh.Cult)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = "none") +
  scale_shape_manual(values=c(19,19))+
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = MH_indiv_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Mh.Cult, fill = Mh.Cult), alpha = c(0.1), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = MH_indiv_centroids.df, aes(x=MH_indiv_NMDS_col1, y=MH_indiv_NMDS_col2), size = 10, shape = 19) +
  geom_text(data = MH_indiv_centroids.df, aes(x=MH_indiv_NMDS_col1, y=MH_indiv_NMDS_col2, label = indiv_col), colour = "white", size = 6, fontface = "bold") +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", size = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, colour = "black"))

ggsave("MH_indiv_nmds.tiff", plot = MH_indiv_nmds_plot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/DATA/MH", width = 90, height = 50, units = "mm")

# stats
MH_indiv.adonis <- adonis2(MH_indiv.dist~Mh.Cult, MHindiv.df, perm = 9999)
MH_indiv.adonis #Sig, P=0.0017
write.csv(MH_indiv.adonis, "PSV_indivprev_adonis.csv")

MH_indiv.disper <- betadisper(MH_indiv.dist, MHindiv.df$Mh.Cult)
plot(MH_indiv.disper)
MH_indiv.permdisp <- permutest(MH_indiv.disper, permutations = 9999, pairwise = F) #Sig, P=0.0265
write.csv(MH_indiv.permdisp[["pairwise"]][["permuted"]],"indivprev_PSV_bray_permdisp.csv")


##PSV script from Enrique
MH_high_sample_factor_by_abund <- MH_high_melt %>%
  dplyr::group_by(PSV) %>%
  dplyr::summarize(median_PSV_perc = median(Abundance)) %>%
  arrange(-median_PSV_perc)

# find Phyla whose rel. abund. is less than 1%
remainder_MHhigh <- MH_high_sample_factor_by_abund[MH_high_sample_factor_by_abund$median_PSV_perc <= 0.5,]$PSV

# change the name of taxa to whose rel. abund. is less than 0.5% to "indiv abundance phyla (<1%)"
MH_high_melt[MH_high_melt$PSV %in% remainder_MHhigh,]$PSV <- 'indiv abundance PSV (<0.5%)'

MH_high_melt_agg <- MH_high_melt %>%
  group_by(Sample, PSV) %>%  # Group by both Sample and PSV
  summarise(Abundance = sum(Abundance)) %>%  # Sum the Abundance for each group
  ungroup()


# Determine the number of taxa
length(unique(MH_high_melt_agg$PSV))

# Cluster using the function hclust() and default settings
MH_high.dist <- vegdist(t(otu_table(MH_high.css)), method = "bray")
MH_high.hclust <- hclust(MH_high.dist)
plot(MH_high.hclust) # example plot

# Extract data as dendrogram
MH_high.dendro <- as.dendrogram(MH_high.hclust)
MH_high.dendro.data <- dendro_data(MH_high.dendro, type = "rectangle")
MH_high.dendro.data #  this object contains $segments and $labels

# Sample names in order based on clustering
sample_names_MH_high <- MH_high.dendro.data$labels$label


# Add metadata #### 
# make sure that the "sample_id" column (or change name) in your phyloseq object matches with this column: ps_mhpp.dendro.data$labels
MH_high.dendro.data$labels <- left_join(MH_high.dendro.data$labels, sample_data(MH_high.css), by = c("label" = "samplename"))

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_MH_high <- with(
  segment(MH_high.dendro.data),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table_MH_high <- with(
  MH_high.dendro.data$labels, 
  data.frame(y_center = x, gene = as.character(label), x = y ,mh_percentage = as.character(PSV_per) , height = 1))


# Table to position the samples
sample_pos_table_MH_high <- data.frame(sample = sample_names_MH_high) %>%
  dplyr::mutate(x_center = (1:dplyr::n()),  width = 1)

##
######## Relative abundance bar plot #########
##

# Use class melted data and add gene locations
joined_ra_PSV_MH_high <-  MH_high_melt_agg %>%
  left_join(gene_pos_table_MH_high, by = c("Sample" = "gene")) %>%
  left_join(sample_pos_table_MH_high, by = c("Sample" = "sample")) 

# Calculate the mean relative abundance of each class taxa, sort by most abundant to least
factor_by_abund_MH_high <- rel_abund_melt %>%
  dplyr::group_by(PSV) %>%
  dplyr::summarize(median_PSV = median(Abundance)) %>%
  arrange(-median_PSV)



plt_rel_PSV_MH_high <- ggplot(joined_ra_PSV_MH_high, 
                      aes(x = x_center, y = Abundance, fill = PSV, 
                          height = height, width = width)) + 
  coord_flip() +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values = col_vector) + #use this if not using color palette from ggthemes of ggsci
  scale_fill_tableau(palette = "Tableau 20")+
  scale_x_continuous(breaks = sample_pos_table_MH_high$x_center, 
                     labels = sample_pos_table_MH_high$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  #scale_y_continuous(breaks = gene_pos_table[, "y_center"], labels = rep("", nrow(gene_pos_table)),limits = gene_axis_limits, expand = c(0, 0)) + 
  labs(x = "", y = "Relative abundance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.75),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.75, lineend = "square", colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())
plt_rel_PSV_MH_high

# Dendrogram plot ####
plt_dendr_PSV_MH_high <- ggplot(segment_data_MH_high) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=.75, lineend="round", linejoin="round") +
  #geom_point(data=gene_pos_table_ps_mhpp, aes(x, y_center, colour=Level, shape=Level),
             #size=8, stroke=0.75, position=position_nudge(x=-0.03)) +
  scale_y_continuous(breaks=gene_pos_table_MH_high$y_center, 
                     labels=gene_pos_table_MH_high$gene, 
                     #limits=ge, 
                     expand=c(0, 0)) + 
  labs(x="Ward's Distance", y="", title="") +
  scale_x_reverse() + 
  #scale_colour_manual(name="16S Mh level", values = c("yellow","light gray", "#6497B1", "#005B96", "#03396C", "#011F4B")) +  # Set legend title and colors
  #scale_shape_manual(values=c(15,15,15,15,15,15), guide="none") +  
  theme_bw() + 
  theme(legend.position="right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        panel.border=element_blank(),
        plot.title=element_text(size=30),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(linewidth=0.75),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12, colour="black"))
plt_dendr_PSV_MH_high

## Make plot showing relative abundance of PSV hits to nonhost counts ####

# Add a column for the remaining percentage to make up 100%
mhpp_PSV.ps.css.df$remainder_percentage <- 100 - mhpp_PSV.ps.css.df$mh_percentage

# Reshape data to long format for ggplot2
long_data <- tidyr::pivot_longer(mhpp_PSV.ps.css.df, cols = c("mh_percentage", "remainder_percentage"), 
                                 names_to = "category", values_to = "value")

long_data$category <- factor(long_data$category, levels = c("remainder_percentage","mh_percentage"))

long_data$sample_id <- factor(long_data$sample_id, levels = as.factor(gene_pos_table_ps_mhpp$gene))

