#############################################################################################
##############################         BETA DIVERSITY         ##############################
#############################################################################################
#############################################################################################
##double check no samples with 0 reads
any(sample_sums(data_50k)==0)
any(sample_sums(subsethigh)==0)
any(sample_sums(subsetlow)==0)
any(sample_sums(randomhigh)==0)
any(sample_sums(randomlow)==0)
any(sample_sums(pools_3v6high)==0)
any(sample_sums(pools_3v6low)==0)
any(sample_sums(individual_hl)==0)
any(sample_sums(hs_pools)==0)
any(sample_sums(ls_pools)==0)
any(sample_sums(hr_pools)==0)
any(sample_sums(lr_pools)==0)
any(sample_sums(high_all)==0)
any(sample_sums(low_all)==0)
any(sample_sums(high_pools)==0)
any(sample_sums(low_pools)==0)
##None

##double check taxa sums
any(taxa_sums(subsethigh)==0)
any(taxa_sums(subsetlow)==0)
any(taxa_sums(randomhigh)==0)
any(taxa_sums(randomlow)==0)
any(taxa_sums(pools_3v6high)==0)
any(taxa_sums(pools_3v6low)==0)
any(taxa_sums(individual_hl)==0)
any(taxa_sums(hs_pools)==0)
any(taxa_sums(ls_pools)==0)
any(taxa_sums(hr_pools)==0)
any(taxa_sums(lr_pools)==0)
any(taxa_sums(high_all)==0)
any(taxa_sums(low_all)==0)
any(taxa_sums(high_pools)==0)
any(taxa_sums(low_pools)==0)
#also false


####################################   CSS TRANSFORM    ###################################
#overall
data_16S.css <- phyloseq_transform_css(data_16S_50k, log = F)
data_16S.css
any(taxa_sums(data_16S.css)==0)
data_16S.css.ra <- transform_sample_counts(data_16S.css, function(x) {x/sum(x)} * 100)

data_16S.css.df <- as(sample_data(data.css), "data.frame")

#all pools
##high
high.css <- phyloseq_transform_css(high_all,log = F)
high.css
high.css.df <- as(sample_data(high.css),"data.frame")

##low
low.css <- phyloseq_transform_css(low_all,log = F)
low.css
low.css.df <- as(sample_data(low.css),"data.frame")



#subset
subset.css  <-phyloseq_transform_css(subsetpooling)

sh_16S.css <- phyloseq_transform_css(HSpool_16S.ps, log = F)
sh_16S.css.ra <- transform_sample_counts(sh_16S.css, function(x) {x/sum(x)}*100)
sh_16s.css.df <- as(sample_data(sh_16S.css),"data.frame")

sl_16S.css <- phyloseq_transform_css(LSpool_16S.ps, log = F)
sl_16S.css.ra <- transform_sample_counts(sl_16S.css, function(x) {x/sum(x)}*100)
sl_16s.css.df <- as(sample_data(sl_16S.css),"data.frame")


##high
subseth.css <- phyloseq_transform_css(subsethigh,log = F)
subseth.css
subseth.css.ra <- transform_sample_counts(subseth.css, function(x) {x/sum(x)}*100)
subseth.css.df <- as(sample_data(subseth.css),"data.frame")

##low
subsetl.css <- phyloseq_transform_css(subsetlow,log = F)
subsetl.css
subsetl.css.ra <- transform_sample_counts(subsetl.css, function(x) {x/sum(x)}*100)
subsetl.css.df <- as(sample_data(subsetl.css),"data.frame")
sample_data(subsetl.css)

#random
##high
HR_16S.css <- phyloseq_transform_css(HRpool_16S.ps,log = F)
HR_16S.css.ra <-transform_sample_counts(HR_16S.css, function(x) {x/sum(x)}*100)
HR_16S.css.df <- as(sample_data(HR_16S.css), "data.frame")

LR_16S.css <- phyloseq_transform_css(LRpool_16S.ps,log = F)
LR_16S.css.ra <-transform_sample_counts(LR_16S.css, function(x) {x/sum(x)}*100)
LR_16S.css.df <- as(sample_data(LR_16S.css), "data.frame")

randomh.css.df <- as(sample_data(randomh.css),"data.frame")

##low
randoml.css <- phyloseq_transform_css(randomlow,log = F)
randoml.css

randoml.css.df <- as(sample_data(randoml.css),"data.frame")

#pools3and6
##high
pools_3v6h.css <- phyloseq_transform_css(pools_3v6high,log = F)
pools_3v6h.css
pools_3v6h.css.ra <- transform_sample_counts(pools_3v6h.css, function(x) {x/sum(x)}*100)
pools_3v6h.css.df <- as(sample_data(pools_3v6h.css),"data.frame")

##low
pools_3v6l.css <- phyloseq_transform_css(pools_3v6low,log = F)
pools_3v6l.css
pools_3v6l.css.ra <- transform_sample_counts(pools_3v6l.css, function(x) {x/sum(x)}*100)
pools_3v6l.css.df <- as(sample_data(pools_3v6l.css),"data.frame")


## calculate UniFrac dists
#all
all_wunifrac.dist <- wunifrac(data.css)
all_gunifrac.dist <- gunifrac(data.css)
all_uwunifrac.dist <- uwunifrac(data.css)

#pools (54 unique individuals in high pools of 6 (at most), should be a way to count this, I did this in Excel)
#.    (31 unique individuals in high pools of 3)
#high
high_wunifrac.dist <- wunifrac(high.css)
high_gunifrac.dist <- gunifrac(high.css)
high_uwunifrac.dist <- uwunifrac(high.css)

#low. 51 low pools 6, 31 in 3 
low_wunifrac.dist <- wunifrac(low.css)
low_gunifrac.dist <- gunifrac(low.css)
low_uwunifrac.dist <- uwunifrac(low.css)


#subsetpooling
#high
subseth_wunifrac.dist <- wunifrac(subseth.css)
subseth_gunifrac.dist <- gunifrac(subseth.css)
subseth_uwunifrac.dist <- uwunifrac(subseth.css)

#low
subsetl_wunifrac.dist <- wunifrac(subsetl.css)
subsetl_gunifrac.dist <- gunifrac(subsetl.css)
subsetl_uwunifrac.dist <- uwunifrac(subsetl.css)

#randompooling
#high
randomh_wunifrac.dist <- wunifrac(randomh.css)
randomh_gunifrac.dist <- gunifrac(randomh.css)
randomh_uwunifrac.dist <- uwunifrac(randomh.css)

#low
randoml_wunifrac.dist <- wunifrac(randoml.css)
randoml_gunifrac.dist <- gunifrac(randoml.css)
randoml_uwunifrac.dist <- uwunifrac(randoml.css)

#pools_3v6pooling
#high
pools_3v6h_wunifrac.dist <- wunifrac(pools_3v6h.css)
pools_3v6h_gunifrac.dist <- gunifrac(pools_3v6h.css)
pools_3v6h_uwunifrac.dist <- uwunifrac(pools_3v6h.css)

#low
pools_3v6l_wunifrac.dist <- wunifrac(pools_3v6l.css)
pools_3v6l_gunifrac.dist <- gunifrac(pools_3v6l.css)
pools_3v6l_uwunifrac.dist <- uwunifrac(pools_3v6l.css)



# ordinate
all_wunifrac.ord <- ordinate(data.css, method = "NMDS", distance = all_wunifrac.dist)
all_gunifrac.ord <- ordinate(data.css, method = "NMDS", distance = all_gunifrac.dist)
all_uwunifrac.ord <- ordinate(data.css, method = "NMDS", distance = all_uwunifrac.dist)

#high
high_wunifrac.ord <- ordinate(high.css, method = "NMDS", distance = high_wunifrac.dist)
high_gunifrac.ord <- ordinate(high.css, method = "NMDS", distance = high_gunifrac.dist)
high_uwunifrac.ord <- ordinate(high.css, method = "NMDS", distance = high_uwunifrac.dist)

#low
low_wunifrac.ord <- ordinate(low.css, method = "NMDS", distance = low_wunifrac.dist)
low_gunifrac.ord <- ordinate(low.css, method = "NMDS", distance = low_gunifrac.dist)
low_uwunifrac.ord <- ordinate(low.css, method = "NMDS", distance = low_uwunifrac.dist)

#subsetpooling
#high
subseth_wunifrac.ord <- ordinate(subseth.css,method = "NMDS",distance = subseth_wunifrac.dist)
subseth_gunifrac.ord <- ordinate(subseth.css,method = "NMDS",distance = subseth_gunifrac.dist)
subseth_uwunifrac.ord <- ordinate(subseth.css,method = "NMDS",distance = subseth_uwunifrac.dist)

#low
subsetl_wunifrac.ord <- ordinate(subsetl.css,method = "NMDS",distance = subsetl_wunifrac.dist)
subsetl_gunifrac.ord <- ordinate(subsetl.css,method = "NMDS",distance = subsetl_gunifrac.dist)
subsetl_uwunifrac.ord <- ordinate(subsetl.css,method = "NMDS",distance = subsetl_uwunifrac.dist)

#randompooling
#high
randomh_wunifrac.ord <- ordinate(randomh.css,method = "NMDS",distance = randomh_wunifrac.dist)
randomh_gunifrac.ord <- ordinate(randomh.css,method = "NMDS",distance = randomh_gunifrac.dist)
randomh_uwunifrac.ord <- ordinate(randomh.css,method = "NMDS",distance = randomh_uwunifrac.dist)

#low
randoml_wunifrac.ord <- ordinate(randoml.css,method = "NMDS",distance = randoml_wunifrac.dist)
randoml_gunifrac.ord <- ordinate(randoml.css,method = "NMDS",distance = randoml_gunifrac.dist)
randoml_uwunifrac.ord <- ordinate(randoml.css,method = "NMDS",distance = randoml_uwunifrac.dist)

#pools_3v6pooling
#high
pools_3v6h_wunifrac.ord <- ordinate(pools_3v6h.css,method = "NMDS",distance = pools_3v6h_wunifrac.dist)
pools_3v6h_gunifrac.ord <- ordinate(pools_3v6h.css,method = "NMDS",distance = pools_3v6h_gunifrac.dist)
pools_3v6h_uwunifrac.ord <- ordinate(pools_3v6h.css,method = "NMDS",distance = pools_3v6h_uwunifrac.dist)

#low
pools_3v6l_wunifrac.ord <- ordinate(pools_3v6l.css,method = "NMDS",distance = pools_3v6l_wunifrac.dist)
pools_3v6l_gunifrac.ord <- ordinate(pools_3v6l.css,method = "NMDS",distance = pools_3v6l_gunifrac.dist)
pools_3v6l_uwunifrac.ord <- ordinate(pools_3v6l.css,method = "NMDS",distance = pools_3v6l_uwunifrac.dist)

#PCoA-wunifrac
subseth_wunifrac.pcoa <- ordinate(subseth.css,method = "PCoA",distance = subseth_wunifrac.dist)
subsetl_wunifrac.pcoa <- ordinate(subsetl.css,method = "PCoA",distance = subsetl_wunifrac.dist)
randomh_wunifrac.pcoa <- ordinate(randomh.css,method = "PCoA",distance = randomh_wunifrac.dist)
randoml_wunifrac.pcoa <- ordinate(randoml.css,method = "PCoA",distance = randoml_wunifrac.dist)

#PCoA-uwunifrac
subseth_uwunifrac.pcoa <- ordinate(subseth.css,method = "PCoA",distance = subseth_uwunifrac.dist)
subsetl_uwunifrac.pcoa <- ordinate(subsetl.css,method = "PCoA",distance = subsetl_uwunifrac.dist)
randomh_uwunifrac.pcoa <- ordinate(randomh.css,method = "PCoA",distance = randomh_uwunifrac.dist)
randoml_uwunifrac.pcoa <- ordinate(randoml.css,method = "PCoA",distance = randoml_uwunifrac.dist)


# plot weighted unifrac

pool_col <- c("1","3","6","12")
Number.in.Pool <- c("1","3","6","12")

#### findingcenters
#high
high_wunifrac_plot1 <-ordiplot(high_wunifrac.ord$points)
high_wunifrac_siteslong <- sites.long(high_wunifrac_plot1,high.css.df)
high_wunifrac_centroids <- envfit(high_wunifrac.ord~high.css.df$Number.in.Pool)
high_wunifrac_centroids

high_wunifrac_NMDS_col1 <- c(0.0032,-0.0082,-0.0154,0.0003)
high_wunifrac_NMDS_col2 <- c(0.0022,-0.0074,-0.0079,-0.0008)
high_wunifrac_centroids_df <-data.frame(pool_col, high_wunifrac_NMDS_col1, high_wunifrac_NMDS_col2)


#hs

subseth_wunifrac_plot1 <- ordiplot(subseth_wunifrac.ord$points)
subseth_wunifrac_pcoaplot <- ordiplot(subseth_wunifrac.pcoa$vectors)
subseth_wunifrac_siteslong <-sites.long(subseth_wunifrac_plot1,subseth.css.df)
subseth_wunifrac_pcoasites <- sites.long(subseth_wunifrac_pcoaplot, subseth.css.df)


subseth_wunifrac_centroids <- envfit(subseth_wunifrac.ord~subseth.css.df$Number.in.Pool)
subseth_wunifrac_centroids
subseth_wunifrac_pcoa_centroids <- envfit(subseth_wunifrac.pcoa$vectors~subseth.css.df$Number.in.Pool)
subseth_wunifrac_pcoa_centroids


subseth_wunifrac_NMDS_col1 <- c(0.0027,-0.0101,-0.0258,-0.0006)
subseth_wunifrac_NMDS_col2 <- c(0.0011,0.0002,-0.0135,-0.0016)
subseth_wunifrac_centroids.df <-data.frame(subseth_wunifrac_col,subseth_wunifrac_NMDS_col1,subseth_wunifrac_NMDS_col2)

subseth_wunifrac_PCoA_col1 <- c(-0.0041, 0.0100, 0.0421, 0.0039)
subseth_wunifrac_PCoA_col2 <- c(-0.0008, 0.0138, 0.0059, -0.0090)
subseth_wunifrac_PCoA_centroids.df <- data.frame(Number.in.Pool, subseth_wunifrac_PCoA_col1, subseth_wunifrac_PCoA_col2)


##ls
subsetl_wunifrac_plot1 <- ordiplot(subsetl_wunifrac.ord$points)
subsetl_wunifrac_siteslong <- sites.long(subsetl_wunifrac_plot1,subsetl.css.df)
subsetl.css.df$Number.in.Pool <- as.factor(subsetl.css.df$Number.in.Pool)
subsetl_wunifrac_centroids <- envfit(subsetl_wunifrac.ord~subsetl.css.df$Number.in.Pool)
subsetl_wunifrac_centroids

subsetl_wunifrac_NMDS_col1 <- c(0.0015,0.027,-0.0033,-0.0353)
subsetl_wunifrac_NMDS_col2 <- c(0.0007,0.0322,0.0101,-0.0508)
subsetl_wunifrac_centroids.df <-data.frame(pool_col,subsetl_wunifrac_NMDS_col1,subsetl_wunifrac_NMDS_col2)

##hr
randomh.css.df$Number.in.Pool <- as.factor(randomh.css.df$Number.in.Pool)
randomh_wunifrac_plot1 <- ordiplot(randomh_wunifrac.ord$points)
randomh_wunifrac_siteslong <- sites.long(randomh_wunifrac_plot1,randomh.css.df)

randomh_wunifrac_centroids <- envfit(randomh_wunifrac.ord~randomh.css.df$Number.in.Pool)
randomh_wunifrac_centroids

randomh_wunifrac_NMDS_col1 <- c( 0.0011,-0.0095,-0.0048,-0.0016)
randomh_wunifrac_NMDS_col2 <- c(0.0016,-0.0196,-0.0019,-0.0017)
randomh_wunifrac_centroids.df <-data.frame(pool_col,randomh_wunifrac_NMDS_col1,randomh_wunifrac_NMDS_col2)

##lr
randoml.css.df$Number.in.Pool <- as.factor(randoml.css.df$Number.in.Pool)
randoml_wunifrac_plot1 <- ordiplot(randoml_wunifrac.ord$points)
randoml_wunifrac_siteslong <- sites.long(randoml_wunifrac_plot1,randoml.css.df)

randoml_wunifrac_centroids <- envfit(randoml_wunifrac.ord~randoml.css.df$Number.in.Pool)
randoml_wunifrac_centroids

randoml_wunifrac_NMDS_col1 <- c(0.0027,-0.0152,0.0056,-0.0224)
randoml_wunifrac_NMDS_col2 <- c(0.0036,0.0071,0.0089,-0.0583)
randoml_wunifrac_centroids.df <-data.frame(pool_col,randoml_wunifrac_NMDS_col1,randoml_wunifrac_NMDS_col2)

##hr
pools_3v6h_wunifrac_plot1 <- ordiplot(pools_3v6h_wunifrac.ord$points)
pools_3v6h_wunifrac_siteslong <- sites.long(pools_3v6h_wunifrac_plot1,pools_3v6h.css.df)

pools_3v6h_wunifrac_centroids <- envfit(pools_3v6h_wunifrac.ord~pools_3v6h.css.df$PoolingType)
pools_3v6h_wunifrac_centroids

pooltypepalette <- distinctColorPalette(k=2) #"#B5D1A1" "#C17FCA"
pooltypepalette

type_col <- c("R","S")

pools_3v6h_wunifrac_NMDS_col1 <- c(0.0088, -0.0072)
pools_3v6h_wunifrac_NMDS_col2 <- c(0.0115, -0.0094)
pools_3v6h_wunifrac_centroids.df <-data.frame(type_col,pools_3v6h_wunifrac_NMDS_col1,pools_3v6h_wunifrac_NMDS_col2)

##lr
pools_3v6l.css.df$PoolingType <- as.factor(pools_3v6l.css.df$PoolingType)
pools_3v6l_wunifrac_plot1 <- ordiplot(pools_3v6l_wunifrac.ord$points)
pools_3v6l_wunifrac_siteslong <- sites.long(pools_3v6l_wunifrac_plot1,pools_3v6l.css.df)

pools_3v6l_wunifrac_centroids <- envfit(pools_3v6l_wunifrac.ord~pools_3v6l.css.df$PoolingType)
pools_3v6l_wunifrac_centroids

pools_3v6l_wunifrac_NMDS_col1 <- c(-0.0136,0.0136)
pools_3v6l_wunifrac_NMDS_col2 <- c(-0.0072,0.0072)
pools_3v6l_wunifrac_centroids.df <-data.frame(type_col,pools_3v6l_wunifrac_NMDS_col1,pools_3v6l_wunifrac_NMDS_col2)


##plot+stats
###high
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Weighted Unifrac") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = high_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,4,4,4)) +
  stat_ellipse(data = high_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = high_wunifrac_centroids_df, aes(x=high_wunifrac_NMDS_col1, y=high_wunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(23,19,19,19)) +
  geom_text(data = high_wunifrac_centroids_df, aes(x=high_wunifrac_NMDS_col1, y=high_wunifrac_NMDS_col2, label = pool_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, colour = "black"),
        plot.title=element_blank())

adonis2(high_wunifrac.dist ~ Number.in.Pool, high.css.df, p.adjust.methods = "BH", nperm = 9999) # NS, P=0.535

high_wunifrac.disper <- betadisper(high_wunifrac.dist, high.css.df$Number.in.Pool)
plot(high_wunifrac.disper)
boxplot(high_wunifrac.disper)
TukeyHSD(high_wunifrac.disper)
anova(high_wunifrac.disper)
high_wunifrac.permdisp <- permutest(high_wunifrac.disper, permutations = 9999, pairwise = T)
high_wunifrac.permdisp #NS


###hs
hs_wunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Weighted Unifrac") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subseth_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2)) +
  stat_ellipse(data = subseth_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subseth_wunifrac_centroids.df, aes(x=subseth_wunifrac_NMDS_col1, y=subseth_wunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subseth_wunifrac_centroids.df, aes(x=subseth_wunifrac_NMDS_col1, y=subseth_wunifrac_NMDS_col2, label = subseth_wunifrac_col), colour = "white", size = 2, fontface = "bold") +
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

hs_wunifrac_PCoA_plot <- ggplot() + theme_bw() +
  #labs(x= "NMDS1", y= "NMDS2", title = "Weighted Unifrac") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subseth_wunifrac_pcoasites, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(2,2,2,2)) +
  stat_ellipse(data = subseth_wunifrac_pcoasites, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subseth_wunifrac_PCoA_centroids.df, aes(x=subseth_wunifrac_PCoA_col1, y=subseth_wunifrac_PCoA_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subseth_wunifrac_PCoA_centroids.df, aes(x=subseth_wunifrac_PCoA_col1, y=subseth_wunifrac_PCoA_col2, label = pool_col), colour = "white", size = c(3,3,3,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"),
        #axis.title.x = element_blank(),
        plot.title=element_blank())

plot_ordination(subseth.css, subseth_wunifrac.pcoa,type = "samples", color = "Number.in.Pool", shape = "Number.in.Pool") +
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

adonis2(subseth_wunifrac.dist ~ Number.in.Pool, subseth.css.df, nperm = 9999)

hs_wunifrac_pa <- pairwise.adonis2(subseth_wunifrac.dist ~ Number.in.Pool, subseth.css.df, nperm = 9999, p.adjust.methods="BH") # NS

write.csv(hs_wunifrac_pa,"hs_wunifrac_permanova.csv")


subseth_wunifrac.disper <- betadisper(subseth_wunifrac.dist, subseth.css.df$Number.in.Pool)
plot(subseth_wunifrac.disper)
boxplot(subseth_wunifrac.disper)
TukeyHSD(subseth_wunifrac.disper)
anova(subseth_wunifrac.disper)
permutest(subseth_wunifrac.disper, permutations = 9999, pairwise=F)
subseth_wunifrac.permdisp <- permutest(subseth_wunifrac.disper, permutations = 9999, pairwise = T)
subseth_wunifrac.permdisp #NS
write.csv(subseth_wunifrac.permdisp[["pairwise"]][["permuted"]],"hs_wunifrac_permdisp.csv")

#PCoA
plot_ordination(subseth.css, subseth_wunifrac.pcoa, type = "samples", color = "Number.in.Pool") +
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
subsetl_wunifrac_siteslong$Number.in.Pool <- as.factor(subsetl_wunifrac_siteslong$Number.in.Pool)
ls_wunifrac_nmdsplot <-ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subsetl_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = subsetl_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subsetl_wunifrac_centroids.df, aes(x=subsetl_wunifrac_NMDS_col1, y=subsetl_wunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = subsetl_wunifrac_centroids.df, aes(x=subsetl_wunifrac_NMDS_col1, y=subsetl_wunifrac_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
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

ggsave("ls_wunifrac.tiff", plot = ls_wunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(subsetl_wunifrac.dist ~ Number.in.Pool, subsetl.css.df, nperm = 9999) # NS
ls_wunifrac_pa <- pairwise.adonis2(subsetl_wunifrac.dist ~ Number.in.Pool, subsetl.css.df, p.adjust.methods = "BH", nperm = 9999)
write.csv(ls_wunifrac_pa,"ls_wunifrac_permanova.csv")


subsetl_wunifrac.disper <- betadisper(subsetl_wunifrac.dist, subsetl.css.df$Number.in.Pool)
plot(subsetl_wunifrac.disper)
boxplot(subsetl_wunifrac.disper)
TukeyHSD(subsetl_wunifrac.disper)
anova(subsetl_wunifrac.disper)
permutest(subsetl_wunifrac.disper, permutations = 9999, pairwise = F)
subsetl_wunifrac.permdisp <- permutest(subsetl_wunifrac.disper, permutations = 9999, pairwise = T)
subsetl_wunifrac.permdisp #no significant difference in homogeneity of dispersion (overall 0.22, pairwise 6 and 1 different)

###hr
randomh_wunifrac_siteslong$Number.in.Pool <- as.factor(randomh_wunifrac_siteslong$Number.in.Pool)
hr_wunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randomh_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randomh_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randomh_wunifrac_centroids.df, aes(x=randomh_wunifrac_NMDS_col1, y=randomh_wunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = randomh_wunifrac_centroids.df, aes(x=randomh_wunifrac_NMDS_col1, y=randomh_wunifrac_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.title.x= element_blank(),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("hr_wunifrac.tiff", plot = hr_wunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(randomh_wunifrac.dist ~ Number.in.Pool, randomh.css.df, p.adjust.methods = "BH", nperm = 9999) # NS

randomh_wunifrac.disper <- betadisper(randomh_wunifrac.dist, randomh.css.df$Number.in.Pool)
plot(randomh_wunifrac.disper)
boxplot(randomh_wunifrac.disper)
TukeyHSD(randomh_wunifrac.disper)
anova(randomh_wunifrac.disper)
permutest(randomh_wunifrac.disper, permutations = 9999, pairwise = F)
randomh_wunifrac.permdisp <- permutest(randomh_wunifrac.disper, permutations = 9999, pairwise = T)
randomh_wunifrac.permdisp #NS, 0.07 overall, 1v12 p=0.052

###lr
randoml_wunifrac_siteslong$Number.in.Pool <- as.factor(randoml_wunifrac_siteslong$Number.in.Pool)
lr_wunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randoml_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randoml_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randoml_wunifrac_centroids.df, aes(x=randoml_wunifrac_NMDS_col1, y=randoml_wunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3, shape = c(18,19,19,19)) +
  geom_text(data = randoml_wunifrac_centroids.df, aes(x=randoml_wunifrac_NMDS_col1, y=randoml_wunifrac_NMDS_col2, label = pool_col), colour = "white", size = 2, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("lr_wunifrac.tiff", plot = lr_wunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm", dpi=600)


adonis2(randoml_wunifrac.dist ~ Number.in.Pool, randoml.css.df, nperm = 9999) # NS, p=0.072

randoml_wunifrac.disper <-betadisper(randoml_wunifrac.dist, randoml.css.df$Number.in.Pool)
plot(randoml_wunifrac.disper)
boxplot(randoml_wunifrac.disper)
TukeyHSD(randoml_wunifrac.disper)
anova(randoml_wunifrac.disper)
permutest(randoml_wunifrac.disper, permutations = 9999, pairwise = F)
randoml_wunifrac.permdisp <- permutest(randoml_wunifrac.disper, permutations = 9999, pairwise = T)
randoml_wunifrac.permdisp #NS overall p=0.22


###pool36h
pools_3v6h_wunifrac_siteslong$PoolingType <- as.factor(pools_3v6h_wunifrac_siteslong$PoolingType)
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19)) +
  geom_point(data = pools_3v6h_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= PoolingType, shape = PoolingType, alpha = PoolingType, size = PoolingType)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_3v6h_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= PoolingType, fill = PoolingType), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = pools_3v6h_wunifrac_centroids.df, aes(x=pools_3v6h_wunifrac_NMDS_col1, y=pools_3v6h_wunifrac_NMDS_col2), fill = pooltypepalette, colour = pooltypepalette, size = 10, shape = c(19,19)) +
  geom_text(data = pools_3v6h_wunifrac_centroids.df, aes(x=pools_3v6h_wunifrac_NMDS_col1, y=pools_3v6h_wunifrac_NMDS_col2, label = type_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values =pooltypepalette) +
  scale_fill_manual(values = pooltypepalette) +
  facet_wrap(~Number.in.Pool)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 14, colour = "black"))

adonis2(pools_3v6h_wunifrac.dist ~ PoolingType, pools_3v6h.css.df,  nperm = 9999) # NS

pools_3v6h_wunifrac.disper <-betadisper(pools_3v6h_wunifrac.dist, pools_3v6h.css.df$PoolingType)
plot(pools_3v6h_wunifrac.disper)
pools_3v6h_wunifrac.permdisp <- permutest(pools_3v6h_wunifrac.disper, permutations = 9999)
pools_3v6h_wunifrac.permdisp #NS

###pool36h
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19,19,18)) +
  geom_point(data = pools_3v6l_wunifrac_siteslong, aes(x=axis1,y=axis2, colour= PoolingType, shape = PoolingType, alpha = PoolingType, size = PoolingType)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_3v6l_wunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= PoolingType, fill = PoolingType), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = pools_3v6l_wunifrac_centroids.df, aes(x=pools_3v6l_wunifrac_NMDS_col1, y=pools_3v6l_wunifrac_NMDS_col2), fill = pooltypepalette, colour = pooltypepalette, size = 14, shape = c(19,19)) +
  geom_text(data = pools_3v6l_wunifrac_centroids.df, aes(x=pools_3v6l_wunifrac_NMDS_col1, y=pools_3v6l_wunifrac_NMDS_col2, label = type_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values =pooltypepalette) +
  scale_fill_manual(values = pooltypepalette) +
  #facet_wrap(~Number.in.Pool)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

adonis2(pools_3v6l_wunifrac.dist ~ PoolingType, pools_3v6l.css.df, nperm = 9999) # NS

pools_3v6l_wunifrac.disper <-betadisper(pools_3v6l_wunifrac.dist, pools_3v6l.css.df$PoolingType)
plot(pools_3v6l_wunifrac.disper)
pools_3v6l_wunifrac.permdisp <- permutest(pools_3v6l_wunifrac.disper, permutations = 9999)
pools_3v6l_wunifrac.permdisp #NS

# plot unweighted unifrac
#### findingcenters
#hs
subseth_uwunifrac_plot1 <- ordiplot(subseth_uwunifrac.ord$points)
subseth_uwunifrac_siteslong <-sites.long(subseth_uwunifrac_plot1,subseth.css.df)

subseth_uwunifrac_centroids <- envfit(subseth_uwunifrac.ord~subseth.css.df$Number.in.Pool)
subseth_uwunifrac_centroids


subseth_uwunifrac_NMDS_col1 <- c(0.0299,-0.1772,-0.0207,-0.1589)
subseth_uwunifrac_NMDS_col2 <- c(0.0145,-0.0631,-0.0651,-0.0542)
subseth_uwunifrac_centroids.df <-data.frame(pool_col,subseth_uwunifrac_NMDS_col1,subseth_uwunifrac_NMDS_col2)

##ls
subsetl_uwunifrac_plot1 <- ordiplot(subsetl_uwunifrac.ord$points)
subsetl_uwunifrac_siteslong <- sites.long(subsetl_uwunifrac_plot1,subsetl.css.df)
subsetl_uwunifrac_centroids <- envfit(subsetl_uwunifrac.ord~subsetl.css.df$Number.in.Pool)
subsetl_uwunifrac_centroids

subsetl_uwunifrac_NMDS_col1 <- c(0.0361,0.0027,-0.1096,-0.3141)
subsetl_uwunifrac_NMDS_col2 <- c(-0.0239,0.1477,0.0811,0.0505)
subsetl_uwunifrac_centroids.df <-data.frame(pool_col,subsetl_uwunifrac_NMDS_col1,subsetl_uwunifrac_NMDS_col2)

##hr
randomh_uwunifrac_plot1 <- ordiplot(randomh_uwunifrac.ord$points)
randomh_uwunifrac_siteslong <- sites.long(randomh_uwunifrac_plot1,randomh.css.df)

randomh_uwunifrac_centroids <- envfit(randomh_uwunifrac.ord~randomh.css.df$Number.in.Pool)


randomh_uwunifrac_NMDS_col1 <- c(0.0192,0.0070,-0.0741,-0.1838)
randomh_uwunifrac_NMDS_col2 <- c(0.0052,0.0173,-0.0234,-0.0606)
randomh_uwunifrac_centroids.df <-data.frame(pool_col,randomh_uwunifrac_NMDS_col1,randomh_uwunifrac_NMDS_col2)

##lr
randoml_uwunifrac_plot1 <- ordiplot(randoml_uwunifrac.ord$points)
randoml_uwunifrac_siteslong <- sites.long(randoml_uwunifrac_plot1,randoml.css.df)

randoml_uwunifrac_centroids <- envfit(randoml_uwunifrac.ord~randoml.css.df$Number.in.Pool)
randoml_uwunifrac_centroids

randoml_uwunifrac_NMDS_col1 <- c(0.0527,-0.1165,-0.2004,-0.2979)
randoml_uwunifrac_NMDS_col2 <- c(-0.0132,0.0628,0.0318,0.0593)
randoml_uwunifrac_centroids.df <-data.frame(pool_col,randoml_uwunifrac_NMDS_col1,randoml_uwunifrac_NMDS_col2)

##poolh
pools_3v6h_uwunifrac_plot1 <- ordiplot(pools_3v6h_uwunifrac.ord$points)
pools_3v6h_uwunifrac_siteslong <- sites.long(pools_3v6h_uwunifrac_plot1,pools_3v6h.css.df)

pools_3v6h_uwunifrac_centroids <- envfit(pools_3v6h_uwunifrac.ord~pools_3v6h.css.df$PoolingType)
pools_3v6h_uwunifrac_centroids

pools_3v6h_uwunifrac_NMDS_col1 <- c(0.0479, -0.0392)
pools_3v6h_uwunifrac_NMDS_col2 <- c(0.0008,-0.0007)
pools_3v6h_uwunifrac_centroids.df <-data.frame(type_col,pools_3v6h_uwunifrac_NMDS_col1,pools_3v6h_uwunifrac_NMDS_col2)

##pooll
pools_3v6l_uwunifrac_plot1 <- ordiplot(pools_3v6l_uwunifrac.ord$points)
pools_3v6l_uwunifrac_siteslong <- sites.long(pools_3v6l_uwunifrac_plot1,pools_3v6l.css.df)

pools_3v6l_uwunifrac_centroids <- envfit(pools_3v6l_uwunifrac.ord~pools_3v6l.css.df$PoolingType)
pools_3v6l_uwunifrac_centroids

pools_3v6l_uwunifrac_NMDS_col1 <- c(0.1629,-0.1629)
pools_3v6l_uwunifrac_NMDS_col2 <- c(0.0385,-0.0385)
pools_3v6l_uwunifrac_centroids.df <-data.frame(type_col,pools_3v6l_uwunifrac_NMDS_col1,pools_3v6l_uwunifrac_NMDS_col2)


##plot+stats
###hs
hs_uwunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous() +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subseth_uwunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = subseth_uwunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subseth_uwunifrac_centroids.df, aes(x=subseth_uwunifrac_NMDS_col1, y=subseth_uwunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3 , shape = c(18,19,19,19)) +
  geom_text(data = subseth_uwunifrac_centroids.df, aes(x=subseth_uwunifrac_NMDS_col1, y=subseth_uwunifrac_NMDS_col2, label = c(1,3,6,12)),colour = "white", size = c(2,2,2,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6, colour = "black"))
ggsave("hs_uwunifrac.tiff", plot = hs_uwunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm", dpi=600)

adonis2(subseth_uwunifrac.dist ~ Number.in.Pool, subseth.css.df,nperm = 9999) # different, 0.001
pairwise.adonis2(subseth_uwunifrac.dist ~ Number.in.Pool, subseth.css.df, p.adjust.methods = "BH", nperm = 9999)
#1 different from all (3=0.003, 6=0.007, 12=0.006), no difference among others


subseth_uwunifrac.disper <- betadisper(subseth_uwunifrac.dist, subseth.css.df$Number.in.Pool)
plot(subseth_uwunifrac.disper)
boxplot(subseth_uwunifrac.disper)
permutest(subseth_uwunifrac.disper, permutations = 9999, pairwise = F)
subseth_uwunifrac.permdisp <- permutest(subseth_uwunifrac.disper, permutations = 9999, pairwise = T, p.adjust.methods="BH")
subseth_uwunifrac.permdisp #significant difference in homogeneity of dispersion (0.0029), 1 different from 3 (0.0012) and 12 (0.0024), not different from 6 (0.076)

###ls
ls_uwunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subsetl_uwunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = subsetl_uwunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = subsetl_uwunifrac_centroids.df, aes(x=subsetl_uwunifrac_NMDS_col1, y=subsetl_uwunifrac_NMDS_col2),  fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3 , shape = c(18,19,19,19)) +
  geom_text(data = subsetl_uwunifrac_centroids.df, aes(x=subsetl_uwunifrac_NMDS_col1, y=subsetl_uwunifrac_NMDS_col2, label = pool_col), label = c(1,3,6,12),colour = "white", size = c(2,2,2,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("ls_uwunifrac.tiff", plot = ls_uwunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm")


adonis2(subsetl_uwunifrac.dist ~ Number.in.Pool, subsetl.css.df, nperm = 9999) # Significant, p=0.001
pairwise.adonis2(subsetl_uwunifrac.dist ~ Number.in.Pool, subsetl.css.df, p.adjust.methods = "BH", nperm = 9999) #1 different from all(3=0.002, 6=0.003, 12=0.001), no difference among others (12v3 close)

subsetl_uwunifrac.disper <- betadisper(subsetl_uwunifrac.dist, subsetl.css.df$Number.in.Pool)
plot(subsetl_uwunifrac.disper)
permutest(subsetl_uwunifrac.disper, permutations = 9999, pairwise = F)
subsetl_uwunifrac.permdisp <- permutest(subsetl_uwunifrac.disper, permutations = 9999, pairwise = T)
subsetl_uwunifrac.permdisp #significant difference in homogeneity of dispersion, (1 different from all, 12 not different from 3 and 6)

###hr
hr_uwunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(breaks = c(-0.3, 0, 0.3, 0.5)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randomh_uwunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randomh_uwunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool),  alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randomh_uwunifrac_centroids.df, aes(x=randomh_uwunifrac_NMDS_col1, y=randomh_uwunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3 , shape = c(18,19,19,19)) +
  geom_text(data = randomh_uwunifrac_centroids.df, aes(x=randomh_uwunifrac_NMDS_col1, y=randomh_uwunifrac_NMDS_col2, label = pool_col), colour = "white", size = c(2,2,2,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6, colour = "black"))
ggsave("hr_uwunifrac.tiff", plot = hr_uwunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm",dpi = 600)

adonis2(randomh_uwunifrac.dist ~ Number.in.Pool, randomh.css.df,  nperm = 9999) #P=0.001
pairwise.adonis2(randomh_uwunifrac.dist ~ Number.in.Pool, randomh.css.df, p.adjust.methods = "BH", nperm = 9999) # Different, 1 different from 6 and 12 (6=0.006, 12=0.008)

randomh_uwunifrac.disper <- betadisper(randomh_uwunifrac.dist, randomh.css.df$Number.in.Pool)
plot(randomh_uwunifrac.disper)
permutest(randomh_uwunifrac.disper, permutations = 9999, pairwise = F)
randomh_uwunifrac.permdisp <- permutest(randomh_uwunifrac.disper, permutations = 9999, pairwise = T)
randomh_uwunifrac.permdisp #1 different from 3 (0.0025), 6 (0.021) and 12 (0.0029)

###lr
lr_uwunifrac_nmdsplot <- ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous() +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randoml_uwunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(3,2,2,2), guide = "none") +
  stat_ellipse(data = randoml_uwunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.1), lty=1, level = 0.95, linewidth = 1) +
  geom_point(data = randoml_uwunifrac_centroids.df, aes(x=randoml_uwunifrac_NMDS_col1, y=randoml_uwunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 3 , shape = c(18,19,19,19)) +
  geom_text(data = randoml_uwunifrac_centroids.df, aes(x=randoml_uwunifrac_NMDS_col1, y=randoml_uwunifrac_NMDS_col2, label = pool_col), label = c(1,3,6,12),colour = "white", size = c(2,2,2,2), fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = "black"))

ggsave("lr_uwunifrac.tiff", plot = lr_uwunifrac_nmdsplot, path="~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Publication/Figures", width = 90, height = 50, units = "mm", dpi = 600)

adonis2(randoml_uwunifrac.dist~Number.in.Pool, randoml.css.df, nperm=9999)
pairwise.adonis2(randoml_uwunifrac.dist ~ Number.in.Pool, randoml.css.df, p.adjust.methods = "BH", nperm = 9999) # Different, 1 v all (3=0.005, 6=0.001, 12=0.002), 12 v 3 (0.038)

randoml_uwunifrac.disper <-betadisper(randoml_uwunifrac.dist, randoml.css.df$Number.in.Pool)
plot(randoml_uwunifrac.disper)
permutest(randoml_uwunifrac.disper, permutations = 9999, pairwise = F)
randoml_uwunifrac.permdisp <- permutest(randoml_uwunifrac.disper, permutations = 9999, pairwise = T)
randoml_uwunifrac.permdisp #1 different all (3 and 6=10-4, 12=0.0179)



###pool36h
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19)) +
  geom_point(data = pools_3v6h_uwunifrac_siteslong, aes(x=axis1,y=axis2, colour= PoolingType, shape = PoolingType, alpha = PoolingType, size = PoolingType)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_3v6h_uwunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= PoolingType, fill = PoolingType), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = pools_3v6h_uwunifrac_centroids.df, aes(x=pools_3v6h_uwunifrac_NMDS_col1, y=pools_3v6h_uwunifrac_NMDS_col2), fill = pooltypepalette, colour = pooltypepalette, size = 14, shape = c(19,19)) +
  geom_text(data = pools_3v6h_uwunifrac_centroids.df, aes(x=pools_3v6h_uwunifrac_NMDS_col1, y=pools_3v6h_uwunifrac_NMDS_col2, label = type_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values =pooltypepalette) +
  scale_fill_manual(values = pooltypepalette) +
  #facet_wrap(~Number.in.Pool)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

adonis2(pools_3v6h_uwunifrac.dist ~ PoolingType, pools_3v6h.css.df,  nperm = 9999) # NS, 0.87

pools_3v6h_uwunifrac.disper <-betadisper(pools_3v6h_uwunifrac.dist, pools_3v6h.css.df$PoolingType)
plot(pools_3v6h_uwunifrac.disper)
pools_3v6h_uwunifrac.permdisp <- permutest(pools_3v6h_uwunifrac.disper, permutations = 9999)
pools_3v6h_uwunifrac.permdisp #NS, p=0.84

###pool36h
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19)) +
  geom_point(data = pools_3v6l_uwunifrac_siteslong, aes(x=axis1,y=axis2, colour= PoolingType, shape = PoolingType, alpha = PoolingType, size = PoolingType)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_3v6l_uwunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= PoolingType, fill = PoolingType), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = pools_3v6l_uwunifrac_centroids.df, aes(x=pools_3v6l_uwunifrac_NMDS_col1, y=pools_3v6l_uwunifrac_NMDS_col2), fill = pooltypepalette, colour = pooltypepalette, size = 14, shape = c(19,19)) +
  geom_text(data = pools_3v6l_uwunifrac_centroids.df, aes(x=pools_3v6l_uwunifrac_NMDS_col1, y=pools_3v6l_uwunifrac_NMDS_col2, label = type_col), colour = "white", size = 9, fontface = "bold") +
  scale_colour_manual(values =pooltypepalette) +
  scale_fill_manual(values = pooltypepalette) +
  #facet_wrap(~Number.in.Pool)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

adonis2(pools_3v6l_uwunifrac.dist ~ PoolingType, pools_3v6l.css.df, nperm = 9999) #difference in low prevalence, p=0.001

pools_3v6l_uwunifrac.disper <-betadisper(pools_3v6l_uwunifrac.dist, pools_3v6l.css.df$PoolingType)
plot(pools_3v6l_uwunifrac.disper)
pools_3v6l_uwunifrac.permdisp <- permutest(pools_3v6l_uwunifrac.disper, permutations = 9999)
pools_3v6l_uwunifrac.permdisp #NS


# plot generalized unifrac
#### findingcenters
#hs
subseth_gunifrac_plot1 <- ordiplot(subseth_gunifrac.ord$points)
subseth_gunifrac_siteslong <-sites.long(subseth_gunifrac_plot1,subseth.css.df)

subseth_gunifrac_centroids <- envfit(subseth_gunifrac.ord~subseth.css.df$Number.in.Pool)
subseth_gunifrac_centroids


subseth_gunifrac_NMDS_col1 <- c(0.0115,-0.0840,0.0114,-0.0616)
subseth_gunifrac_NMDS_col2 <- c(-0.0094,0.0425,0.0439,0.0316)
subseth_gunifrac_centroids.df <-data.frame(pool_col,subseth_gunifrac_NMDS_col1,subseth_gunifrac_NMDS_col2)

##ls
subsetl_gunifrac_plot1 <- ordiplot(subsetl_gunifrac.ord$points)
subsetl_gunifrac_siteslong <- sites.long(subsetl_gunifrac_plot1,subsetl.css.df)
subsetl_gunifrac_centroids <- envfit(subsetl_gunifrac.ord~subsetl.css.df$Number.in.Pool)
subsetl_gunifrac_centroids

subsetl_gunifrac_NMDS_col1 <- c(0.0202,0.0635,-0.0049,-0.2949)
subsetl_gunifrac_NMDS_col2 <- c(-0.0031,0.0608,0.0200,-0.0451)
subsetl_gunifrac_centroids.df <-data.frame(pool_col,subsetl_gunifrac_NMDS_col1,subsetl_gunifrac_NMDS_col2)

##hr
randomh_gunifrac_plot1 <- ordiplot(randomh_gunifrac.ord$points)
randomh_gunifrac_siteslong <- sites.long(randomh_gunifrac_plot1,randomh.css.df)

randomh_gunifrac_centroids <- envfit(randomh_gunifrac.ord~randomh.css.df$Number.in.Pool)
randomh_gunifrac_centroids

randomh_gunifrac_NMDS_col1 <- c(-0.0021,0.0691,0.0541,-0.0685)
randomh_gunifrac_NMDS_col2 <- c(-0.0010,0.0252,-0.0605,0.0314)
randomh_gunifrac_centroids.df <-data.frame(pool_col,randomh_gunifrac_NMDS_col1,randomh_gunifrac_NMDS_col2)

##lr
randoml_gunifrac_plot1 <- ordiplot(randoml_gunifrac.ord$points)
randoml_gunifrac_siteslong <- sites.long(randoml_gunifrac_plot1,randoml.css.df)

randoml_gunifrac_centroids <- envfit(randoml_gunifrac.ord~randoml.css.df$Number.in.Pool)
randoml_gunifrac_centroids

randoml_gunifrac_NMDS_col1 <- c(0.0292,0.0175,-0.0717,-0.2870)
randoml_gunifrac_NMDS_col2 <- c(0.0074,-0.0501,0.0032,-0.0397)
randoml_gunifrac_centroids.df <-data.frame(pool_col,randoml_gunifrac_NMDS_col1,randoml_gunifrac_NMDS_col2)

##poolh
pools_3v6h_gunifrac_plot1 <- ordiplot(pools_3v6h_gunifrac.ord$points)
pools_3v6h_gunifrac_siteslong <- sites.long(pools_3v6h_gunifrac_plot1,pools_3v6h.css.df)

pools_3v6h_gunifrac_centroids <- envfit(pools_3v6h_gunifrac.ord~pools_3v6h.css.df$PoolingType)
pools_3v6h_gunifrac_centroids

pools_3v6h_gunifrac_NMDS_col1 <- c(0.0940, -0.0769)
pools_3v6h_gunifrac_NMDS_col2 <- c(-0.0214,0.0175)
pools_3v6h_gunifrac_centroids.df <-data.frame(type_col,pools_3v6h_gunifrac_NMDS_col1,pools_3v6h_gunifrac_NMDS_col2)

##pooll
pools_3v6l_gunifrac_plot1 <- ordiplot(pools_3v6l_gunifrac.ord$points)
pools_3v6l_gunifrac_siteslong <- sites.long(pools_3v6l_gunifrac_plot1,pools_3v6l.css.df)

pools_3v6l_gunifrac_centroids <- envfit(pools_3v6l_gunifrac.ord~pools_3v6l.css.df$PoolingType)
pools_3v6l_gunifrac_centroids

pools_3v6l_gunifrac_NMDS_col1 <- c(0.0915,-0.0915)
pools_3v6l_gunifrac_NMDS_col2 <- c(0.0206,-0.0206)
pools_3v6l_gunifrac_centroids.df <-data.frame(type_col,pools_3v6l_gunifrac_NMDS_col1,pools_3v6l_gunifrac_NMDS_col2)


##plot+stats
###hs
subseth_gunifrac_siteslong$Number.in.Pool <- as.factor(subseth_gunifrac_siteslong$Number.in.Pool)
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2", title = "Beta Diversity, High Prevalence Subset Pooling (Weighted Unifrac)") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subseth_gunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5,5,5), guide = "none") +
  stat_ellipse(data = subseth_gunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = subseth_gunifrac_centroids.df, aes(x=subseth_gunifrac_NMDS_col1, y=subseth_gunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(23,19,19,19)) +
  geom_text(data = subseth_gunifrac_centroids.df, aes(x=subseth_gunifrac_NMDS_col1, y=subseth_gunifrac_NMDS_col2, label = pool_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 8, colour = "black"),
    plot.title=element_blank())

adonis2(subseth_gunifrac.dist ~ Number.in.Pool, subseth.css.df, p.adjust.methods = "BH", nperm = 9999) # different, 0.007
pairwise.adonis2(subseth_gunifrac.dist ~ Number.in.Pool, subseth.css.df, p.adjust.methods = "BH", nperm = 9999)
#1 different from 12 (0.013) and 3 (0.022)


subseth_gunifrac.disper <- betadisper(subseth_gunifrac.dist, subseth.css.df$Number.in.Pool)
plot(subseth_gunifrac.disper)
subseth_gunifrac.permdisp <- permutest(subseth_gunifrac.disper, permutations = 9999, pairwise = T)
subseth_gunifrac.permdisp #NS

###ls
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = subsetl_gunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(6,5,5,5), guide = "none") +
  stat_ellipse(data = subsetl_gunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = subsetl_gunifrac_centroids.df, aes(x=subsetl_gunifrac_NMDS_col1, y=subsetl_gunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(23,19,19,19)) +
  geom_text(data = subsetl_gunifrac_centroids.df, aes(x=subsetl_gunifrac_NMDS_col1, y=subsetl_gunifrac_NMDS_col2, label = pool_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

adonis2(subsetl_gunifrac.dist ~ Number.in.Pool, subsetl.css.df, nperm = 9999) # Signficant, 0.001
pairwise.adonis2(subsetl_gunifrac.dist ~ Number.in.Pool, subsetl.css.df, p.adjust.methods = "BH", nperm = 9999) #1 different from 6(p=0.042) and 12 (p=0.002), no difference among others (12v3 close)

subsetl_gunifrac.disper <- betadisper(subsetl_gunifrac.dist, subsetl.css.df$Number.in.Pool)
plot(subsetl_gunifrac.disper)
subsetl_gunifrac.permdisp <- permutest(subsetl_gunifrac.disper, permutations = 9999, pairwise = T)
subsetl_gunifrac.permdisp #significant difference in homogeneity of dispersion, 12 different from 3 (p=0.017) and 6 (p=0.034)

###hr
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randomh_gunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(6,5,5,5), guide = "none") +
  stat_ellipse(data = randomh_gunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = randomh_gunifrac_centroids.df, aes(x=randomh_gunifrac_NMDS_col1, y=randomh_gunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(19,19,19,23)) +
  geom_text(data = randomh_gunifrac_centroids.df, aes(x=randomh_gunifrac_NMDS_col1, y=randomh_gunifrac_NMDS_col2, label = pool_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

pairwise.adonis2(randomh_gunifrac.dist ~ Number.in.Pool, randomh.css.df, p.adjust.methods = "BH", nperm = 9999) # Different (0.02), 1 different from  and 12 (0.014)

randomh_gunifrac.disper <- betadisper(randomh_gunifrac.dist, randomh.css.df$Number.in.Pool)
plot(randomh_gunifrac.disper)
randomh_gunifrac.permdisp <- permutest(randomh_gunifrac.disper, permutations = 9999, pairwise = T)
randomh_gunifrac.permdisp #NS

###lr
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = randoml_gunifrac_siteslong, aes(x=axis1,y=axis2, colour= Number.in.Pool, shape = Number.in.Pool, alpha = Number.in.Pool, size = Number.in.Pool)) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.5), guide = F) +
  scale_size_manual(values = c(6,5,5,5), guide = "none") +
  stat_ellipse(data = randoml_gunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= Number.in.Pool, fill = Number.in.Pool), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = randoml_gunifrac_centroids.df, aes(x=randoml_gunifrac_NMDS_col1, y=randoml_gunifrac_NMDS_col2), fill = numberinpoolpalette, colour = numberinpoolpalette, size = 10, shape = c(23,19,19,19)) +
  geom_text(data = randoml_gunifrac_centroids.df, aes(x=randoml_gunifrac_NMDS_col1, y=randoml_gunifrac_NMDS_col2, label = pool_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = numberinpoolpalette) +
  scale_fill_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

pairwise.adonis2(randoml_gunifrac.dist ~ Number.in.Pool, randoml.css.df, p.adjust.methods = "BH", nperm = 9999) # Different, 1 v 6 (0.042) & 12 (0.002)

randoml_gunifrac.disper <-betadisper(randoml_gunifrac.dist, randoml.css.df$Number.in.Pool)
plot(randoml_gunifrac.disper)
randoml_gunifrac.permdisp <- permutest(randoml_gunifrac.disper, permutations = 9999, pairwise = T)
randoml_gunifrac.permdisp #12 different from 3 (0.0179) and 6 (0.0362)


###pool36h
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(18,19,19,19)) +
  geom_point(data = pools_3v6h_gunifrac_siteslong, aes(x=axis1,y=axis2, colour= PoolingType, shape = PoolingType, alpha = PoolingType, size = PoolingType)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_3v6h_gunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= PoolingType, fill = PoolingType), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = pools_3v6h_gunifrac_centroids.df, aes(x=pools_3v6h_gunifrac_NMDS_col1, y=pools_3v6h_gunifrac_NMDS_col2), fill = pooltypepalette, colour = pooltypepalette, size = 10, shape = c(19,19)) +
  geom_text(data = pools_3v6h_gunifrac_centroids.df, aes(x=pools_3v6h_gunifrac_NMDS_col1, y=pools_3v6h_gunifrac_NMDS_col2, label = type_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values =pooltypepalette) +
  scale_fill_manual(values = pooltypepalette) +
  #facet_wrap(~Number.in.Pool)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

adonis2(pools_3v6h_gunifrac.dist ~ PoolingType, pools_3v6h.css.df,  nperm = 9999) # NS (0.087)

pools_3v6h_gunifrac.disper <-betadisper(pools_3v6h_gunifrac.dist, pools_3v6h.css.df$PoolingType)
plot(pools_3v6h_gunifrac.disper)
pools_3v6h_gunifrac.permdisp <- permutest(pools_3v6h_gunifrac.disper, permutations = 9999)
pools_3v6h_gunifrac.permdisp #NS, 0.0632

###pool36h
ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_shape_manual(values = c(19,19)) +
  geom_point(data = pools_3v6l_gunifrac_siteslong, aes(x=axis1,y=axis2, colour= PoolingType, shape = PoolingType, alpha = PoolingType, size = PoolingType)) +
  scale_alpha_manual(values = c(0.5,0.5), guide = F) +
  scale_size_manual(values = c(5,5), guide = "none") +
  stat_ellipse(data = pools_3v6l_gunifrac_siteslong, geom = "polygon", aes(x=axis1,y=axis2, colour= PoolingType, fill = PoolingType), alpha = c(0.06), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = pools_3v6l_gunifrac_centroids.df, aes(x=pools_3v6l_gunifrac_NMDS_col1, y=pools_3v6l_gunifrac_NMDS_col2), fill = pooltypepalette, colour = pooltypepalette, size = 10, shape = c(19,19)) +
  geom_text(data = pools_3v6l_gunifrac_centroids.df, aes(x=pools_3v6l_gunifrac_NMDS_col1, y=pools_3v6l_gunifrac_NMDS_col2, label = type_col), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values =pooltypepalette) +
  scale_fill_manual(values = pooltypepalette) +
  #facet_wrap(~Number.in.Pool)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8, colour = "black"))

adonis2(pools_3v6l_gunifrac.dist ~ PoolingType, pools_3v6l.css.df, nperm = 9999) #difference in low prevalence (0.012)

pools_3v6l_gunifrac.disper <-betadisper(pools_3v6l_gunifrac.dist, pools_3v6l.css.df$PoolingType)
plot(pools_3v6l_gunifrac.disper)
pools_3v6l_gunifrac.permdisp <- permutest(pools_3v6l_gunifrac.disper, permutations = 9999)
pools_3v6l_gunifrac.permdisp #NS

######RELATIVE ABUNDANCE



