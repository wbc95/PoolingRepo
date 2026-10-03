##Upset Plot

library(UpSetR)
library(MicrobiotaProcess)
library(MicEco)
library(ggupset)

BiocManager::install("MicrobiotaProcess")

BiocManager::install("randomForest")

#Subset
upsetda_16S_sh <- get_upset(sh_16S.css.ra , factorNames="Number.in.Pool") ## ASV
sh_melt_asv <- psmelt(subseth.css.ra)

upsetda_16S_sl <- get_upset(sl_16S.css.ra, factorNames="Number.in.Pool")

#random
upsetda_16S_HR <- get_upset(HR_16S.css.ra , factorNames="Number.in.Pool") ## ASV
upsetda_16S_LR <- get_upset(LR_16S.css.ra, factorNames="Number.in.Pool")

##All
upset_data <- get_upset(data.css.ra, factorNames="Number.in.Pool")

#Split
upset_all <- upset(upset_data, sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.2)


upset_16S_sh <- upset(upsetda_16S_sh , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="no", mainbar.y.label="", text.scale=1.2)
upset_16S_sl <- upset(upsetda_16S_sl, sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.2)


upset_16S_rh <- upset(upsetda_16S_HR , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.2)
upset_16S_rl <- upset(upsetda_16S_LR, sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.2)

#Compute prevalence of each feature, store as data.frame for each sample type 
prevdf_16S_HS  <- as.data.frame(prevalence(sh_16S.css.ra, detection = 0))
prevdf_16S_LS  <- as.data.frame(prevalence(sl_16S.css.ra, detection = 0))
prevdf_16S_HR  <- as.data.frame(prevalence(HR_16S.css.ra, detection = 0))
prevdf_16S_LR  <- as.data.frame(prevalence(LR_16S.css.ra, detection = 0))



#Extract relative abundance
abundancedf_16S_HS <- as.data.frame(taxa_sums(sh_16S.css.ra))
write.csv(prevdf_16S_HS, "OTU_prevalence_16S_HS.csv") #Did this to check they were in same order
write.csv(abundancedf_16S_HS, "OTU_ra_16S_HS.csv") #Did this to check they were in the same order

abundancedf_16S_LS <- as.data.frame(taxa_sums(sl_16S.css.ra))
write.csv(prevdf_16S_LS, "OTU_prevalence_16S_LS.csv") #Did this to check they were in same order
write.csv(abundancedf_16S_LS, "OTU_ra_16S_LS.csv") #Did this to check they were in the same order

abundancedf_16S_HR <- as.data.frame(taxa_sums(HR_16S.css.ra))
write.csv(prevdf_16S_HR, "OTU_prevalence_16S_HR.csv") #Did this to check they were in same order
write.csv(abundancedf_16S_HR, "OTU_ra_16S_HR.csv") #Did this to check they were in the same order

abundancedf_16S_LR <- as.data.frame(taxa_sums(LR_16S.css.ra))
write.csv(prevdf_16S_LR, "OTU_prevalence_16S_LR.csv") #Did this to check they were in same order
write.csv(abundancedf_16S_LR, "OTU_ra_16S_LR.csv") #Did this to check they were in the same order



##Idea..subset taxa, with RA=0...repeat each pool...will have to take this to Excel and manually add using =vlookup or =if..maybe too many rows
otus_16S_HS <- as.data.frame(otu_table(sh_16S.css.ra))
colnames(otus_16S_HS)
otus_16S_HS<-subset(otus_16S_HS,`H_12-01`==0&`H_12-02`==0&`H_12-03`==0&`H_12-04`==0&`H_12-05`==0&`H_12-06`==0&`HP_03-01`==0&`HP_03-02`==0&`HP_03-03`==0&`HP_03-04`==0&`HP_03-05`==0&`HP_03-06`==0&`HP_06-01`==0&`HP_06-02`==0&`HP_06-03`==0&`HP_06-04`==0&`HP_06-06`==0)

sample_names(sh_16S.css.ra) <- sample_data(sh_16S.css.ra)$sample.id
otus_16S_LS <- as.data.frame(otu_table(sl_16S.css.ra))
colnames(otus_16S_LS)
otus_16S_LS <- subset(otus_16S_LS,`L_12-01`==0&`L_12-02`==0&`L_12-03`==0&`L_12-04`==0&`L_12-05`==0&`L_12-06`==0&`LP_03-01`==0&`LP_03-02`==0&`LP_03-03`==0&`LP_03-04`==0&`LP_03-05`==0&`LP_03-06`==0&`LP_06-01`==0&`LP_06-02`==0&`LP_06-03`==0&`LP_06-04`==0&`LP_06-05`==0&`LP_06-06`==0)

otus_16S_HR <- as.data.frame(otu_table(HR_16S.css.ra))
colnames(otus_16S_HR)
otus_16S_HR <- subset(otus_16S_HR,`H_12-01`==0&`H_12-02`==0&`H_12-03`==0&`H_12-04`==0&`H_12-05`==0&`H_12-06`==0&`HR_03-01`==0&`HR_03-02`==0&`HR_03-03`==0&`HR_03-04`==0&`HR_03-05`==0&`HR_06-01`==0&`HR_06-02`==0&`HR_06-03`==0&`HR_06-04`==0)

otus_16S_LR <- as.data.frame(otu_table(LR_16S.css.ra))
colnames(otus_16S_LR)
otus_16S_LR <- subset(otus_16S_LR,`L_12-01`==0&`L_12-02`==0&`L_12-03`==0&`L_12-04`==0&`L_12-05`==0&`L_12-06`==0&`LR_03-01`==0&`LR_03-02`==0&`LR_03-03`==0&`LR_03-04`==0&`LR_03-05`==0&`LR_03-06`==0&`LR_06-01`==0&`LR_06-02`==0&`LR_06-03`==0&`LR_06-04`==0&`LR_06-05`==0&`LR_06-06`==0)

##I think that worked
write.csv(otus_16S_HS,"OTU_inviduals_HS.csv")
write.csv(otus_16S_LS,"OTU_inviduals_LS.csv")
write.csv(otus_16S_HR,"OTU_inviduals_HR.csv")
write.csv(otus_16S_LR,"OTU_inviduals_LR.csv")

#combine df, used =countif and =if and copied values into column in csv with relative abundance
abundancewindividual_16S_HS <- read.csv("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/RProjects/OTU_ra_16S_HS_windi.csv")
prev_abund_16S_HS_df <- cbind(prevdf_16S_HS, abundancewindividual_16S_HS)
colnames(prev_abund_16S_HS_df) <- c("Proportion","ASV","RelativeAbundance", "IndividualOnly")
prev_abund_16S_HS_df <- transform.data.frame(prev_abund_16S_HS_df, RA=(RelativeAbundance/88))
prev_abund_16S_HS_df <- transform.data.frame(prev_abund_16S_HS_df, Prop=(Proportion*100))
prev_abund_16S_HS_df <- transform.data.frame(prev_abund_16S_HS_df, NumberofSamples=(Proportion*88))


abundancewindividual_16S_LS <- read.csv("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/RProjects/OTU_ra_16S_LS_windi.csv")
prev_abund_16S_LS_df <- cbind(prevdf_16S_LS, abundancewindividual_16S_LS)
colnames(prev_abund_16S_LS_df) <- c("Proportion","ASV","RelativeAbundance", "IndividualOnly")
prev_abund_16S_LS_df <- transform.data.frame(prev_abund_16S_LS_df, RA=(RelativeAbundance/88))
prev_abund_16S_LS_df <- transform.data.frame(prev_abund_16S_LS_df, Prop=(Proportion*100))
prev_abund_16S_LS_df <- transform.data.frame(prev_abund_16S_LS_df, NumberofSamples=(Proportion*88))

abundancewindividual_16S_HR <- read.csv("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/RProjects/OTU_ra_16S_HR_windi.csv")
prev_abund_16S_HR_df <- cbind(prevdf_16S_HR, abundancewindividual_16S_HR)
colnames(prev_abund_16S_HR_df) <- c("Proportion","ASV","RelativeAbundance", "IndividualOnly")
prev_abund_16S_HR_df <- transform.data.frame(prev_abund_16S_HR_df, RA=(RelativeAbundance/86))
prev_abund_16S_HR_df <- transform.data.frame(prev_abund_16S_HR_df, Prop=(Proportion*100))
prev_abund_16S_HR_df <- transform.data.frame(prev_abund_16S_HR_df, NumberofSamples=(Proportion*86))

abundancewindividual_16S_LR <- read.csv("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/Sequencing/RProjects/OTU_ra_16S_LR_windi.csv")
prev_abund_16S_LR_df <- cbind(prevdf_16S_LR, abundancewindividual_16S_LR)
colnames(prev_abund_16S_LR_df) <- c("Proportion","ASV","RelativeAbundance", "IndividualOnly")
prev_abund_16S_LR_df <- transform.data.frame(prev_abund_16S_LR_df, RA=(RelativeAbundance/88))
prev_abund_16S_LR_df <- transform.data.frame(prev_abund_16S_LR_df, Prop=(Proportion*100))
prev_abund_16S_LR_df <- transform.data.frame(prev_abund_16S_LR_df, NumberofSamples=(Proportion*88))

#### OK let's try this

scatter_16S_HS <- ggplot(prev_abund_16S_HS_df, aes(x=Prop,y=RA, color=IndividualOnly, alpha=IndividualOnly))+
  theme_bw()+
  geom_jitter()+
  labs(y="Relative Abundance (log10)", x="Proportion of Samples with ASVs Identified")+
  scale_y_log10(breaks=c(1e-6,1e-5,1e-3,1e-1,1e1,1e2), labels=c(-6,-5,-3,-1,1,""))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105))+
  scale_color_manual(values=c("black","red"))+
  scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "none",
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank())

sh_prev_abund_io_df <- subset(sh_prev_abund_df, IndividualOnly=="Y")
sh_prev_abund_all_df <- subset(sh_prev_abund_df, IndividualOnly=="N")

summary(sh_prev_abund_io_df$Prev)
sum(sh_prev_abund_io_df$Prev>1.137)#2043 in more than 1 sample, > 80% only in one sample
sum(sh_prev_abund_io_df$Prev>2.274) #806 in more than 2 samples, > 90 % only 2 or less samples
sum(sh_prev_abund_io_df$Prev>3.411) #374 in more than 3 samples, > 95 % in 3 or less samples
sum(sh_prev_abund_io_df$NumberofSamples>6) #54, 54/10372*100=0.52, less than 1% unshared ASVs in more than 10 % of animals sampled

summary(sh_prev_abund_all_df$Prevalence)

scatter_16S_LS <- ggplot(prev_abund_16S_LS_df, aes(x=Prop,y=RA, color=IndividualOnly, alpha=IndividualOnly))+
  theme_bw()+
  geom_jitter()+
  labs(y="Relative Abundance (log10)", x="Proportion of Samples with ASVs Identified")+
  scale_y_log10(breaks=c(1e-6,1e-5,1e-3,1e-1,1e1,1e2), labels=c(-6,-5,-3,-1,1,""))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105))+
  scale_color_manual(values=c("black","red"))+
  scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

sl_prev_abund_io_df <- subset(sl_prev_abund_df, IndividualOnly=="Y")
sl_prev_abund_all_df <- subset(sl_prev_abund_df, IndividualOnly=="N")

summary(sl_prev_abund_io_df$Prevalence)
summary(sl_prev_abund_all_df$Prevalence)

sum(sl_prev_abund_io_df$Prev>1.137) #739 in more than 1 sample, > 95% only in one sample
sum(sl_prev_abund_io_df$Prev>2.274) #198 in more than 2 samples, 99 % only 2 or less samples
sum(sl_prev_abund_io_df$Prev>3.411) #374 in more than 3 samples, > 99 % in 3 or less samples
sum(sl_prev_abund_io_df$Number)

#Plot binned with low bin containing 95% of all ASVs in individuals only (from code below)
scatter_16S_HR <- ggplot(prev_abund_16S_HR_df, aes(x=Prop,y=RA, color=IndividualOnly, alpha=IndividualOnly))+
  theme_bw()+
  geom_jitter()+
  labs(y="Relative Abundance (log10)", x="Proportion of Samples with ASVs Identified")+
  scale_y_log10(breaks=c(1e-6,1e-5,1e-3,1e-1,1e1,1e2), labels=c(-6,-5,-3,-1,1,""))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105))+
  scale_color_manual(values=c("black","red"))+
  scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

ggplot(rh_prev_abund_io_df, aes(x=Prev,y=RA, color=IndividualOnly,alpha=IndividualOnly))+
  theme_bw()+
  geom_jitter()+
  labs(y="Relative Abundance (log10)", x="Prevalence (%)")+
  scale_y_log10(breaks=c(1e-6,1e-5,1e-3,1e-1,1e1,1e2), labels=c(-6,-5,-3,-1,1,""))+
  scale_x_continuous(breaks = c(0,4,20,40,60,80,100))+
  scale_color_manual(values=c("red"))+
  scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "none",
        axis.text.y = element_text())

rh_prev_abund_io_df <- subset(rh_prev_abund_df, IndividualOnly=="Y")
rh_prev_abund_all_df <- subset(rh_prev_abund_df, IndividualOnly=="N")

summary(rh_prev_abund_io_df$Prevalence)
summary(rh_prev_abund_all_df$Prevalence)
sum(rh_prev_abund_io_df$Prev>1.164) #2049 in more than 1 sample, > 80% only in one sample
sum(rh_prev_abund_io_df$Prev>2.328) #793 in more than 2 samples, > 90 % only 2 or less samples
sum(rh_prev_abund_io_df$Prev>3.492) #366 in more than 3 samples, > 95 % in 3 or less samples
sum(rh_prev_abund_io_df$Prev>4) #366 in more than 3 samples, > 95 % in 3 or less samples

##Plot with prevalence cut off where 95% of ASVs only in individual are to the left of line---not sure which I like better
scatter_16S_LR <- ggplot(prev_abund_16S_LR_df, aes(x=Prop,y=RA, color=IndividualOnly, alpha=IndividualOnly))+
  theme_bw()+
  geom_jitter()+
  labs(y="Relative Abundance (log10)", x="Proportion of Samples with ASVs Identified")+
  scale_y_log10(breaks=c(1e-6,1e-5,1e-3,1e-1,1e1,1e2), labels=c(-6,-5,-3,-1,1,""))+
  scale_x_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105))+
  scale_color_manual(values=c("black","red"))+
  scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

rl_prev_abund_io_df <- subset(rl_prev_abund_df, IndividualOnly=="Y")
rl_prev_abund_all_df <- subset(rl_prev_abund_df, IndividualOnly=="N")

summary(rl_prev_abund_io_df$Prevalence)
summary(rl_prev_abund_all_df$Prevalence)
sum(rl_prev_abund_io_df$Prev>1.164) #735 in more than 1 sample, > 95% only in one sample
sum(rl_prev_abund_io_df$Prev>20)

###Look at genus level
upsetda_subset_high_genus <- get_upset(subseth_genus , factorNames="Number.in.Pool")

upset_sh_genus<-upset(upsetda_subset_high_genus  , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.8)
upset_sh_genus <- as_grob(upset_sh_genus)

upsetda_subset_low_genus <- get_upset(subsetl_genus , factorNames="Number.in.Pool")
upset_sl_genus<-upset(upsetda_subset_low_genus  , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9", keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.8)

upsetda_random_high_genus <- get_upset(randomh_genus , factorNames="Number.in.Pool")
upset_rh_genus<-upset(upsetda_random_high_genus  , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9", keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.8)

upsetda_random_low_genus <- get_upset(randoml_genus , factorNames="Number.in.Pool")
upset_rl_genus<-upset(upsetda_random_low_genus  , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.8)


combinedscatter <- ggarrange(scatter_16S_HS, scatter_16S_LS, scatter_16S_HR, scatter_16S_LR, ncol = 1)

ggsave("ASVscatter_20240602.tiff", path = "../DATA", plot = combinedscatter,device = "tiff", dpi = 600, width = 90, units = "mm")


#####MH UPSET PLOTS
#Subset
upsetda_MH_high <- get_upset(MH_high.css , factorNames="Number.in.Pool") ## ASV
MH_high_melt_asv <- psmelt(MH_high.css)

upsetda_MH_low <- get_upset(MH_low.css , factorNames="Number.in.Pool") ## ASV
MH_low_melt_asv <- psmelt(MH_low.css)

upset_MH_high <- upset(upsetda_MH_high , sets=c("1","3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.2)
upset_MH_low <- upset(upsetda_MH_low , sets=c("3" , "6", "12"),sets.bar.color = "#56B4E9",keep.order=T, order.by="freq", show.numbers="yes", mainbar.y.label="", text.scale=1.2)







