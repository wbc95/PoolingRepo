library(UpSetR)
library(MicrobiotaProcess)
library(MicEco)

#Import
upsetARG_high<- get_upset(ARG_high, factorNames="Number.in.Pool") ## ASV
upsetARG_low <- get_upset(ARG_low, factorNames="Number.in.Pool")


upsetARG_sh <- upset(upsetARG_high, sets = c("1","3","6","12"), sets.bar.color = "#56B4E9",order.by = "freq",keep.order=T, show.numbers="no", mainbar.y.label="", text.scale=1.8)
upsetARG_sl <- upset(upsetARG_low, sets = c("3","6","12"), sets.bar.color = "#56B4E9",order.by = "freq",keep.order=T, show.numbers="no", mainbar.y.label="", text.scale=1.8)

upsetARG_sh <- as_grob(upsetARG_sh)

upsetARG_both <- ggarrange(upsetARG_sh,upsetARG_sl, ncol = 1)

upsetARG_low <- upset(upsetARG_low, sets = c("1","3","6","12"), sets.bar.color = "#56B4E9",order.by = "freq",keep.order=T, show.numbers="no", mainbar.y.label="", text.scale=1.8)

#Compute prevalence of each feature, store as data.frame for each sample type 
prevdf_high_trimmed  <- as.data.frame(prevalence(rel_abund_high_trimmed, detection = 0))
#Extract relative abundance
abundancedf_high_trimmed <- as.data.frame(taxa_sums(rel_abund_high_trimmed))
write.csv(prevdf_high_trimmed , "OTU_prevalence_high_trimmed.csv") #Did this to check they were in same order
write.csv(abundancedf_high_trimmed, "OTU_ra_high_trimmed.csv") #Did this to check they were in the same order

##Idea..subset taxa, with RA=0...repeat each pool...will have to take this to Excel and manually add using =vlookup or =if..maybe too many rows
otus_high_trimmed <- as.data.frame(otu_table(rel_abund_high_trimmed))
colnames(otus_high_trimmed)
otus_high_trimmed<-subset(otus_high_trimmed, `H_12-01`==0&`HP_03-04`==0&`HP_03-01`==0&`HP_06-04`==0&`HP_06-01`==0&`HP_06-06`==0&`H_12-03`==0&`H_12-05`==0&`H_12-02`==0&`HP_03-02`==0&`HP_03-05`==0&`HP_03-03`==0&`HP_03-06`==0&`HP_06-05`==0&`H_12-06`==0)


##I think that worked
write.csv(otus_high_trimmed,"OTU_inviduals_high_trimmed.csv")

#combine df, used =countif and =if and copied values into column in csv with relative abundance
abundancewindividual_high_trimmed <- read.csv("~/Documents/Research/Projects/USDA_NIFA_AMR/Pooling Pilot/TE Pooling/OTU_ra_high_trimmed_windividuals.csv")
high_trimmed_prev_abund_df <- cbind(prevdf_high_trimmed,abundancewindividual_high_trimmed)
colnames(high_trimmed_prev_abund_df) <- c("Prevalence","ASV","RelativeAbundance", "IndividualOnly")
sum(high_trimmed_prev_abund_df$RelativeAbundance)
high_trimmed_prev_abund_df <- transform.data.frame(high_trimmed_prev_abund_df, RA=(RelativeAbundance/28))
high_trimmed_prev_abund_df <- transform.data.frame(high_trimmed_prev_abund_df, Prev=(Prevalence*100))
high_trimmed_prev_abund_df <- transform.data.frame(high_trimmed_prev_abund_df, NumberofSamples=(Prevalence*28))


#### OK let's try this

high_trimmed_scatter <- ggplot(high_trimmed_prev_abund_df, aes(x=Prev,y=RA, color=IndividualOnly, alpha=IndividualOnly))+
  theme_bw()+
  geom_jitter()+
  geom_jitter()+
  labs(y="Relative Abundance (log10)", x="Prevalence (%)")+
  scale_y_log10()+
  scale_x_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,105))+
  scale_color_manual(values=c("black","red"))+
  scale_alpha_manual(values=c(0.5, 1))+
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

high_trimmed_prev_abund_io_df <- subset(high_trimmed_prev_abund_df, IndividualOnly=="Y")
high_trimmed_prev_abund_all_df <- subset(high_trimmed_prev_abund_df, IndividualOnly=="N")




summary(high_trimmed_prev_abund_io_df$Prev)
sum(high_trimmed_prev_abund_io_df$Prev>3.571)#23 in more than 1 sample, all in more than one sample
sum(high_trimmed_prev_abund_io_df$Prev>7.142) #12 in more than 2 samples,< 50 % only 2 or less samples
sum(high_trimmed_prev_abund_io_df$Prev>10.713) #1 in more than 3 samples, >95 % in 3 or less samples
sum(sh_prev_abund_io_df$NumberofSamples>6) #54, 54/10372*100=0.52, less than 1% unshared ASVs in more than 10 % of animals sampled