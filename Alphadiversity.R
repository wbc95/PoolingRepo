###Alphadiversity
pool_alpha_div1 <- estimate_richness(pool_16S, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
pool_alpha_div2 <- estimate_pd(pool_16S) # calculating Faith's PD

pool_alpha_div <- cbind(pool_alpha_div1, pool_alpha_div2) # combining Faith's and other metrics
pool_alpha_div # looks good, but have some duplicates so lets trim it a bit
pool_alpha_div <- pool_alpha_div[,c(1:5)]
pool_alpha_div
pool_alpha_div.df <- as(sample_data(pool_16S), "data.frame") # making into DF for metadata
pool_alpha_div_meta <- cbind(pool_alpha_div, pool_alpha_div.df)
pool_alpha_div_meta
min(pool_alpha_div_meta$Observed)
max(pool_alpha_div_meta$Observed)


indiv_alpha_div1 <- estimate_richness(indiv_16S, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # calculating most metrics
indiv_alpha_div2 <- estimate_pd(indiv_16S) # calculating Faith's PD

indiv_alpha_div <- cbind(indiv_alpha_div1, indiv_alpha_div2) # combining Faith's and other metrics
indiv_alpha_div # looks good, but have some duplicates so lets trim it a bit
indiv_alpha_div <- indiv_alpha_div[,c(1:5)]
indiv_alpha_div
indiv_alpha_div.df <- as(sample_data(indiv_16S), "data.frame") # making into DF for metadata
indiv_alpha_div_meta <- cbind(indiv_alpha_div, indiv_alpha_div.df)
indiv_alpha_div_meta
min(indiv_alpha_div_meta$Observed)
max(indiv_alpha_div_meta$Observed)

###
alpha_div_16S <- estimate_richness(data_16S_50k, measures = c("Observed","Shannon"))
alpha_div_16S.df <- as(sample_data(data_16S_50k),"data.frame")
alpha_div_16S_meta <- cbind(alpha_div_16S,alpha_div_16S.df)
min(alpha_div_16S_meta$Observed)


###ARGdata
ARG_high_alpha_div <- estimate_richness(ARG_high_trimmed, measures = c("Observed","Shannon"))
ARG_high_div.df <- as(sample_data(ARG_high_trimmed), "data.frame")
ARG_high_div_meta <- cbind(ARG_high_alpha_div, ARG_high_div.df)

ARG_low_alpha_div <- estimate_richness(ARG_low, measures = c("Observed","Shannon"))
ARG_low_div.df <- as(sample_data(ARG_low), "data.frame")
ARG_low_div_meta <- cbind(ARG_low_alpha_div, ARG_low_div.df)


#Poolsize-High
ARG_high_div_meta %>%
  wilcox_test(Observed~Contains.CP) #P=0.52
ARG_high_div_meta %>%
  kruskal_test(Observed~Number.in.Pool) #P=0.000494


ggplot(ARG_high_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "High", y= "Richness") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

kruskal_test(ARG_high_div_meta, Observed~Number.in.Pool) #P=0.000289, sig
dunn_test(ARG_high_div_meta, Observed~Number.in.Pool, p.adjust.method = "BH") #individuals=lower richness

ARG_high_div_meta %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Observed), median=median(Observed), n=n(),sd=sd(Observed), min=min(Observed),max=max(Observed))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

ARG_high_div_meta %>%
  group_by(Number.in.Pool)%>%
  reframe(mean=mean(Shannon), median=median(Shannon), n=n(),sd=sd(Shannon), min=min(Shannon),max=max(Shannon))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

hs_observed <- ggplot(ARG_high_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Richness") +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

kruskal_test(ARG_high_div_meta, Shannon~Number.in.Pool) #P=0.00116, NS
dunn_test(ARG_high_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH") #A lower

hs_shannon <- ggplot(ARG_high_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "high_trimmed", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

#Poolsize-Low
kruskal_test(ARG_low_div_meta, Observed~Number.in.Pool) #P=0.977, NS
kruskal_test(ARG_low_div_meta, Shannon~Number.in.Pool) #P=0.895, NS

ls_observed <- ggplot(ARG_low_div_meta, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "No. unique ARGs") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        title = element_blank(),       
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ls_shannon <- ggplot(ARG_low_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "No. unique ARGs") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggarrange(hs_observed, ls_observed, hs_shannon, ls_shannon, nrow=2, ncol=2)

kruskal_test(low_alpha_div_meta, Shannon~Number.in.Pool) #P=0.080, NS
dunn_test(low_alpha_div_meta, Shannon~Number.in.Pool, p.adjust.method = "BH")

kruskal_test(low_trimmed_alpha_div_meta, Shannon~Number.in.Pool) #P=0.863



ggplot(low_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggplot(low_trimmed_alpha_div_meta, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

##Compare differences faceted by poolsize


####MH####
alpha_div_MHhigh <- estimate_richness(MH_high, measures = c("Observed"))
MHhigh.df <- as(sample_data(MH_high), "data.frame")
alpha_div_MHhigh <- cbind(alpha_div_MHhigh, MHhigh.df)

alpha_div_MHlow <-estimate_richness(MH_low,measures = c("Observed"))
MHlow.df <- as(sample_data(MH_low),"data.frame")
alpha_div_MHlow <- cbind(alpha_div_MHlow, MHlow.df)

kruskal_test(alpha_div_MHhigh, Observed~Number.in.Pool) #NS, 0.491

ggplot(alpha_div_MHhigh, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "High_MH", y= "Richness") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggplot(alpha_div_MHhigh, aes(x= Number.in.Pool, y = Shannon, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "High_MH", y= "Shannon's Diversity Index") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = numberinpoolpalette) +
  scale_colour_manual(values = numberinpoolpalette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

kruskal_test(alpha_div_MHlow, Observed~Number.in.Pool) #NS, 0.831

ggplot(alpha_div_MHlow, aes(x= Number.in.Pool, y = Observed, fill =  Number.in.Pool, colour =  Number.in.Pool)) + 
  theme_bw() + 
  labs(title= "Low_MH", y= "Richness") +
  geom_boxplot(alpha = 0.5, size = 1) +
  geom_point(size = 2.5) +
  scale_fill_manual(values = poolonly_palette) +
  scale_colour_manual(values = poolonly_palette) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())


alpha_div_MHindiv <- estimate_richness(MH_indiv, measures = c("Observed","Shannon"))
MHindiv.df <- as(sample_data(MH_indiv), "data.frame")
alpha_div_MHindiv <- cbind(alpha_div_MHindiv, MHindiv.df)

wilcox_test(alpha_div_MHindiv, Observed~Mh.Cult) #Sig, 0.00172
wilcox_test(alpha_div_MHindiv, Shannon~Mh.Cult) #Sig, 0.00312

ggplot(alpha_div_MHindiv, aes(x= Number.in.Pool, y = Observed, fill =  Mh.Cult, colour =  Mh.Cult)) + 
  theme_bw() + 
  labs(title= "indiv_MH", y= "Richness") +
  geom_boxplot(alpha = 0.5, size = 1) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

ggplot(alpha_div_MHindiv, aes(x= Number.in.Pool, y = Shannon, fill =  Mh.Cult, colour =  Mh.Cult)) + 
  theme_bw() + 
  labs(title= "indiv_MH", y= "Richness") +
  geom_boxplot(alpha = 0.5, size = 1) +
  theme(#legend.position = "right",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 28, vjust = 1.75),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(),
    panel.border = element_rect(colour = "black", linewidth = 1.0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

alpha_div_MHindiv %>%
  group_by(Mh.Cult)%>%
  reframe(mean=mean(Observed), median=median(Observed), n=n(),sd=sd(Observed), min=min(Observed),max=max(Observed))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)

alpha_div_MHindiv %>%
  group_by(Mh.Cult)%>%
  reframe(mean=mean(Shannon), median=median(Shannon), n=n(),sd=sd(Shannon), min=min(Shannon),max=max(Shannon))%>%
  mutate(se=sd/sqrt(n), lower_ci=median-qt(1-(0.05/2),n-1)*se,
         upper_ci=median+qt(1-(0.05/2),n-1)*se)
