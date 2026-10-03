###Calculating relative abundance on css
ARG.css <- phyloseq_transform_css(ARG_noSNP_ANCOMBC, log = F)
ARG.ra <- transform_sample_counts(ARG.css, function(x) {x/sum(x)}*100)
head(tax_table(ARG.ra))

ARG_type <- tax_glom(ARG.ra, taxrank = "Phylum")
ARG_type_melt <- psmelt(ARG_type)

ARG_class <- tax_glom(ARG.ra, taxrank = "Class")
ARG_class_melt <- psmelt(ARG_class)

ARG_SNV <- tax_glom(ARG.ra, taxrank = "Species")
ARG_SNV_melt <- psmelt(ARG_SNV)

ARG_class_ra_sd <- ARG_class_melt %>%
  group_by(Number.in.Pool,Class)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))

View(ARG_class_ra_sd)

ARG_SNV_ra_sd <- ARG_SNV_melt %>%
  group_by(Species)%>%
  reframe(meanRA=mean(Abundance),sd=sd(Abundance))

View(ARG_SNV_ra_sd)
