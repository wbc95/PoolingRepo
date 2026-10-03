####Venn Diagrams

#322 total taxa

ARG_indiv.ps <- subset_samples(data_trimmed_middle, Number.in.Pool==1)
sum(taxa_sums(ARG_indiv.ps)==0)
ARG_indiv.ps <- prune_taxa(taxa_sums(ARG_indiv.ps)>0, ARG_indiv.ps) #116 taxa

ARG_pools <- subset_samples(data_trimmed_middle, Number.in.Pool!=1)
sum(taxa_sums(ARG_pools)==0)
ARG_pools <- prune_taxa(taxa_sums(ARG_pools)>0, ARG_pools) #286 taxa

ARG_venn_list=list(A=1:116, B=37:322)

ggVennDiagram(ARG_venn_list)
