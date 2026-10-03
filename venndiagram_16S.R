###For venn diagrams

install.packages("ggVennDiagram")
library(ggVennDiagram)

install.packages("VennDiagram")
library(VennDiagram)

data_16S_50k #63542 taxa
data_16S_50k_indiv <- subset_samples(data_16S_50k, Number.in.Pool==1)
sum(taxa_sums(data_16S_50k_indiv)==0) #28129
data_16S_50k_indiv <- prune_taxa(taxa_sums(data_16S_50k_indiv)>0,data_16S_50k_indiv)
data_16S_50k_indiv #35413 taxa


data_16S_50k_pool <-subset_samples(data_16S_50k, Number.in.Pool!=1)
sum(taxa_sums(data_16S_50k_pool)==0) #18346
data_16S_50k_pool <- prune_taxa(taxa_sums(data_16S_50k_pool)>0,data_16S_50k_pool)
data_16S_50k_pool #45196

venn_16S_list = list(A=1:35413,B=18347:63542)

venn_16S_fig <- ggVennDiagram(venn_16S_list)

ggsave("venn16S_fig.tiff", plot=venn_16S_fig,device = "tiff", width = 100, units = "mm", dpi = 300, path = "Results/16S/Figures" )
