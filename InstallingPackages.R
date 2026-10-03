##Keep this somewhere for installing packages after R update
install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("metagenomeSeq")
BiocManager::install("MicrobiotaProcess")
##Phyloseq all I need to do? Nopen retried library on other script...
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("metagMisc")
install_github("vmikk/metagMisc")
BiocManager::install("DESeq2") #Have to do this for btools
BiocManager::install("genefilter") #Have to do this for btools too
install_github("twbattaglia/btools")
install.packages("ggbreak")

BiocManager::install("microbiome")

install.packages("GUniFrac")
install_github("Russel88/MicEco")

