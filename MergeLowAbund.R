merge_low_abundance <- function(data, threshold=1){
  transformed <- transform_sample_counts(data, function(x) {x/sum(x)}*100)
  otu.table <- as.data.frame(otu_table(transformed))
  otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "zzzOther"
      tax_table(merged)[i,1:7] <- "zzzOther "}
  }
  return(merged)
}

###merge anything less than top code
merge_less_than_top <- function(data, top=20){
  transformed <- transform_sample_counts(data, function(x) {x/sum(x)}*100)
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "zzzOther"
      tax_table(merged)[i,1:7] <- "zzzOther"} # change to represent the number of levels
  }
  return(merged)
}

merge_less_than_top_amr <- function(data, top=20){
  transformed <- transform_sample_counts(data, function(x) {x/sum(x)}*100)
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "zzzOther"
      tax_table(merged)[i,1:5] <- "zzzOther"} # change to represent the number of levels
  }
  return(merged)
}
