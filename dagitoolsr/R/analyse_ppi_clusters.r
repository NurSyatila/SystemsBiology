#============= PPI_Network ====================
#Function: Analyse gene modules resulting from PPI clustering
#Usage: get_PPI_clusters(dfClusters)
#Arguments:
# - dfClusters: Data frame of clusters derived from clustering
# - summaryFile: summary file to record output
# Value: none

analyse_ppi_clusters <- function(dfClusters,summaryFile){
  library(STRINGdb)
  library(org.Hs.eg.db)
  library(stringr)
  gz_filename <- "9606.protein.info.v11.5.txt.gz"
  annot_read <- readLines(gz_filename)
  annot_read2 <- strsplit(annot_read,'\t')
  annot_table <- lapply(annot_read2, function(x) {x[c(1,2)]})
  string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory=".")
  dir.create("GeneClusters")
  for(i in seq_along(dfClusters$cluster_id)){
      clustergenes <- unlist(str_split(dfClusters$genes[i],","),use.names=FALSE)
      cluster_aliases <- Filter(function(x) grepl(paste(clustergenes,collapse = "|"), x[1]), annot_table)
      cluster_aliases2 <- unique(unlist(lapply(cluster_aliases, function(x) {x[2]})))
      cluster_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = cluster_aliases2, c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
      write(paste("**** CLUSTER ",i," (",length(cluster_aliases2)," genes"," )"," ****",sep=""),file=summaryFile,append=TRUE)
      write(paste(cluster_aliases2,collapse="\n"),file=paste("GeneClusters/cluster_",i,".txt",sep=""))
      write(paste("Gene Symbols: ",paste(cluster_aliases2,collapse=", "),sep=""),file=summaryFile,append=TRUE)
      write(paste("Gene IDs: ",paste(cluster_entrez$ENTREZID,collapse=", "),sep=""),file=summaryFile,append=TRUE)
      get_kegg_enrichment(cluster_entrez$ENTREZID,summaryFile)
  }

}