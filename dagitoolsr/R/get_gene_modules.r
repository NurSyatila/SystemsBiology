#============= PPI_Network ====================
#Function: Get gene modules from network clustering
#Usage: get_gene_modules(geneNetwork,summaryFile)
# - geneNetwork: Gene network for plotting (PPI network derived from STRINGdb)
# - summaryFile: summary file to record output
#Value: dfClusters

get_gene_modules <- function(geneNetwork,summaryFile){
  library(igraph)
  library(STRINGdb)
  library(stringr)
  print ("Get modules from clustering..")
  string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory=".")
  g <- geneNetwork
  c <- walktrap.community(g,modularity=TRUE)
  c_keep_ids <- as.numeric(names(sizes(c)[sizes(c) >= 4]))
  c_keep_v_idxs <- which(c$membership %in% c_keep_ids)

  g_sub <- induced_subgraph(g, V(g)[c_keep_v_idxs])
  c_sub <- c
  c_sub$names <- c$names[c_keep_v_idxs]
  c_sub$membership <- c$membership[c_keep_v_idxs]
  c_sub$vcount <- length(c_sub$names)
  c_sub$modularity <- modularity(g_sub, c_sub$membership, E(g_sub)$weight)
  #sort by edges
  df_levels <- unique(membership(c_sub))
  df_edges <- list()
  df_genes <- list()
  for (x in seq_along(df_levels)){
	edges <- gsize(induced_subgraph(g_sub,unlist(c_sub[x],use.names=FALSE)))
	df_edges[[names(c_sub[x])]] <- edges
    df_genes[[names(c_sub[x])]] <- paste(unlist(c_sub[x],use.names=FALSE),collapse=",")
  }
  df_clusters1 <- data.frame(cluster_id=names(df_edges),edges=unlist(df_edges,use.names=FALSE))
  df_clusters2 <- data.frame(cluster_id=names(df_edges),genes=unlist(df_genes,use.names=FALSE))
  df_clusters <- merge(df_clusters1,df_clusters2,by="cluster_id")
  df_clusters <- df_clusters[order(df_clusters$edges,decreasing=TRUE),]
  if (nrow(df_clusters) > 3) {
    par(mfrow=c(2,2))
    for(i in seq(1:4)){
      string_db$plot_network(unlist(str_split(df_clusters$genes[i],","),use.names=FALSE), add_link=FALSE, add_summary=TRUE)
    }
  }
  else {
    par(mfrow=c(2,2))
    for(i in seq_along(df_clusters$cluster_id)){
      string_db$plot_network(unlist(str_split(df_clusters$genes[i],","),user.names=FALSE), add_link=FALSE, add_summary=TRUE)
    }
  }
  write(paste("\n** Number of clusters in main network: ",length(df_edges)," **",collapse=" "),file=summaryFile,append=TRUE)
  print ("Get pathways for individual clusters..
  ")
  analyse_ppi_clusters(df_clusters,summaryFile)
  return (df_clusters)
}