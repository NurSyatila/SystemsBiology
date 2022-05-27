#============= PPI_Network ====================
#Function: Analyse protein-protein interaction using STRINGdb (retrieve PPI, top genes and gene modules from clustering)
#Usage: analyse_PPI_network(geneList,summaryFile)
#Arguments:
# - geneList: input gene list (entrez identifiers separated by new line)
# - summaryFile: summary file to record output
# Value: none
analyse_ppi_network <- function(geneList){
  library(STRINGdb)
  library(org.Hs.eg.db)
  library(igraph)
  summaryFile <- "PPINetworkAnalysis.txt"
  options("download.file.method" = "libcurl")
  string_db <- suppressMessages(STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory="."))
  write(paste("** Number of input genes: ",length(geneList),"gene(s) **",collapse=" "),file=summaryFile,append=TRUE)
  genelist <- as.data.frame(geneList)
  mapped_genes <- suppressWarnings(string_db$map(genelist, "geneList", removeUnmappedRows = TRUE, quiet = TRUE ))
  write(paste("** Number of mapped genes (Entrez -> STRINGdb ids): ",length(mapped_genes$STRING_id),"gene(s) **",collapse=" "),file=summaryFile,append=TRUE)
  print ("Get PPI network..")
  pdf('STRING_PPI_Network.pdf')
  string_db$plot_network( mapped_genes$STRING_id, add_link=TRUE, add_summary=TRUE)

  mapped_genes <- suppressMessages(string_db$map(genelist, "geneList", removeUnmappedRows = TRUE ))
  gz_filename <- "9606.protein.info.v11.5.txt.gz"
  annot_read <- readLines(gz_filename)
  dat_annot <- as.data.frame(do.call(rbind, strsplit(annot_read, split="\t")))
  dat_annot <- dat_annot[,c(1,2)]
  #string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory=".")
  gene_network <- suppressMessages(string_db$get_subnetwork(mapped_genes$STRING_id))
  gene_network_new <- delete.vertices(gene_network, degree(gene_network)==0)
  gene_network_new2 <- igraph::simplify(gene_network_new, remove.multiple=TRUE, remove.loops=TRUE)
  filtered_genes <- as.data.frame(sort(degree(gene_network_new2),decreasing=TRUE))
  filtered_genes <- cbind(rownames(filtered_genes), data.frame(filtered_genes, row.names=NULL))
  colnames(filtered_genes) <- c("ensembl_id","degree")
  filtered_aliases <- dplyr::filter(dat_annot, grepl(paste(filtered_genes$ensembl_id,collapse = "|"),V1))
  colnames(filtered_aliases) <- c("ensembl_id","gene_symbol")

  filtered_genes2 <- merge(filtered_genes, filtered_aliases,all = TRUE)
  filtered_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = filtered_genes2$gene_symbol, c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
  colnames(filtered_entrez) <- c("gene_symbol","gene_id")
  filtered_genes3 <- merge(filtered_genes2, filtered_entrez,all = TRUE)
  filtered_genes3 <- filtered_genes3[order(filtered_genes3$degree,decreasing=TRUE),]

  write(paste("** Total genes from PPI network (isolated nodes discarded): ", length(filtered_genes3$ensembl_id)," **\n",sep=""),file=summaryFile,append=TRUE)
  string_db$plot_network( filtered_genes3$ensembl_id, add_link=TRUE, add_summary=TRUE)
  print ("Functional enrichment analysis..")
  write(paste("** FilteredGenes: Functional enrichment using clusterProfiler [",length(filtered_genes3$ensembl_id),"gene(s)] (isolated nodes discarded) **\n",collapse=" "),file=summaryFile,append=TRUE)
  suppressWarnings(write.table(filtered_genes3,file="FilteredGenes.txt",quote=FALSE,append=TRUE, col.names = TRUE, row.names = FALSE, sep = "\t"))
  dir.create("FilteredGenes")
  setwd("FilteredGenes")
  ORAnames <- c("BP","MF","CC","KEGG","DO")
  for (ORAname in ORAnames){
    print (paste("**", ORAname, "**", sep=" "))
    suppressMessages(get_enrichment_results(filtered_genes3$gene_id,ORAname,"../PPINetworkAnalysis.txt"))
  }
  setwd("../")
  if (length(filtered_genes$ensembl_id)>9){
    print ("Get top 10/20 genes..")
    if (length(filtered_genes$ensembl_id)>19){ top_genes_degrees <- head(as.data.frame(sort(degree(gene_network_new2),decreasing=TRUE)),20)}
    else {top_genes_degrees <- head(as.data.frame(sort(degree(gene_network_new2),decreasing=TRUE)),10)}
    top_genes_degrees <- cbind(rownames(top_genes_degrees), data.frame(top_genes_degrees, row.names=NULL))
    colnames(top_genes_degrees) <- c("ensembl_id","degree")
    top_aliases <- dplyr::filter(dat_annot, grepl(paste(top_genes_degrees$ensembl_id,collapse = "|"),V1))
    colnames(top_aliases) <- c("ensembl_id","gene_symbol")
    top_genes_degrees2 <- merge(top_genes_degrees, top_aliases,all = TRUE)
    top_entrez <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = top_genes_degrees2$gene_symbol, c("SYMBOL","ENTREZID"), keytype = "SYMBOL"))
    colnames(top_entrez) <- c("gene_symbol","gene_id")
    top_genes_degrees3 <- merge(top_genes_degrees2, top_entrez,all = TRUE)
    top_genes_degrees3 <- top_genes_degrees3[order(top_genes_degrees3$degree,decreasing=TRUE),]
    gene_network_top <- string_db$get_subnetwork(top_genes_degrees3$ensembl_id)
    suppressWarnings(get_ppi_plot(gene_network_top,"PPI Network (Top genes with the highest degree)",dat_annot))

    write("** Top genes from PPI network (with highest degree) **",file=summaryFile,append=TRUE)
    top_20_genes_table <- knitr::kable(top_genes_degrees3,row.names=TRUE,caption="List of top genes from PPI network","simple")
    suppressMessages(write(top_20_genes_table,file=summaryFile,append=TRUE))

    dir.create("TopGenes")
    setwd("TopGenes")
    write(paste(top_genes_degrees3$gene_symbol,collapse="\n"),file="../Top_PPI_genes.txt")
    write(paste("** Top genes: Functional enrichment using clusterProfiler **\n",collapse=" "),file="../PPINetworkAnalysis.txt",append=TRUE)
    for (ORAname in ORAnames){
      suppressMessages(get_enrichment_results(top_genes_degrees3$gene_id,ORAname,"../PPINetworkAnalysis.txt"))
    }
    setwd("../")
    #top_gene_id_list <- top_genes_degrees3$gene_symbol
  } else {
    write(paste("** Top genes cannot be generated **\n",collapse=" "),file=summaryFile,append=TRUE)
    #top_gene_id_list <- list()
  }
  print ("Get gene modules/clusters..")
  get_gene_modules(gene_network_new2,summaryFile)
  dev.off()
}

gene_file <- "database_genes.txt"
geneList <- readLines(gene_file)
analyse_ppi_network(geneList)