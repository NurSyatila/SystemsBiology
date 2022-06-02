#============= PPI_Network ====================
#Function: Plot PPI network
#Usage: get_PPI_plot(geneNetwork,graphTitle,geneAnnotation)
# Arguments:
# - geneNetwork: Gene network for plotting (PPI network derived from STRINGdb)
# - graphTitle: Graph title
# - geneAnnotation: Gene annotation
# Value: graph plot
check_libraries <- function(pack_cran, pack_bioc) {
    for (pack in pack_cran){
        if (!pack %in% installed.packages()) {
            print(paste("installing",pack))
            #install.packages(pack)
            install.packages(pack,repos='http://cran.us.r-project.org')
        } else {print(paste(pack,"is already installed"))}
    }
    for (pack in pack_bioc){
        if (!pack %in% installed.packages()) {
            print(paste("installing",pack))
            if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
            BiocManager::install(pack,ask=FALSE)
        } else {print(paste(pack,"is already installed"))}
    }
    print("Loading required libraries, this might takes a while..")
    all_packages <- c(pack_cran,pack_bioc)
    for (pack in all_packages){
        suppressPackageStartupMessages({library(pack, character.only = TRUE)})
    }
}
get_ppi_plot <- function(geneNetwork,graphTitle,geneAnnotation){
  geneAnnotation_sorted <- geneAnnotation[match(V(geneNetwork)$name, geneAnnotation$V1),]
  V(geneNetwork)$name <- geneAnnotation_sorted$V2
  V(geneNetwork)$vertex.frame.color <- "white"
  edgeweights <- E(geneNetwork)$weight * 2.0
  suppressWarnings(plot(
    geneNetwork,
    layout=layout.fruchterman.reingold,
    edge.curved=FALSE,
    vertex.label.color="black",
    asp=FALSE,
    vertex.label.cex=0.6,
    vertex.shape="circle",
    vertex.color=rainbow(betweenness(geneNetwork), start=0, end=2/6),
    edge.width=edgeweights,
    edge.arrow.mode=0,
    main=graphTitle))
}

get_kegg_enrichment <- function(geneList,summaryFile) {
    ORAresults <- enrichKEGG(gene = geneList, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", use_internal_data = FALSE, minGSSize = 3,maxGSSize = 10000)
    sumORA <- as.data.frame(ORAresults)
    titlepic <- "**** Enriched KEGG Pathways ****"
    if (nrow(sumORA)>1){
        if (nrow(sumORA)>4){
          sumORA <- sumORA[1:5,c(1,2,3,6)]
          sumORA$p.adjust <- sprintf("%e",sumORA$p.adjust)
          sumORA2 <- knitr::kable(sumORA,row.names=TRUE,caption=titlepic,"simple")
          write(sumORA2,file=summaryFile,append=TRUE)
          write("\n",file=summaryFile,append=TRUE)
        } else {
          sumORA <- sumORA[,c(1,2,3,6)]
          sumORA$p.adjust <- sprintf("%e",sumORA$p.adjust)
          sumORA2 <- knitr::kable(sumORA,row.names=TRUE,caption=titlepic,"simple")
          write(sumORA2,file=summaryFile,append=TRUE)
          write("\n",file=summaryFile,append=TRUE)
        }
    } else {
        write(paste(titlepic, " - No results",sep=""),file=summaryFile,append=TRUE)
    }
}

get_gene_modules <- function(geneNetwork,summaryFile){
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
      suppressMessages(string_db$plot_network(unlist(str_split(df_clusters$genes[i],","),use.names=FALSE), add_link=FALSE, add_summary=TRUE))
    }
  }
  else {
    par(mfrow=c(2,2))
    for(i in seq_along(df_clusters$cluster_id)){
      suppressMessages(string_db$plot_network(unlist(str_split(df_clusters$genes[i],","),user.names=FALSE), add_link=FALSE, add_summary=TRUE))
    }
  }
  write(paste("\n** Number of clusters in main network: ",length(df_edges)," **",collapse=" "),file=summaryFile,append=TRUE)
  print ("Get pathways for individual clusters..
  ")
  suppressMessages(analyse_ppi_clusters(df_clusters,summaryFile))
  return (df_clusters)
}
get_enrichment_results <- function(geneList,ORAterm, summaryFile) {
    if (ORAterm=="BP"){
        ORAresults <- enrichGO(gene = geneList, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000)
        #ORAresults <- simplify(ORAresults, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    }
    if (ORAterm=="MF"){
        ORAresults <- enrichGO(gene = geneList, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "MF",pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000)
        #ORAresults <- simplify(ORAresults, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    }
    if (ORAterm=="CC"){
        ORAresults <- enrichGO(gene = geneList, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000)
        #ORAresults <- simplify(ORAresults, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    }
    if (ORAterm=="KEGG"){
        ORAresults <- enrichKEGG(gene = geneList, organism = "hsa", keyType = "ncbi-geneid", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", use_internal_data = FALSE, minGSSize = 3,maxGSSize = 10000)
    }
    if (ORAterm=="DO"){
        ORAresults <- enrichDO(gene = geneList, ont = "DO", pvalueCutoff = 0.01, qvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 10000)
    }
    cluster_summary <- data.frame(ORAresults)
    write.csv(cluster_summary , file = paste(ORAterm,"_Results.csv",sep=''), row.names=TRUE)
    if (ORAterm=="BP"){titlepic <- "**** Enriched GO Biological Process ****"}
    if (ORAterm=="MF"){titlepic <- "**** Enriched GO Molecular Function ****"}
    if (ORAterm=="CC"){titlepic <- "**** Enriched GO Cellular Component ****"}
    if (ORAterm=="KEGG"){titlepic <- "**** Enriched KEGG Pathways ****"}
    if (ORAterm=="DO"){titlepic <- "**** Enriched Disease Ontology ****"}
    sumORA <- as.data.frame(ORAresults)
    if (nrow(sumORA)>1){
        sumORA <- sumORA[1:10,c(1,2,3,6)]
        sumORA$p.adjust <- sprintf("%e",sumORA$p.adjust)
        sumORA2 <- knitr::kable(sumORA,row.names=TRUE,caption=titlepic,"simple")
        write(sumORA2,file=summaryFile,append=TRUE)
        write("\n",file=summaryFile,append=TRUE)
        pdf(file=paste(ORAterm,"_Plots.pdf",sep=''))
        plot1 <- suppressWarnings(barplot(ORAresults, showCategory=10))
        plot2 <- suppressWarnings(dotplot(ORAresults, showCategory=10))
        ORAresultsx <- setReadable(ORAresults, 'org.Hs.eg.db', 'ENTREZID')
        plot3 <- suppressWarnings(cnetplot(ORAresultsx,node_label="category", cex_label_category = 1.2))
        plot4 <- suppressWarnings(cnetplot(ORAresultsx, node_label="none",color_category='firebrick', color_gene='steelblue'))
        plot5 <- suppressWarnings(heatplot(ORAresultsx,showCategory=5))
        ORAresultsx2 <- pairwise_termsim(ORAresultsx)
        suppressWarnings(options(ggrepel.max.overlaps = Inf))
        plot6 <- suppressWarnings(emapplot(ORAresultsx2))
        print(plot1)
        print(plot2)
        print(plot3)
        print(plot4)
        print(plot5)
        print(plot6)
        dev.off()
        if (ORAterm=="KEGG"){
            print ("pathview")
            enriched_pathways <- sumORA[1:5,1]
            dir.create("Pathview")
            setwd("Pathview")
            for (pth_id in enriched_pathways){
                pathview(gene.data=geneList, pathway.id=pth_id, species = "hsa")
            }
            setwd("../")

        }
    } else {
        print(paste(titlepic, " - No results",sep=""))
        write(paste(titlepic, " - No results",sep=""),file=summaryFile,append=TRUE)
    }
}
analyse_ppi_clusters <- function(dfClusters,summaryFile){
  gz_filename <- "9606.protein.info.v11.5.txt.gz"
  annot_read <- readLines(gz_filename)
  annot_read2 <- strsplit(annot_read,'\t')
  annot_table <- lapply(annot_read2, function(x) {x[c(1,2)]})
  string_db <- suppressMessages(STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory="."))
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
      suppressMessages(get_kegg_enrichment(cluster_entrez$ENTREZID,summaryFile))
  }
}

analyse_ppi_network <- function(geneList){
  summaryFile <- "PPINetworkAnalysis.txt"
  options("download.file.method" = "libcurl")
  string_db <- suppressMessages(STRINGdb$new( version="11.5", species=9606, score_threshold=700, input_directory="."))
  write(paste("** Number of input genes: ",length(geneList),"gene(s) **",collapse=" "),file=summaryFile,append=TRUE)
  genelist <- as.data.frame(geneList)
  mapped_genes <- suppressWarnings(string_db$map(genelist, "geneList", removeUnmappedRows = TRUE, quiet = TRUE ))
  write(paste("** Number of mapped genes (Entrez -> STRINGdb ids): ",length(mapped_genes$STRING_id),"gene(s) **",collapse=" "),file=summaryFile,append=TRUE)
  print ("Get PPI network..")
  pdf('STRING_PPI_Network.pdf')
  suppressMessages(string_db$plot_network( mapped_genes$STRING_id, add_link=TRUE, add_summary=TRUE))

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
  suppressMessages(string_db$plot_network( filtered_genes3$ensembl_id, add_link=TRUE, add_summary=TRUE))
  print ("Functional enrichment analysis..")
  write(paste("** FilteredGenes: Functional enrichment using clusterProfiler [",length(filtered_genes3$ensembl_id),"gene(s)] (isolated nodes discarded) **\n",collapse=" "),file=summaryFile,append=TRUE)
  suppressWarnings(write.table(filtered_genes3,file="FilteredGenes.txt",quote=FALSE,append=TRUE, col.names = TRUE, row.names = FALSE, sep = "\t"))
  dir.create("FilteredGenes")
  setwd("FilteredGenes")
  #ORAnames <- c("BP","MF","CC","KEGG","DO")
  ORAnames <- c("BP","MF","CC","DO")
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
  suppressMessages(get_gene_modules(gene_network_new2,summaryFile))
  dev.off()
}


#gene_file <- "database_genes.txt"
#geneList <- readLines(gene_file)
#analyse_ppi_network(geneList)
#dir.create("enrichment_analysis")
#    setwd("enrichment_analysis")

cat("Please enter the path/full path to gene list file (Entrez IDs separated by newline): ");
gene_file <- readLines("stdin",n=1);
cat( "\n" )
if (file.exists(gene_file)){
    geneList <- readLines(gene_file)
    if (length(geneList)>3){
        dir.create("FunctionalEnrichment")
        setwd("FunctionalEnrichment")
        pack_cran = c("tools","ggplot2","umap","ggnewscale","igraph","stringr")
        pack_bioc = c("clusterProfiler","STRINGdb","enrichplot","pathview","org.Hs.eg.db","DOSE")
        print("Check required libraries and dependencies...")
        check_libraries(pack_cran, pack_bioc)
        analyse_ppi_network(geneList)
        setwd("../")
        print ("Analysis completed.")
    } else {print("Number of genes < 3. Analysis cannot be conducted.")}
} else {
    print ("The input gene file cannot be specified. Please enter the path/full path to the gene list file.")
}
