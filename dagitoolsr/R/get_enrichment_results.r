#============= Functional_Enrichment====================
#library(dagitoolsr)
#Function: : Perform functional enrichment and generate visualization plots for a given gene list
#Usage: get_enrichment_results(geneList,ORAterm, summaryFile)
#Arguments:
# - geneList: input gene list (entrez identifiers separated by new line)
# - ORAterm: Functional enrichment terms ("BP","MF","CC","KEGG","DO")
# - summaryFile: summary file to record output
# Value: none
get_enrichment_results <- function(geneList,ORAterm, summaryFile) {
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(umap)
    library(ggnewscale)
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
        library(DOSE)
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
            library(pathview)
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