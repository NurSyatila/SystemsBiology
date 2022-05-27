#============= PPI_Network ====================
#Function: Perform KEGG enrichment given a gene list using ClusterProfiler
#Usage: get_KEGG_enrichment(geneList,summaryFile)
# - geneList: input gene list (entrez identifiers separated by new line)
# - summaryFile: summary file to record output
#Value: none

get_kegg_enrichment <- function(geneList,summaryFile) {
    library(clusterProfiler)
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
