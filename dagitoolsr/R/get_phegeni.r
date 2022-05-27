#============= Gene_Search: ================
#Function: Get genes for a given query disease in PheGenI
#Usage: get_PheGenI(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
# - summaryFile: summary file to record output
#Value: List of disease-associated genes

get_phegeni <- function(queryTerms) {
    library(stringr)
    library(org.Hs.eg.db)
    #options(download.file.method="curl",timeout = 300)
    #options(download.file.method="curl")
    gene_file <- "Genes_PheGenI.txt"
    pcgene_file <- "PCGenes_PheGenI.txt"
    url_download <- 'http://www.ncbi.nlm.nih.gov/projects/gap/eqtl/EpiViewBE.cgi?type=dl.tab'
    suppressMessages(download.file(url_download, destfile='phegeni.tab'))
    if (file.exists('phegeni.tab') && file.size('phegeni.tab') > 0) {
        print ("Processing, get association... ")
        phegeni <- readLines(file('phegeni.tab'))
        phegeni_tab <- strsplit(phegeni,'\t')
        #phegeni_tab_short <- lapply(phegeni_tab, function(x) {x[c(1:8)]})
        patterns <- tolower(paste(queryTerms, collapse = "|"))
        summarylist <- Filter(function(x) grepl(patterns, tolower(x[2])), phegeni_tab)
        if (length(summarylist)>0){
            summarytable <- do.call(rbind,summarylist)
            colnames(summarytable) <- c("#", "Trait", "SNP rs", "Context", "Gene", "Gene ID", "Gene 2", "Gene ID 2", "Chromosome", "Location", "P-Value", "Source", "PubMed", "Analysis ID", "Study ID")
            write.csv(summarytable, file="SummaryPheGenI.csv",row.names = FALSE)
            print ("Get summary of query searches, output file: SummaryPheGenI.csv ")
            genelist <- unique(unlist(lapply(summarylist, function(x) {x[c(6,8)]})))
            genelist <- unique(unlist(genelist))
            Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = genelist, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "ENTREZID")
            write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
            print ("Get disease-related genes, output file: Genes_PheGenI.txt ")
            if (length(genelist)>0){
                PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
                PCGenes <- unique(PCGenes[,1])
                if (length(PCGenes)>0){
                       MapGenes <- paste(PCGenes, collapse="\n")
                       write(MapGenes,file=pcgene_file)
                       print ("Get protein-coding genes, output file: PCGenes_PheGenI.txt ")
                       MappedGenesx <- PCGenes
                   }
                   else {
                       print ("No disease-related protein-coding genes from PheGenI.")
                       MappedGenesx <- list()
                   }
            } else {
                print ("No disease-related genes from PheGenI.")
                MappedGenesx <- list()
            }
        } else {
            print ("No results from PheGenI.")
            MappedGenesx <- list()
        }
    } else {
        print ("Error: Unable to download file (ncbi/gap/phegeni).")
        MappedGenesx <- list()
    }
    return(MappedGenesx)
}


#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#[1] "#"           "Trait"       "SNP rs"      "Context"     "Gene"
# [6] "Gene ID"     "Gene 2"      "Gene ID 2"   "Chromosome"  "Location"
#[11] "P-Value"     "Source"      "PubMed"      "Analysis ID" "Study ID"
#[16] "Study Name"
#get_phegeni(queryTerms)