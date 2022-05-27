#============= Gene_Search: ================
#Function: Get genes for a given query disease in DrugBank
#Usage: get_DrugBank(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#XMLpath - path to drugnak XML
#Value: List of disease-associated genes

get_drugbank <- function(queryTerms) {
    suppressPackageStartupMessages({library(Autoseed)})
    suppressPackageStartupMessages({library(org.Hs.eg.db)})
    data("drugbank")
    gene_file <- "Genes_DrugBank.txt"
    pcgene_file <- "PCGenes_DrugBank.txt"
    drugbank_list <- character()
    summarylist <- list()
    for (qt in queryTerms){
        results_drugbank <- suppressMessages(drugbank_disease_gene(qt))
        if (identical(results_drugbank,character(0))==FALSE) {
            drugbank_symbols <- sapply(strsplit(results_drugbank,"::"), `[`, 1)
            drugbank_geneids <- AnnotationDbi::select(org.Hs.eg.db, keys = drugbank_symbols, c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
            drugbank_geneids2 <- drugbank_geneids[, 2]
            drugbank_list <- c(drugbank_list, drugbank_geneids2)
        }
    }
    drugbank_list2 <- na.omit(drugbank_list)
    drugbank_list2 <- unique(drugbank_list2)
    if (length(drugbank_list2)>0){
        drugbank_protein_coding_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = drugbank_list2, c("ENTREZID", "GENETYPE"), keytype = "ENTREZID")
        drugbank_protein_coding_genes <- dplyr::filter(drugbank_protein_coding_genes, grepl("protein",tolower(GENETYPE)))
        drugbank_protein_coding_genes <- unique(drugbank_protein_coding_genes[,1])
        if (length(drugbank_protein_coding_genes)>0){
            to_write_drugbank <- paste(drugbank_protein_coding_genes, collapse="\n")
            write(to_write_drugbank,file=gene_file)
            drugbank_protein_coding_genesx <- drugbank_protein_coding_genes
        } else {
            print ("No results from Drugbank.")
            drugbank_protein_coding_genesx <- list()
        }
    } else {
        print ("No results from drugbank.")
        drugbank_protein_coding_genesx <- list()
    }
    return(drugbank_protein_coding_genesx)
}
queryTerms <- c("endometriosis")
qt <- "endometriosis"



drugbank xml
curl -Lfv -o filename.zip -u syatila@mfrlab.org:135kemahiranhidup https://go.drugbank.com/releases/5-1-9/downloads/all-full-database
