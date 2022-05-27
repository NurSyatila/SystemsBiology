#============= Gene_Search: ================
#Function: Get genes for a given query disease in GWAS catalog using gwassrapid package
#Usage: get_GWAScatalog(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of disease-associated genes

get_gwascatalog <- function(queryTerms) {
    library(gwasrapidd)
    library(org.Hs.eg.db)
    gene_file <- "Genes_GWASCatalog.txt"
    pcgene_file <- "PCGenes_GWASCatalog.txt"
    gwas_list = character()
    summarylist <- list()
    print ("Processing, get association... ")
    for (qt in queryTerms){
        print (paste("Search for ",qt,collapse=""))
        my_associations <- get_associations(efo_trait = qt)
        if (length(my_associations@associations$association_id)>0){
            dplyr::filter(my_associations@associations, pvalue < 1e-6) %>% # Filter by p-value
            tidyr::drop_na(pvalue) %>%
            dplyr::pull(association_id) -> association_ids # Extract column association_id
            my_associations2 <- my_associations[association_ids]
            if (length(my_associations2@associations$association_id)>0){
                gwas_geneids <- my_associations2@entrez_ids$entrez_id # Get Entrez IDs
                gwas_list <- c(gwas_list, gwas_geneids)
                for (i in 1:length(my_associations2@associations$association_id)){
                    var1 <- my_associations2[i,]@associations$association_id
                    var2 <- my_associations2[i,]@associations$pvalue
                    var3 <- my_associations2[i,]@risk_alleles$variant_id
                    var4 <- my_associations2[i,]@entrez_ids$gene_name
                    var5 <- my_associations2[i,]@entrez_ids$entrez_id
                    getvar6 <- function(){
                        tryCatch({ var6 <- get_traits(association_id = my_associations2[125,]@entrez_ids$association_id)@traits$efo_id }, error=function(e)
                        { var6 <- "NA" })
                    }
                    getvar7 <- function(){
                        tryCatch({ var7 <- get_traits(association_id = my_associations2[i,]@entrez_ids$association_id)@traits$trait }, error=function(e)
                        { var6 <- "NA" })
                    }
                    var6 <- getvar6()
                    var7 <- getvar7()
                    df <- data.frame(paste(var1,collapse=', '),paste(var2,collapse=', '),paste(var3,collapse=', '),paste(var4,collapse=', '),paste(var5,collapse=', '),paste(var6,collapse=', '),paste(var7,collapse=', '))
                    summarylist[[paste(var1,collapse=', ')]] <- df
                }
            }
        }
    }
    if (length(summarylist) > 0){
        summarytable <- unique(do.call(rbind,summarylist))
        colnames(summarytable) <- c("AssociationID","pValue","VariantID","GeneName","EntrezID","EFOID","Trait")
        write.csv(summarytable, file="SummaryGWASCatalog.csv",row.names = FALSE)
        print ("Get summary of query searches, output file: SummaryGWASCatalog.csv ")
        gwas_list2 <- na.omit(gwas_list)
        gwas_list2 <- unique(gwas_list2)
        if (length(gwas_list2)>0){
            Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = gwas_list2, c("ENTREZID", "SYMBOL", "GENETYPE"), keytype = "ENTREZID")
            write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
            print ("Get disease-related genes, output file: Genes_GWASCatalog.txt ")
            PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
            gwas_protein_coding_genes <- unique(PCGenes[,1])
            if (length(gwas_protein_coding_genes)>0){
                print ("Get protein-coding genes, output file: PCGenes_GWASCatalog.txt ")
                to_write_gwas <- paste(gwas_protein_coding_genes, collapse="\n")
                write(to_write_gwas,file=pcgene_file)
                gwas_protein_coding_genesx <- gwas_protein_coding_genes
            } else {
                print ("No disease-related protein-coding genes from GWAS Catalog.")
                gwas_protein_coding_genesx <- list()
            }
        } else {
            print ("No disease-related genes from GWAS Catalog.")
            gwas_protein_coding_genesx <- list()
        }
    } else {
        print ("No results from GWAS Catalog.")
        gwas_protein_coding_genesx <- list()
    }

    return(gwas_protein_coding_genesx)
}
#queryTerms <- c("endometriosis")
queryTerms <- c("colorectal cancer")
#,"endometrioma")
#get_gwascatalog(queryTerms)
#https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId=28530673
#https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByEfoTrait?efoTrait=colorectal%20cancer