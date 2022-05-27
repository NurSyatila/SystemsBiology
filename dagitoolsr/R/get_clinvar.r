#============= Gene_Search: ================
#Function: Get genes for a given query disease in Clinvar using rentrez package
#Usage: get_ClinVar(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of disease-associated genes

get_clinvar <- function(queryTerms) {
   library(rentrez)
   library(stringr)
   library(org.Hs.eg.db)
   gene_file <- "Genes_ClinVar.txt"
   pcgene_file <- "PCGenes_ClinVar.txt"
   new_query_terms <- character()
   for (var in queryTerms){
       query_term <- paste('"',var,'"[DIS]',sep='')
       new_query_terms <- c(new_query_terms, query_term)
   }
   qTerms <- paste('(((',new_query_terms,') AND 9606[TID]) NOT "clinsig conflicts"[FILT]) NOT "clinsig vus"[FILT]',sep='')
    dbSearch <- suppressMessages(entrez_search(db='clinvar',term=qTerms,use_history=TRUE,retmax=10000))
   if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing...")
               summarylist <- list()
               genelist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(entrez_summary(db="clinvar",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(extract_from_esummary(summary, c("uid","accession","title","clinical_significance","genes", "trait_set")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$clinical_significance$description,collapse=', '),paste(esummaries[,i]$genes$geneid,collapse=', '),paste(esummaries[,i]$trait_set$trait_name,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    }
                    dbLinks <- suppressMessages(entrez_link(dbfrom="clinvar", db="gene", cmd="neighbor_history", web_history=dbSearch$web_history))
                    if( !is.null(dbLinks$web_histories$clinvar_gene)) {
                      dbFetch <- suppressMessages(entrez_fetch(db="gene", web_history=dbLinks$web_histories$clinvar_gene,rettype="uilist",retmode="text"))
                      dbList <- str_split(dbFetch,"\n")
                      dbList <- unlist(dbList,recursive=TRUE)
                      dbList <- head(dbList,-1)
                    } else {dbList <- list()}
               } else {
                    summary <- suppressMessages(entrez_summary(db="clinvar",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(extract_from_esummary(summary, c("uid","accession","title","clinical_significance","genes", "trait_set")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$clinical_significance$description,collapse=', '),paste(esummaries[,i]$genes$geneid,collapse=', '),paste(esummaries[,i]$trait_set$trait_name,collapse=', '))
                        summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    } else {
                       df <- data.frame(paste(esummaries$uid,collapse=', '),paste(esummaries$accession,collapse=', '),paste(esummaries$title,collapse=', '),paste(esummaries$clinical_significance$description,collapse=', '),paste(esummaries$genes$geneid,collapse=', '),paste(esummaries$trait_set$trait_name,collapse=', '))
                        summarylist[[paste(esummaries$uid,collapse=', ')]] <- df
                    }
                    dbLinks <- suppressMessages(entrez_link(dbfrom="clinvar", db="gene", id=dbSearch$ids))
                    if( !is.null(dbLinks$links$clinvar_gene)) {
                      dbFetch <- suppressMessages(entrez_fetch(db="gene", id=dbLinks$links$clinvar_gene,rettype="uilist",retmode="text"))
                      dbList <- str_split(dbFetch,"\n")
                      dbList <- unlist(dbList,recursive=TRUE)
                      dbList <- head(dbList,-1)
                    } else {dbList <- list()}
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("UID","Accession","Title","Clinical Significance","GeneIDs","Traits")
                write.csv(summarytable, file="SummaryClinvar.csv",row.names = FALSE)
                print ("Get summary of query searches, output file: SummaryClinvar.csv ")
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){
                    Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = dbList, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "ENTREZID")
                    write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
                   print ("Get disease-related genes, output file: Genes_ClinVar.txt ")
                   PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
                   PCGenes <- unique(PCGenes[,1])
                   if (length(PCGenes)>0){
                       MapGenes <- paste(PCGenes, collapse="\n")
                       write(MapGenes,file=pcgene_file)
                       print ("Get protein-coding genes, output file: PCGenes_ClinVar.txt ")
                       MappedGenesx <- PCGenes
                   }
                   else {
                       print ("No disease-related protein-coding genes from ClinVar.")
                       MappedGenesx <- list()
                   }

               } else {
                   print ("No disease-related genes from ClinVar.")
                   MappedGenesx <- list()
               }
            }, error=function(e)
                {
                    print("Please check the output files.")
                    MappedGenesx <- list()
                }
                )
            }
            MappedGenesx <- getData()
   } else {
       print ("No results from ClinVar.")
       MappedGenesx <- list()

    }
   return(MappedGenesx)
}
#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#get_clinvar(queryTerms)
