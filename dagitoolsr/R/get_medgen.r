#============= Gene_Search: ================
#Function: Get genes for a given query disease in MedGen using rentrez package
#Usage: get_medgen(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of disease-associated genes

get_medgen <- function(queryTerms) {
   library(rentrez)
   library(stringr)
   library(org.Hs.eg.db)
   gene_file <- "Genes_MedGen.txt"
   pcgene_file <- "PCGenes_MedGen.txt"
   new_query_terms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[TITL]',sep='')
        new_query_terms <- c(new_query_terms, query_term)
    }
    new_query_terms <- paste(new_query_terms, collapse=' OR ')
    qTerms <- paste('(',new_query_terms,')',sep='')
    dbSearch <- suppressMessages(entrez_search(db='medgen',term=qTerms,use_history=TRUE,retmax=10000))
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
                      summary <- suppressMessages(entrez_summary(db="medgen",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <-suppressMessages(extract_from_esummary(summary, c("uid","semanticid","semantictype","conceptid","title","definition")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$semanticid,collapse=', '),paste(esummaries[,i]$semantictype$value,collapse=', '),paste(esummaries[,i]$conceptid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$definition$value,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    }
                    dbLinks <- suppressMessages(entrez_link(dbfrom="medgen", db="gene", cmd="neighbor_history", web_history=dbSearch$web_history))
                    if( !is.null(dbLinks$web_histories$medgen_gene)) {
                        dbFetch <- suppressMessages(entrez_fetch(db="gene", web_history=dbLinks$web_histories$medgen_gene,rettype="uilist",retmode="text"))
                        dbList <- str_split(dbFetch,"\n")
                        dbList <- unlist(dbList,recursive=TRUE)
                        dbList <- head(dbList,-1)
                    } else { dbList <- list() }
               } else {
                    summary <- suppressMessages(entrez_summary(db="medgen",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <-suppressMessages(extract_from_esummary(summary, c("uid","semanticid","semantictype","conceptid","title","definition")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$semanticid,collapse=', '),paste(esummaries[,i]$semantictype$value,collapse=', '),paste(esummaries[,i]$conceptid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$definition$value,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    } else {
                       df <- data.frame(paste(esummaries$uid,collapse=', '),paste(esummaries$semanticid,collapse=', '),paste(esummaries$semantictype$value,collapse=', '),paste(esummaries$conceptid,collapse=', '),paste(esummaries$title,collapse=', '),paste(esummaries$definition$value,collapse=', '))
                         summarylist[[paste(esummaries$uid,collapse=', ')]] <- df
                    }
                    dbLinks <- suppressMessages(entrez_link(dbfrom="medgen", db="gene", id=dbSearch$ids))
                    if( !is.null(dbLinks$links$medgen_gene)) {
                        dbFetch <- suppressMessages(entrez_fetch(db="gene", id=dbLinks$links$medgen_gene,rettype="uilist",retmode="text"))
                        dbList <- str_split(dbFetch,"\n")
                        dbList <- unlist(dbList,recursive=TRUE)
                        dbList <- head(dbList,-1)
                    } else { dbList <- list() }
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("UID","SemanticID","SemanticType","ConceptID","Title","Definition")
                write.csv(summarytable, file="SummaryMedGen.csv",row.names = FALSE)
                print ("Get summary of query searches, output file: SummaryMedGen.csv ")
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){
                    Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = dbList, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "ENTREZID")
                    write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
                   print ("Get disease-related genes, output file: Genes_MedGen.txt ")
                   PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
                   PCGenes <- unique(PCGenes[,1])
                   if (length(PCGenes)>0){
                       MapGenes <- paste(PCGenes, collapse="\n")
                       write(MapGenes,file=pcgene_file)
                       print ("Get protein-coding genes, output file: PCGenes_MedGen.txt ")
                       MappedGenesx <- PCGenes
                   }
                   else {
                       print ("No disease-related protein-coding genes from MedGen.")
                       MappedGenesx <- list()
                   }

               } else {
                   print ("No disease-related genes from MedGen.")
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
       print ("No results from MedGen.")
       MappedGenesx <- list()
    }
   return(MappedGenesx)
}
#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#get_medgen(queryTerms)
