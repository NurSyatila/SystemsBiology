#============= Gene_Search: ================
#Function: Get genes for a given query disease in OMIM using rentrez package
#Usage: get_omim(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of disease-associated genes

get_omim <- function(queryTerms) {
   library(rentrez)
   library(stringr)
   library(org.Hs.eg.db)
   gene_file <- "Genes_OMIM.txt"
   pcgene_file <- "PCGenes_OMIM.txt"
   new_queryTerms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[DSDR]',sep='')
        new_queryTerms <- c(new_queryTerms, query_term)
    }
    new_queryTerms <- paste(new_queryTerms, collapse=' OR ')
    qTerms <- paste('(',new_queryTerms,')',sep='')
    dbSearch <- suppressMessages(entrez_search(db='omim',term=qTerms,use_history=TRUE,retmax=10000))
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
                      summary <- suppressMessages(entrez_summary(db="omim",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(extract_from_esummary(summary, c("uid","title","alttitles")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$alttitles,collapse=', '))
                         summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    }
                    dbLinks <- suppressMessages(entrez_link(dbfrom="omim", db="gene", cmd="neighbor_history", web_history=dbSearch$web_history))
                    if( !is.null(dbLinks$web_histories$omim_gene)) {
                      dbFetch <- suppressMessages(entrez_fetch(db="gene", web_history=dbLinks$web_histories$omim_gene,rettype="uilist",retmode="text"))
                      dbList <- str_split(dbFetch,"\n")
                      dbList <- unlist(dbList,recursive=TRUE)
                      dbList <- head(dbList,-1)
                    } else {dbList <- list()}
               } else {
                    summary <- suppressMessages(entrez_summary(db="omim",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(extract_from_esummary(summary, c("uid","title","alttitles")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$uid,collapse=', '),paste(esummaries[,i]$title,collapse=', '),paste(esummaries[,i]$alttitles,collapse=', '))
                          summarylist[[paste(esummaries[,i]$uid,collapse=', ')]] <- df
                      }
                    } else {
                        df <- data.frame(paste(esummaries$uid,collapse=', '),paste(esummaries$title,collapse=', '),paste(esummaries$alttitles,collapse=', '))
                        summarylist[[paste(esummaries$uid,collapse=', ')]] <- df
                    }
                    dbLinks <- suppressMessages(entrez_link(dbfrom="omim", db="gene", id=dbSearch$ids))
                    if( !is.null(dbLinks$links$omim_gene)) {
                      dbFetch <- suppressMessages(entrez_fetch(db="gene", id=dbLinks$links$omim_gene,rettype="uilist",retmode="text"))
                      dbList <- str_split(dbFetch,"\n")
                      dbList <- unlist(dbList,recursive=TRUE)
                      dbList <- head(dbList,-1)
                    } else {dbList <- list()}
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("Accession","Title","AltTitles")
                write.csv(summarytable, file="SummaryOMIM.csv",row.names = FALSE)
                print ("Get summary of query searches, output file: SummaryOMIM.csv ")
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){
                    Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = dbList, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "ENTREZID")
                    write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
                   print ("Get disease-related genes, output file: Genes_OMIM.txt ")
                   PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
                   PCGenes <- unique(PCGenes[,1])
                   if (length(PCGenes)>0){
                       MapGenes <- paste(PCGenes, collapse="\n")
                       write(MapGenes,file=pcgene_file)
                       print ("Get protein-coding genes, output file: PCGenes_OMIM.txt ")
                       MappedGenesx <- PCGenes
                   }
                   else {
                       print ("No disease-related protein-coding genes from OMIM.")
                       MappedGenesx <- list()
                   }

               } else {
                   print ("No disease-related genes from OMIM.")
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
       print ("No results from OMIM.")
       MappedGenesx <- list()

    }
   return(MappedGenesx)
}
#queryTerms <- c("colorectal cancer")
#queryTerms <- c("endometriosis")
#get_omim(queryTerms)
