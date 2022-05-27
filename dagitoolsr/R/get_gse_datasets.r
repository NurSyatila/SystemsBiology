#============= Differential_Expression ====================
#Function: Get GEO datasets for given MeSH terms
#Usage: get_GEO_datasets(MeSHterms)
#Arguments:
# - MeSHterms = query in list e.g. c("endometriosis"), or c("colorectal cancer","colorectal carcinoma")
#Value: List of GEO accession numbers

#MeSHTerms <- get_meshterms(queryTerms,"term")

get_gse_datasets <- function(MeSHTerms) {
    new_query_terms <- character()
    for (var in MeSHTerms){
        query_term <- paste('"',var,'"[MESH]',sep='')
        new_query_terms <- c(new_query_terms, query_term)
    }
    new_query_terms <- paste(new_query_terms, collapse=' OR ')
    qTerms <- paste('(',new_query_terms,') AND "Homo sapiens"[ORGN] AND "gds"[FILT]',sep='')
    dbSearch <- entrez_search(db='gds',term=qTerms,use_history=TRUE,retmax=10000)
    if (length(dbSearch$ids)>0){
        getData <- function(){
           tryCatch({
               print ("Processing...")
               summarylist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(entrez_summary(db="gds",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(extract_from_esummary(summary, c("accession","gpl","gse","seriestitle","subsetinfo")))
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$gpl,collapse=', '),paste(esummaries[,i]$gse,collapse=', '),paste(esummaries[,i]$seriestitle,collapse=', '),paste(esummaries[,i]$subsetinfo,collapse=', '))
                          summarylist[[paste(esummaries[,i]$accession,collapse=', ')]] <- df
                      }
                    }
               } else {
                    summary <- suppressMessages(entrez_summary(db="gds",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(extract_from_esummary(summary, c("accession","gpl","gse","seriestitle","subsetinfo")))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          df <- data.frame(paste(esummaries[,i]$accession,collapse=', '),paste(esummaries[,i]$gpl,collapse=', '),paste(esummaries[,i]$gse,collapse=', '),paste(esummaries[,i]$seriestitle,collapse=', '),paste(esummaries[,i]$subsetinfo,collapse=', '))
                          summarylist[[paste(esummaries[,i]$accession,collapse=', ')]] <- df
                      }
                    } else {
                       df <- data.frame(paste(esummaries$accession,collapse=', '),paste(esummaries[,i]$gpl,collapse=', '),paste(esummaries$gse,collapse=', '),paste(esummaries$seriestitle,collapse=', '),paste(esummaries$subsetinfo,collapse=', '))
                          summarylist[[paste(esummaries$accession,collapse=', ')]] <- df
                    }
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("accession","gpl","gse","seriestitle","subsetinfo")
                summarytable$gpl <- paste("GPL",summarytable$gpl,sep="")
                summarytable$gse <- paste("GSE",summarytable$gse,sep="")
                write.csv(summarytable, file="SummaryGSE.csv",row.names = FALSE)
                Filteredsummarytable <- dplyr::filter(summarytable, grepl("disease",tolower(summarytable$subsetinfo)))
                write.csv(Filteredsummarytable, file="FilteredGSEDiseases.csv",row.names = FALSE)
                print ("Get summary of query searches, output file: SummaryGSE.csv, Filtered output file: FilteredGSEDiseases.csv ")
                GSElist <- unlist(Filteredsummarytable$gse,recursive=TRUE)
            }, error=function(e)
                {
                    print("Please check the output files.")
                    GSElist <- list()
                }
                )
            }
            GSElist <- getData()
   } else {
       print ("No results from GEO Datasets.")
       GSElist <- list()
    }
    return (GSElist)
}

queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
MeSHTerms <- get_meshterms(queryTerms,"term")
get_gse_datasets(MeSHTerms)
