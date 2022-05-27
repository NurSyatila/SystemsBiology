#============= Gene_Search: ================
#Function: Get MeSH terms for a given query disease
#Usage: get_MeSH_terms(queryTerms)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of MeSH terms

get_meshterms <- function(queryTerms) {
    library(rentrez)
    new_queryTerms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[MESH]',sep='')
        new_queryTerms <- c(new_queryTerms, query_term)
    }
    new_queryTerms <- paste(new_queryTerms, collapse=' OR ')
    qTerms <- paste('(("main heading"[TYPE]) AND (',new_queryTerms,')',sep='')
    dbSearch <- suppressMessages(entrez_search(db='mesh',term=qTerms,use_history=TRUE,retmax=10000))
    if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing...")
               MESHList <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(entrez_summary(db="mesh",id=unlist(chunks[1],recursive=TRUE)))
                      #MESHList <- c(MESHList, summary$ds_meshterms)
                      MESHList <- c(MESHList, gsub('^68','D',summary$uid))
                    }
               } else {
                    summary <- suppressMessages(entrez_summary(db="mesh",id=unlist(dbSearch$ids,recursive=TRUE)))
                    MESHList <- c(MESHList, gsub('^68','D',summary$uid))
                }
                MESHList <- unique(unlist(MESHList,recursive=TRUE))
            }, error=function(e)
                {
                    MESHList <- list()
                }
                )
            }
            MESHList <- getData()
   } else {
       print ("MESH terms cannot be identified.")
       MESHList <- list()

    }
    return(MESHList)
}

get_meshterms <- function(queryTerms,type) {
    #type == "term","id"
    library(rentrez)
    new_queryTerms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[MESH]',sep='')
        new_queryTerms <- c(new_queryTerms, query_term)
    }
    new_queryTerms <- paste(new_queryTerms, collapse=' OR ')
    qTerms <- paste('(("main heading"[TYPE]) AND (',new_queryTerms,')',sep='')
    dbSearch <- suppressMessages(entrez_search(db='mesh',term=qTerms,use_history=TRUE,retmax=10000))
    if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing...")
               MESHids <- list()
               MESHterms <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(entrez_summary(db="mesh",id=unlist(chunks[1],recursive=TRUE)))
                      MESHids <- c(MESHids, gsub('^68','D',summary$uid))
                      MESHterms <- c(MESHterms, summary$ds_meshterms)
                    }
               } else {
                    summary <- suppressMessages(entrez_summary(db="mesh",id=unlist(dbSearch$ids,recursive=TRUE)))
                    MESHids <- c(MESHids, gsub('^68','D',summary$uid))
                    MESHterms <- c(MESHterms, summary$ds_meshterms)
                }
                MESHids <- unique(unlist(MESHids,recursive=TRUE))
                MESHterms <- unique(unlist(MESHterms,recursive=TRUE))
                MESHList <- list(MESHids,MESHterms)
            }, error=function(e)
                {
                    MESHList <- list()
                }
                )
            }
            MESHList <- getData()
   } else {
       print ("MESH terms cannot be identified.")
       MESHList <- list()
    }
    if (type=="id"){return(MESHList[[1]])}
    if (type=="term"){return(MESHList[[2]])}
}

#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#MeSHTerms <- get_meshterms(queryTerms,"term")
#MeSHTerms <- get_meshterms(queryTerms,"id")
