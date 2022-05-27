#============= Differential_Expression ====================
#Function: Analyse sample groupings and perform differential expression analysis
#Usage: get_DEG(GSEaccession)
#Arguments:
# - GSEaccession: GSE dataset accession number
# Value: List of differentially expressed genes from analyse_DEG function

#MeSHTerms <- get_meshterms(queryTerms,"term")
#GSEids <- get_gse_datasets(MeSHTerms)
#1] "GSE25628" "GSE23339" "GSE7846"  "GSE6364"  "GSE7305"  "GSE5108"
#GSEaccession <- "GSE25628"

get_deg <- function(GSEaccession) {
    library(GEOquery)
    summaryFile <- paste(GSEaccession,"_Analysis.txt",sep="")
    parse_gsefile <- function(GSEaccession){
        tryCatch({
            gset <- getGEO(GSEaccession, GSEMatrix =TRUE, AnnotGPL=TRUE,destdir=".")
            gse_file <- paste(GSEaccession,"_series_matrix.txt.gz",sep="")
            if (file.exists(gse_file)){
                gselines <- readLines(gse_file)
                platformid <- gselines[c(grep(pattern = "!Series_platform_id", x = gselines))]
                platformid <- gsub('"','', str_split(platformid[1],"\t")[[1]][2])
                if (length(gset) > 1) idx <- grep(platformid, attr(gset, "names")) else idx <- 1
                gset <- gset[[idx]]
                samples <- pData(gset)$characteristics_ch1
                groups <- unique(samples)
                print (GSEaccession)
                print (groups)
                if (length(groups)==2){
                    sample_groupings <- gsub(groups[1],"0",samples)
                    sampleGrouping <- gsub(groups[2],"1",sample_groupings)
                    groupDescription <- paste("--GroupA: ",groups[1],"\n","--GroupB: ",groups[2],sep="")
                    degsx_list <- analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile)
                } else {
                    patterns <- tolower(paste(query_terms, collapse = "|"))
                    new_groups <- grepl(patterns, tolower(groups))
                    group1 <- character()
                    group2 <- character()
                    for(i in seq(1:length(new_groups))){
                        if (new_groups[i]==TRUE){group1 <- c(group1,groups[i])}
                        if (new_groups[i]==FALSE){group2 <- c(group2,groups[i])}
                    }
                    conds_normal <- c("healthy","normal","control")
                    if ((length(group1) !=0 && length(group2) !=0) && unique(grepl(paste(conds_normal,collapse="|"), tolower(group2)))==TRUE){
                        sample_groupings <- gsub(paste(group1,collapse="|"),"0",samples)
                        sampleGrouping <- gsub(paste(group2,collapse="|"),"1",sample_groupings)
                        groupDescription <- paste("--GroupA: ",paste(group1,collapse=" | "),"\n","--GroupB: ",paste(group2,collapse=" | "),sep="")
                        degsx_list <- analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile)
                    }
                    else {
                        new_groups <- grepl("disease", tolower(groups))
                        group1 <- character()
                        group2 <- character()
                        for(i in seq(1:length(new_groups))){
                            if (new_groups[i]==TRUE){group1 <- c(group1,groups[i])}
                            if (new_groups[i]==FALSE){group2 <- c(group2,groups[i])}
                        }
                        conds_normal <- c("healthy","normal","control")
                        if ((length(group1) !=0 && length(group2) !=0) && unique(grepl(paste(conds_normal,collapse="|"), tolower(group2)))==TRUE){
                            sample_groupings <- gsub(paste(group1,collapse="|"),"0",samples)
                            sampleGrouping <- gsub(paste(group2,collapse="|"),"1",sample_groupings)
                            groupDescription <- paste("--GroupA: ",paste(group1,collapse=" | "),"\n","--GroupB: ",paste(group2,collapse=" | "),sep="")
                            degsx_list <- analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile)
                        } else {
                            print ("Error: More/less than two sample types detected. The GSE data set will not be used.")
                            write(paste("** ",GSEaccession,": More than two sample types detected. The GSE data set will not be used. **\n--GSE samples:",collapse=" "),file=summaryFile,append=TRUE)
                            write(paste(groups,collapse=" | "),file=summaryFile,append=TRUE)
                            degsx_list <- list()
                        }
                    }
                }
            } else {
                print("Error (2): File not exists.")
                write(paste("** ",GSEaccession,": Error - Unable to download / parse the GSE file.",collapse=" "),file=summaryFile,append=TRUE)
                degsx_list <- list()
            }
        }, error=function(e)
                {
                    print("Error (1): Unable to download / parse the GSE file.")
                    write(paste("** ",GSEaccession,": Error - Unable to download / parse the GSE file.",collapse=" "),file=summaryFile,append=TRUE)
                    degsx_list <- list()
                }
        )
    }
    deg_list <- parse_gsefile(GSEaccession)
    return (deg_list)
}
