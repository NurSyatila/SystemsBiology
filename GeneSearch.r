#============= Gene_Search: ================
#Function: Get genes for a given query disease in databases
#Usage: GeneSearch(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of disease-associated genes
#Check libraries
check_libraries <- function(pack_cran, pack_bioc, pack_github) {
    for (pack in pack_cran){
        if (!pack %in% installed.packages()) {
            print(paste("installing",pack))
            #install.packages(pack)
            install.packages(pack,repos='http://cran.us.r-project.org')
        } else {print(paste(pack,"is already installed"))}
    }
    for (pack in pack_bioc){
        if (!pack %in% installed.packages()) {
            print(paste("installing",pack))
            if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager",,repos='http://cran.us.r-project.org')
            BiocManager::install(pack,ask=FALSE)
        } else {print(paste(pack,"is already installed"))}
    }
    for (pack in pack_github){
        if (!pack %in% installed.packages()) {
            if (pack == "disgenet2r"){
                print(paste("installing",pack))
                library(devtools)
                install_bitbucket("ibi_group/disgenet2r")
            }
            if (pack == "gwasrapidd"){
                print(paste("installing",pack))
                install.packages("remotes",repos='http://cran.us.r-project.org')
                remotes::install_github("ramiromagno/gwasrapidd")
            }
        } else {print(paste(pack,"is already installed"))}
    }
    print("Loading required libraries, this might takes a while..")
    all_packages <- c(pack_cran,pack_bioc,pack_github)
    for (pack in all_packages){
        suppressPackageStartupMessages({library(pack, character.only = TRUE)})
    }
}

get_PCGenes <- function(genelist, gene_file,pcgene_file,key_type){
  Genes <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = genelist, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = key_type))
  suppressWarnings(write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE))
  PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
  PCGenes <- unique(PCGenes$ENTREZID)
  if (length(PCGenes)>0){
    MapGenes <- paste(PCGenes, collapse="\n")
    write(MapGenes,file=pcgene_file)
    return (PCGenes)
  } else {
    PCGenes <- list()
    return (PCGenes)
  }
}

get_clinvar <- function(queryTerms) {
   suppressPackageStartupMessages(library(rentrez))
   suppressPackageStartupMessages(library(stringr))
   suppressPackageStartupMessages(library(org.Hs.eg.db))
   gene_file <- "ClinVar_Genes.txt"
   pcgene_file <- "ClinVar_PCGenes.txt"
   new_query_terms <- character()
   for (var in queryTerms){
       query_term <- paste('"',var,'"[DIS]',sep='')
       new_query_terms <- c(new_query_terms, query_term)
   }
   new_query_terms <- paste(new_query_terms, collapse=' OR ')
   qTerms <- paste('(((',new_query_terms,') AND 9606[TID]) NOT "clinsig conflicts"[FILT]) NOT "clinsig vus"[FILT]',sep='')
    dbSearch <- suppressMessages(entrez_search(db='clinvar',term=qTerms,use_history=TRUE,retmax=10000))
   if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (ClinVar)...")
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
                      dbList <- strsplit(dbFetch,"\n")
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
                      dbList <- strsplit(dbFetch,"\n")
                      dbList <- unlist(dbList,recursive=TRUE)
                      dbList <- head(dbList,-1)
                    } else {dbList <- list()}
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("UID","Accession","Title","Clinical Significance","GeneIDs","Traits")
                write.csv(summarytable, file="ClinvarSummary.csv",row.names = FALSE)
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){ MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID") }
               else {
                 print ("No results from Clinvar.")
                 MappedGenesx <- list() }
            }, error=function(e)
                {
             print ("Error: Please check the output files.")
             MappedGenesx <- list()
           }) }
            MappedGenesx <- getData()
   } else {
       print ("No results from ClinVar.")
       MappedGenesx <- list()
    }
   return(MappedGenesx)
}
get_diseasesdb <- function(queryTerms) {
    gene_file <- "DISEASES_Genes.txt"
    pcgene_file <- "DISEASES_PCGenes.txt"
    #options(download.file.method="curl",timeout = 300)
    #options(download.file.method="curl")
    patterns <- paste(queryTerms, collapse = "|")
    print ("Processing DISEASES (Text Mining category)...")
    url_text_mining <- 'http://download.jensenlab.org/human_disease_textmining_filtered.tsv'
    suppressMessages(download.file(url_text_mining, destfile='DISEASES_human_disease_textmining_filtered.tsv'))
    if (file.exists('DISEASES_human_disease_textmining_filtered.tsv') && file.size('DISEASES_human_disease_textmining_filtered.tsv') > 0) {
        text_table <- read.csv(file = 'DISEASES_human_disease_textmining_filtered.tsv', sep = '\t',header=FALSE)
        text_results1 <- text_table[grepl(patterns,tolower(text_table$V4)), ]
        text_results2 <- text_results1[as.numeric(text_results1$V6)>2.0, ]
        if(nrow(text_results2) > 0){
            text_results2['Type']='Text mining'
            colnames(text_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','ZScore','ConfidentScore','URL','Type')
            summarytable1 <- text_results2[,c(ncol(text_results2),1:(ncol(text_results2)-1))]
            suppressMessages(write.csv(summarytable1, file="DISEASESTextSummary.csv",row.names = FALSE))
            text_gene_list <- na.omit(summarytable1$GeneSymbol)
        } else { text_gene_list <- list() }
    } else {
        print ("Error: Unable to download file (DISEASES_human_disease_textmining_filtered.tsv).")
        text_gene_list <- list()
    }
    print ("Processing DISEASES (Knowledge category)...")
    url_knowledge <- 'http://download.jensenlab.org/human_disease_knowledge_filtered.tsv'
    suppressMessages(download.file(url_knowledge, destfile='DISEASES_human_disease_knowledge_filtered.tsv'))
    if (file.exists('DISEASES_human_disease_knowledge_filtered.tsv') && file.size('DISEASES_human_disease_knowledge_filtered.tsv') > 0) {
        knowledge_table <- read.csv(file = 'DISEASES_human_disease_knowledge_filtered.tsv', sep = '\t',header=FALSE)
        knowledge_results1 <- knowledge_table[grepl(patterns,tolower(knowledge_table$V4)), ]
        knowledge_results2 <- knowledge_results1[as.numeric(knowledge_results1$V7)>2, ]
        if(nrow(knowledge_results2) > 0){
            knowledge_results2['Type']='Knowledge'
            colnames(knowledge_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','SourceDB','EvidenceType','ConfidenceScore','Type')
            summarytable2 <- knowledge_results2[,c(ncol(knowledge_results2),1:(ncol(knowledge_results2)-1))]
            suppressMessages(write.csv(summarytable2, file="DISEASESKnowledgeSummary.csv",,row.names = FALSE))
            knowledge_gene_list <- na.omit(summarytable2$GeneSymbol)
        } else { knowledge_gene_list <- list() }
    } else {
        print ("Error: Unable to download file (DISEASES_human_disease_knowledge_filtered.tsv).")
        knowledge_gene_list <- list()
    }
    print ("Processing DISEASES (Experiments category)...")
    url_experiments <- 'http://download.jensenlab.org/human_disease_experiments_filtered.tsv'
    suppressMessages(download.file(url_experiments, destfile='DISEASES_human_disease_experiments_filtered.tsv'))
    if (file.exists('DISEASES_human_disease_experiments_filtered.tsv') && file.size('DISEASES_human_disease_experiments_filtered.tsv') > 0) {
        exps_table <- read.csv(file = 'DISEASES_human_disease_experiments_filtered.tsv', sep = '\t',header=FALSE)
        exps_results1 <- exps_table[grepl(patterns,tolower(exps_table$V4)), ]
        exps_results2 <- exps_results1[as.numeric(exps_results1$V7)>2, ]
        if(nrow(exps_results2) > 0){
            exps_results2['Type']='Experiments'
            colnames(exps_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','SourceDB','RankScore','ConfidenceScore','Type')
            summarytable3 <- exps_results2[,c(ncol(exps_results2),1:(ncol(exps_results2)-1))]
            suppressMessages(write.csv(summarytable3, file="DISEASESExpSummary.csv",,row.names = FALSE))
            exps_gene_list <- na.omit(summarytable3$GeneSymbol)
        } else { exps_gene_list <- list() }
    } else {
        print ("Error: Unable to download file (DISEASES_human_disease_experiments_filtered.tsv).")
        exps_gene_list <- list()
    }
    gene_diseases_list <- unique(unlist(c(text_gene_list,knowledge_gene_list,exps_gene_list),recursive=TRUE))
    if (length(gene_diseases_list)>0){
        PCGenesx <- get_PCGenes(gene_diseases_list, gene_file,pcgene_file,"SYMBOL")
    } else {
      print ("No results from DISEASES")
      PCGenesx <- list() }
    return (PCGenesx)
}

get_meshterms <- function(queryTerms,type) {
    #type == "term","id"
    suppressPackageStartupMessages(library(rentrez))
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
               print ("Processing (MESH)...")
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
            }, error=function(e) {
             print ("Error: Please check the output files.")
             MESHList <- list() }) }
            MESHList <- getData()
   } else {
       print ("MESH terms cannot be identified.")
       MESHList <- list()
    }
    if (type=="id"){return(MESHList[[1]])}
    if (type=="term"){return(MESHList[[2]])}
}
get_disgenet <- function(queryTerms,userEmail,userPassword) {
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(disgenet2r))
    gene_file <- "DisGeNET_Genes.txt"
    pcgene_file <- "DisGeNET_PCGenes.txt"
    MeSHTerms <- get_meshterms(queryTerms,"id")
    if (length(MeSHTerms)>0){
        get_disgenetx <- function(){
            tryCatch({
                print ('Processing (DisGeNET)..')
                genelist <- list()
                disgenet_api_key <- get_disgenet_api_key(email = userEmail,password = userPassword )
                Sys.setenv(DISGENET_API_KEY= disgenet_api_key)
                results <- disease2gene( disease = MeSHTerms, vocabulary = "MESH", database = "CURATED", score = c( 0.4,1 ))
                summarytable <- extract(results)
                write.csv(summarytable, file="DisGeNETSummary.csv",row.names = FALSE)
                dbList <- as.character(unique(summarytable$geneid))
                if (length(dbList)>0){
                   MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID")
                } else {
                  print ("No results from DisGeNET.")
                  MappedGenesx <- list() }
            }, error=function(e)
              {
                print ("Error: Unable to get data from DisGeNET.")
                MappedGenesx <- list()
              }
            )
        }
        MappedGenesx <- get_disgenetx()
    } else { MappedGenesx <- list() }
    return(MappedGenesx)
}

get_gtr <- function(queryTerms) {
   suppressPackageStartupMessages(library(rentrez))
   suppressPackageStartupMessages(library(stringr))
   suppressPackageStartupMessages(library(rvest))
   suppressPackageStartupMessages(library(org.Hs.eg.db))
   gene_file <- "GTR_Genes.txt"
   pcgene_file <- "GTR_PCGenes.txt"
   new_query_terms <- character()
   for (var in queryTerms){
       query_term <- paste('"',var,'"[DISNAME]',sep='')
       new_query_terms <- c(new_query_terms, query_term)
   }
   qTerms <- paste('(',paste(new_query_terms,collapse=' OR '),')',sep="")
    dbSearch <- suppressMessages(entrez_search(db='gtr',term=qTerms,use_history=TRUE,retmax=10000))
   if (length(dbSearch$ids)>0){
       getData <- function(){
           tryCatch({
               print ("Processing (GTR)...")
               cuilist <- list()
               if (length(dbSearch$ids)>250){
                    chunks <- split(dbSearch$ids, ceiling(seq_along(dbSearch$ids)/250))
                    #loop over
                    for (i in 1:length(chunks)){
                      summary <- suppressMessages(entrez_summary(db="gtr",id=unlist(chunks[1],recursive=TRUE)))
                      esummaries <- suppressMessages(extract_from_esummary(summary, "conditionlist"))
                      for (i in 1:length(esummaries[1,])){
                          dfcond <- data.frame(cbind(esummaries[,i]$name, esummaries[,i]$cui))
                          dfcondcui <- dplyr::filter(dfcond, grepl(paste(queryTerms,collapse="|"),tolower(X1)))
                          cuilist <- c(cuilist,dfcondcui$X2)
                      }
                    }
               } else {
                    summary <- suppressMessages(entrez_summary(db="gtr",id=unlist(dbSearch$ids,recursive=TRUE)))
                    esummaries <- suppressMessages(extract_from_esummary(summary, "conditionlist"))
                    if (length(unlist(dbSearch$ids,recursive=TRUE))>1){
                      for (i in 1:length(esummaries[1,])){
                          dfcond <- data.frame(cbind(esummaries[,i]$name, esummaries[,i]$cui))
                          dfcondcui <- dplyr::filter(dfcond, grepl(paste(queryTerms,collapse="|"),tolower(X1)))
                          cuilist <- c(cuilist,dfcondcui$X2)
                      }
                    } else {
                        dfcond <- data.frame(cbind(esummaries$name, esummaries$cui))
                        dfcondcui <- dplyr::filter(dfcond, grepl(paste(queryTerms,collapse="|"),tolower(X1)))
                        cuilist <- c(cuilist,dfcondcui$X2)
                    }
               }
               cuilist <- unique(unlist(cuilist,recursive=TRUE))
               if (length(cuilist)>0){
                  genelist <- list()
                  summarylist <- list()
                  for (cuid in cuilist){
                    print (cuid)
                    url <- gsub(' ','',paste('https://www.ncbi.nlm.nih.gov/gtr/conditions/',cuid,collapse=''))
                    page_html <- gsub('\n','',toString(read_html(url)))
                    gtrmatches <- str_match(page_html,'<ul class="associate_genes gtr-reset-list">.*</ul>')[[1]]
                    gtrmatches2 <- strsplit(strsplit(gtrmatches,"</ul>")[[1]][1], '<li>')[[1]][-1]
                    for (g in gtrmatches2){
                      geneid <- gsub('/gtr/genes/','',str_match_all(g,'/gtr/genes/\\d+')[1])
                      genesymbol <- gsub("<.*?>", "", strsplit(g,"</a>")[[1]][1])
                      gtemp <- gsub("<.*?>", "", strsplit(g,"</a>")[[1]][3])
                      genename <- strsplit(gtemp,"Summary: ")[[1]][2]
                      synonym <- gsub("Also known as: ","",strsplit(gtemp,"Summary: ")[[1]][1])
                      namex <- paste(cuid,"_",geneid,sep="")
                      df <- data.frame(cuid, geneid, genesymbol, genename, synonym)
                      summarylist[[namex]] <- df
                      genelist <- c(genelist,geneid)
                    }
                  }
                  summarytable <- do.call(rbind,summarylist)
                  colnames(summarytable) <- c("CUID","GeneID","GeneSymbol","GeneName","Synonym")
                  write.csv(summarytable, file="GTRSummary.csv",row.names = FALSE)
                  genelist <- unique(unlist(genelist,recursive=TRUE))
                  if (length(genelist)>0){
                    MappedGenesx <- get_PCGenes(genelist, gene_file,pcgene_file,"ENTREZID")
                  } else {
                    print ("No results from GTR.")
                    MappedGenesx <- list() }
               } else {
                 print ("No results from GTR.")
                 MappedGenesx <- list()
               }
            }, error=function(e)
                {
                    print("Error: Please check the output files.")
                    MappedGenesx <- list()
                }) }
       MappedGenesx <- getData()
   } else {
       print ("No results from GTR.")
       MappedGenesx <- list()
    }
   return(MappedGenesx)
}

get_gwascatalog <- function(queryTerms) {
    suppressPackageStartupMessages(library(gwasrapidd))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    gene_file <- "GWASCatalog_Genes.txt"
    pcgene_file <- "GWASCatalog_PCGenes.txt"
    gwas_list = character()
    summarylist <- list()
    print ("Processing (GWAS Catalog)... ")
    for (qt in queryTerms){
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
                        tryCatch({ var6 <- get_traits(association_id = my_associations2[125,]@entrez_ids$association_id)@traits$efo_id }, error=function(e) { var6 <- "NA" })
                    }
                    getvar7 <- function(){
                        tryCatch({ var7 <- get_traits(association_id = my_associations2[i,]@entrez_ids$association_id)@traits$trait }, error=function(e) { var6 <- "NA" })
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
        write.csv(summarytable, file="GWASCatalogSummary.csv",row.names = FALSE)
        gwas_list2 <- na.omit(gwas_list)
        gwas_list2 <- unique(gwas_list2)
        if (length(gwas_list2)>0){
            gwas_protein_coding_genesx <- get_PCGenes(gwas_list2, gene_file,pcgene_file,"ENTREZID")
        } else {
          print ("No results from GWAS Catalog.")
          gwas_protein_coding_genesx <- list() }
    } else {
      gwas_protein_coding_genesx <- list()
      print ("No results from GWAS Catalog.")}
    return(gwas_protein_coding_genesx)
}

get_medgen <- function(queryTerms) {
   suppressPackageStartupMessages(library(rentrez))
   suppressPackageStartupMessages(library(stringr))
   suppressPackageStartupMessages(library(org.Hs.eg.db))
   gene_file <- "MedGen_Genes.txt"
   pcgene_file <- "MedGen_PCGenes.txt"
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
               print ("Processing (MedGen)...")
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
                        dbList <- strsplit(dbFetch,"\n")
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
                        dbList <- strsplit(dbFetch,"\n")
                        dbList <- unlist(dbList,recursive=TRUE)
                        dbList <- head(dbList,-1)
                    } else { dbList <- list() }
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("UID","SemanticID","SemanticType","ConceptID","Title","Definition")
                write.csv(summarytable, file="MedGenSummary.csv",row.names = FALSE)
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){ MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID") }
               else {
                 print ("No results from MedGen.")
                 MappedGenesx <- list() }
            }, error=function(e) {
             print ("Error: Please check the output files.")
             MappedGenesx <- list() })
            }
            MappedGenesx <- getData()
   } else {
     print ("No results from MedGen.")
     MappedGenesx <- list()
    }
   return(MappedGenesx)
}
get_omim <- function(queryTerms) {
   suppressPackageStartupMessages(library(rentrez))
   suppressPackageStartupMessages(library(stringr))
   suppressPackageStartupMessages(library(org.Hs.eg.db))
   gene_file <- "OMIM_Genes.txt"
   pcgene_file <- "OMIM_PCGenes.txt"
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
               print ("Processing (OMIM)...")
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
                      dbList <- strsplit(dbFetch,"\n")
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
                      dbList <- strsplit(dbFetch,"\n")
                      dbList <- unlist(dbList,recursive=TRUE)
                      dbList <- head(dbList,-1)
                    } else {dbList <- list()}
                }
                summarytable <- do.call(rbind,summarylist)
                colnames(summarytable) <- c("Accession","Title","AltTitles")
                write.csv(summarytable, file="OMIMSummary.csv",row.names = FALSE)
                #Get Protein Coding Gene IDs
               if (length(dbList)>0){ MappedGenesx <- get_PCGenes(dbList, gene_file,pcgene_file,"ENTREZID") }
               else {
                 print ("No results from OMIM.")
                 MappedGenesx <- list() }
            }, error=function(e) {
             print ("Error: Please check the output files.")
             MappedGenesx <- list() })
       }
            MappedGenesx <- getData()
   } else { MappedGenesx <- list() }
   return(MappedGenesx)
}
get_phegeni <- function(queryTerms) {
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    gene_file <- "PheGenI_Genes.txt"
    pcgene_file <- "PheGenI_PCGenes.txt"
    url_download <- 'http://www.ncbi.nlm.nih.gov/projects/gap/eqtl/EpiViewBE.cgi?type=dl.tab'
    suppressMessages(download.file(url_download, destfile='phegeni.tab'))
    if (file.exists('phegeni.tab') && file.size('phegeni.tab') > 0) {
        print ("Processing (PheGenI)... ")
        phegeni <- readLines(file('phegeni.tab'))
        phegeni_tab <- strsplit(phegeni,'\t')
        phegeni_tab_short <- lapply(phegeni_tab, function(x) {x[c(1:13)]})
        patterns <- tolower(paste(queryTerms, collapse = "|"))
        summarylist <- Filter(function(x) grepl(patterns, tolower(x[2])), phegeni_tab_short)
        summarytable <- do.call(rbind,summarylist)
        if (length(summarytable[,1])>0){
            colnames(summarytable) <- c("ID","Trait", "SNP rs", "Context", "Gene", "Gene ID", "Gene 2", "Gene ID 2", "Chromosome", "Location", "P-Value", "Source", "PubMed")
            write.csv(summarytable, file="PheGenISummary.csv",row.names = FALSE)
            genelist <- unique(unlist(c(summarytable[,6],summarytable[,8]),recursive=TRUE))
            if (length(genelist)>0){
              MappedGenesx <- get_PCGenes(genelist, gene_file,pcgene_file,"ENTREZID")
            } else {
              print ("No results from PheGenI.")
              MappedGenesx <- list() }
        } else {
          print ("No results from PheGenI.")
          MappedGenesx <- list() }
    } else {
        print ("Error: Unable to download file (ncbi/gap/phegeni).")
        MappedGenesx <- list()
    }
    return(MappedGenesx)
}

GeneSearchDisGeNETYes <- function(queryTerms,EmailID, PasswordID){
  pack_cran = c("rentrez", "stringr", "rvest","devtools")
  pack_bioc = c("AnnotationDbi","org.Hs.eg.db")
  pack_github = c("disgenet2r","gwasrapidd")
  suppressMessages(check_libraries(pack_cran, pack_bioc, pack_github))
  queryTerms <- unlist(strsplit(queryTerms,"\\|"),recursive=TRUE)
  ClinVarGenes <- get_clinvar(queryTerms)
  DISEASESGenes <- get_diseasesdb(queryTerms)
  DisGeNETGenes <- get_disgenet(queryTerms,EmailID, PasswordID)
  GTRGenes <- get_gtr(queryTerms)
  GWASGenes <- get_gwascatalog(queryTerms)
  MedGenGenes <- get_medgen(queryTerms)
  OMIMGenes <- get_omim(queryTerms)
  PheGenIGenes <- get_phegeni(queryTerms)
  print ("Preparing a summary...")
  GeneCombined <- c(ClinVarGenes, DISEASESGenes, DisGeNETGenes, GTRGenes, GWASGenes,MedGenGenes,OMIMGenes,PheGenIGenes)
  GeneCombined <- paste(GeneCombined, collapse="\n")
  write(GeneCombined,file="CombinedGenes.txt")
  pdf('GeneSearch.pdf',paper="a4r")
  colors <- c("darkblue","red","green","orange","cyan","magenta","black","yellow")
  dtb <- c("ClinVar", "DISEASES", "DisGeNET","GTR","GWAS","MedGen","OMIM","PheGenI")
  nogenes <- c(length(ClinVarGenes),length(DISEASESGenes),length(DisGeNETGenes),length(GTRGenes),length(GWASGenes),length(MedGenGenes),length(OMIMGenes),length(PheGenIGenes))
    barp <- barplot(nogenes,
        main = "Protein-coding genes from database search",
        ylab = "No. of genes",
        border = "black",
        col = colors,
        names=dtb,
        las=3 )
    text(barp, nogenes + 6.0, labels = nogenes)
    dev.off()

}
GeneSearchDisGeNETNo <- function(queryTerms){
  pack_cran = c("rentrez", "stringr", "rvest","devtools")
  pack_bioc = c("AnnotationDbi","org.Hs.eg.db")
  pack_github = c("gwasrapidd")
  suppressMessages(check_libraries(pack_cran, pack_bioc, pack_github))
  queryTerms <- unlist(strsplit(queryTerms,"\\|"),recursive=TRUE)
  ClinVarGenes <- get_clinvar(queryTerms)
  DISEASESGenes <- get_diseasesdb(queryTerms)
  GTRGenes <- get_gtr(queryTerms)
  GWASGenes <- get_gwascatalog(queryTerms)
  MedGenGenes <- get_medgen(queryTerms)
  OMIMGenes <- get_omim(queryTerms)
  PheGenIGenes <- get_phegeni(queryTerms)
  print ("Preparing a summary...")
  GeneCombined <- c(ClinVarGenes, DISEASESGenes, GTRGenes, GWASGenes,MedGenGenes,OMIMGenes,PheGenIGenes)
  GeneCombined <- paste(GeneCombined, collapse="\n")
  write(GeneCombined,file="CombinedGenes.txt")
  pdf('GeneSearch.pdf',paper="a4r")
  colors <- c("darkblue","red","orange","cyan","magenta","black","yellow")
  dtb <- c("ClinVar", "DISEASES","GTR","GWAS","MedGen","OMIM","PheGenI")
  nogenes <- c(length(ClinVarGenes),length(DISEASESGenes),length(GTRGenes),length(GWASGenes),length(MedGenGenes),length(OMIMGenes),length(PheGenIGenes))
    barp <- barplot(nogenes,
        main = "Protein-coding genes from database search",
        ylab = "No. of genes",
        border = "black",
        col = colors,
        names=dtb,
        las=3 )
    text(barp, nogenes + 6.0, labels = nogenes)
    dev.off()

}

cat("Please enter a disease or keywords separated by '|' (e.g. colorectal cancer|colorectal carcinoma): ");
query_terms <- readLines("stdin",n=1);
cat( "\n" )

cat("Do you want to include DisGeNET in gene searches (email and password to DisGeNET account must be provided) (yes/no): ");
yes_no <- readLines("stdin",n=1);
if (yes_no != "yes" && yes_no != "no"){
    cat("Please enter 'yes' OR 'no'")
    cat( "\n" )
}
if (yes_no =="yes"){
    cat("You've entered 'yes'.Please provide email and password to DisGENET account:")
    cat( "\n" )
    cat("Email address (DisGENET account): ");
    emailx <- readLines("stdin",n=1);
    cat("Password (DisGENET account): ");
    passwordx <- readLines("stdin",n=1);
    cat( "\n" )
    dir.create("GeneSearch")
    setwd("GeneSearch")
    GeneSearchDisGeNETYes(query_terms,emailx,passwordx)
    setwd("../")
    #print ("Run enrichment analysis.")
}
if (yes_no =="no"){
    cat("You've entered 'no'. DisGeNET will not be included in gene searches.")
    cat( "\n" )
    dir.create("GeneSearch")
    setwd("GeneSearch")
    GeneSearchDisGeNETNo(query_terms)
    setwd("../")
}



#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#get_clinvar(queryTerms)
