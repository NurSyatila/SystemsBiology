#============= Gene_Search: ================
#Function: Get genes for a given query disease in GTR
#Usage: get_GTR(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
# - summaryFile: summary file to record output
#Value: List of disease-associated genes

get_gtr <- function(queryTerms) {
    library(stringr)
    library(rvest)
    library(org.Hs.eg.db)
    gene_file <- "Genes_GTR.txt"
    pcgene_file <- "PCGenes_GTR.txt"
    new_queryTerms <- character()
    for (var in queryTerms){
        query_term <- paste('"',var,'"[DISNAME]',sep='')
        query_term <- gsub(' ','+',query_term)
        new_queryTerms <- c(new_queryTerms, query_term)
    }
    new_queryTerms <- paste(new_queryTerms, collapse='+OR+')
    queryTerms_gtr <- paste('(',new_queryTerms,')',sep='')
    print ("Processing...")
    url <- paste('https://www.ncbi.nlm.nih.gov/gtr/all/genes/?term=',queryTerms_gtr,collapse='')
    url <- gsub(' ','',url)
    #load the starting page
    firstpage <- html_text(read_html(url))
    nopage <- str_extract(firstpage,'1 of \\d+')
    if (is.na(nopage)==FALSE){
        summarylist <- list()
        genelist <- list()
        nopage <- str_split(nopage,' of ')
        nopage <- unlist(nopage,recursive=TRUE)
        lastpage <- strtoi(tail(nopage,n=1))
        #get ids in every page
        for (var in 1:lastpage){
            newurl <- paste(url,'&page=',var,sep='')
            newurl <- gsub(' ','',newurl)
            page_read <- read_html(newurl)
            page_html <- gsub('\n','',toString(page_read))
            gtrmatches <- str_match_all(page_html,'<table class="jig-ncbigrid rprts".*></table>')[[1]]
            gtrmatches2 <- head(strsplit(gtrmatches,'</tr>')[[1]][-1],-1)
            for (g in gtrmatches2){
                gx <- strsplit(g,'<td>')[[1]][-1]
                geneid <- gsub('/gtr/genes/','',str_match_all(gx[2],'/gtr/genes/\\d+')[1])
                genesymbol <- gsub('ref=\"link_id=gene_symbol\">','',str_match_all(gx[2],'ref=\"link_id=gene_symbol\">[A-Za-z0-9]+')[1])
                genenamex <- str_match_all(gx[2],"destText:.*ref=\"link_id=gene_symbol\">")
                genename <- gsub("'\" ref=\"link_id=gene_symbol\">", '', gsub("destText: '",'',genenamex))
                condlist <- strsplit(strsplit(gx[3],"ref=\"link_id=associated_cond\">")[[1]][-1],'</a>')
                condlist <- paste(sapply(condlist, function(x) x[1]),collapse='; ')
                df <- data.frame(geneid, genesymbol, genename, condlist)
                summarylist[[geneid]] <- df
                genelist <- c(genelist,geneid)
            }
        }
        summarytable <- do.call(rbind,summarylist)
        colnames(summarytable) <- c("GeneID","GeneSymbol","GeneName","Conditions (ID)")
        write.csv(summarytable, file="SummaryGTR.csv",row.names = FALSE)
        print ("Get summary of query searches, output file: SummaryGTR.csv ")
        #Get Protein Coding Gene IDs
        genelist <- unlist(genelist,recursive=TRUE)
        if (length(genelist)>0){
            Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = genelist, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "ENTREZID")
            write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
            print ("Get disease-related genes, output file: Genes_GTR.txt ")
            PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
            PCGenes <- unique(PCGenes[,1])
            if (length(PCGenes)>0){
                MapGenes <- paste(PCGenes, collapse="\n")
                write(MapGenes,file=pcgene_file)
                print ("Get protein-coding genes, output file: PCGenes_GTR.txt ")
                MappedGenesx <- PCGenes
            }
            else {
                print ("No disease-related protein-coding genes from GT.")
                MappedGenesx <- list()
            }
        } else {
            print ("No disease-related genes from GTR.")
            MappedGenesx <- list()
        }

    } else {
        print ("No results from GTR.")
        MappedGenesx <- list()
    }
    return(MappedGenesx)
}


#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#get_gtr(queryTerms)
