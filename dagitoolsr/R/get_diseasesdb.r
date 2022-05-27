#============= Gene_Search: ================
#Function: Get genes for a given query disease in DISEASES database
#Usage: get_DISEASESdb(queryTerms,summaryFile)
#Arguments:
# - queryTerms: Query terms in list
#Value: List of disease-associated genes

get_diseasesdb <- function(queryTerms) {
    gene_file <- "Genes_DISEASES.txt"
    pcgene_file <- "PCGenes_DISEASES.txt"
    #options(download.file.method="curl",timeout = 300)
    #options(download.file.method="curl")
    patterns <- paste(queryTerms, collapse = "|")
    print ("Processing Text Mining category...")
    url_text_mining <- 'https://download.jensenlab.org/human_disease_textmining_filtered.tsv'
    suppressMessages(download.file(url_text_mining, destfile='human_disease_textmining_filtered.tsv'))
    if (file.exists('human_disease_textmining_filtered.tsv') && file.size('human_disease_textmining_filtered.tsv') > 0) {
        text_table <- read.csv(file = 'human_disease_textmining_filtered.tsv', sep = '\t',header=FALSE)
        text_results1 <- text_table[grepl(patterns,tolower(text_table$V4)), ]
        text_results2 <- text_results1[as.numeric(text_results1$V6)>2.0, ]
        if(nrow(text_results2) > 0){
            text_results2['Type']='Text mining'
            colnames(text_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','ZScore','ConfidentScore','URL','Type')
            summarytable1 <- text_results2[,c(ncol(text_results2),1:(ncol(text_results2)-1))]
            print ("Get summary of query searches in Text Mining category, output file: SummaryDISEASES.csv ")
            suppressMessages(write.table(summarytable1, file="SummaryDISEASES.csv",append=TRUE, quote=FALSE, sep=",",row.names = FALSE))
            text_gene_list <- na.omit(summarytable1$GeneSymbol)
        } else {
            print ("No results from Text Mining category.")
            text_gene_list <- list()
        }

    } else {
        print ("Error: Unable to download file (human_disease_textmining_filtered.tsv).")
        text_gene_list <- list()
    }
    print ("Processing Knowledge category...")
    url_knowledge <- 'https://download.jensenlab.org/human_disease_knowledge_filtered.tsv'
    suppressMessages(download.file(url_knowledge, destfile='human_disease_knowledge_filtered.tsv'))
    if (file.exists('human_disease_knowledge_filtered.tsv') && file.size('human_disease_knowledge_filtered.tsv') > 0) {
        knowledge_table <- read.csv(file = 'human_disease_knowledge_filtered.tsv', sep = '\t',header=FALSE)
        knowledge_results1 <- knowledge_table[grepl(patterns,tolower(knowledge_table$V4)), ]
        knowledge_results2 <- knowledge_results1[as.numeric(knowledge_results1$V7)>2, ]
        if(nrow(knowledge_results2) > 0){
            knowledge_results2['Type']='Knowledge'
            colnames(knowledge_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','SourceDB','EvidenceType','ConfidenceScore','Type')
            summarytable2 <- knowledge_results2[,c(ncol(knowledge_results2),1:(ncol(knowledge_results2)-1))]
            print ("Get summary of query searches in Knowledge category, output file: SummaryDISEASES.csv ")
            suppressMessages(write.table(summarytable2, file="SummaryDISEASES.csv",append=TRUE, quote=FALSE, sep=",",row.names = FALSE))
            knowledge_gene_list <- na.omit(summarytable2$GeneSymbol)
        } else {
            print ("No results from Knowledge category.")
            knowledge_gene_list <- list()
        }
    } else {
        print ("Error: Unable to download file (human_disease_knowledge_filtered.tsv).")
        knowledge_gene_list <- list()
    }
    print ("Processing Experiments category...")
    url_experiments <- 'https://download.jensenlab.org/human_disease_experiments_filtered.tsv'
    suppressMessages(download.file(url_experiments, destfile='human_disease_experiments_filtered.tsv'))
    if (file.exists('human_disease_experiments_filtered.tsv') && file.size('human_disease_experiments_filtered.tsv') > 0) {
        exps_table <- read.csv(file = 'human_disease_experiments_filtered.tsv', sep = '\t',header=FALSE)
        exps_results1 <- exps_table[grepl(patterns,tolower(exps_table$V4)), ]
        exps_results2 <- exps_results1[as.numeric(exps_results1$V7)>2, ]
        if(nrow(exps_results2) > 0){
            exps_results2['Type']='Experiments'
            colnames(exps_results2) <- c('EnsemblID','GeneSymbol','DO','Disease','SourceDB','RankScore','ConfidenceScore','Type')
            summarytable3 <- exps_results2[,c(ncol(exps_results2),1:(ncol(exps_results2)-1))]
            print ("Get summary of query searches in Experiments category, output file: SummaryDISEASES.csv ")
            suppressMessages(write.table(summarytable3, file="SummaryDISEASES.csv",append=TRUE, quote=FALSE, sep=",",row.names = FALSE))
            exps_gene_list <- na.omit(summarytable3$GeneSymbol)
        } else {
            print ("No results from Experiments category.")
            exps_gene_list <- list()
        }
    } else {
        print ("Error: Unable to download file (human_disease_experiments_filtered.tsv).")
        exps_gene_list <- list()
    }
    gene_diseases_list <- unique(unlist(c(text_gene_list,knowledge_gene_list,exps_gene_list),recursive=TRUE))
    if (length(gene_diseases_list)>0){
        Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_diseases_list, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "SYMBOL")
        suppressMessages(write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE))
        print ("Get disease-related genes, output file: Genes_DISEASES.txt ")
        PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
        PCGenes <- unique(PCGenes[,2])
        if (length(PCGenes)>0){
            to_write_diseases <- paste(PCGenes, collapse="\n")
            write(to_write_diseases,file=pcgene_file)
            PCGenesx <- PCGenes
            print ("Get protein-coding genes, output file: PCGenes_DISEASES.txt ")
        }
    } else {
        print ("No results from DISEASE db (Jensen).")
        PCGenesx <- list()
    }
    return (PCGenesx)
}

#queryTerms <- c("breast cancer")
#queryTerms <- c("colorectal cancer")
#patterns <- paste(queryTerms, collapse = "|")
#get_diseasesdb(queryTerms)
#queryTerms <- c("glaucoma")