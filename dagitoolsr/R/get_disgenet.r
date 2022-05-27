#============= Gene_Search: ================
#Function: Get genes for a given query disease in DisGeNET database using DisGeNET2r package
#Usage: get_DisGeNET(queryTerms,summaryFile)
#Arguments:
# - MeSHTerms: MeSH terms in list
# - userEmail: DisGeNET account ID (nsag@ukm.edu.my)
# - userPassword: DisGeNET account password (ns65ns65)
# - summaryFile: summary file to record output
#Value: List of disease-associated genes

get_disgenet <- function(MeSHTerms,userEmail,userPassword) {
    library(stringr)
    library(org.Hs.eg.db)
    library(disgenet2r)
    gene_file <- "Genes_DisGeNET.txt"
    pcgene_file <- "PCGenes_DisGeNET.txt"
    if (length(MeSHTerms)>0){
        get_disgenetx <- function(){
            tryCatch({
                print ('Processing..')
                genelist <- list()
                disgenet_api_key <- get_disgenet_api_key(email = userEmail,password = userPassword )
                Sys.setenv(DISGENET_API_KEY= disgenet_api_key)
                results <- disease2gene( disease = MeSHTerms, vocabulary = "MESH", database = "CURATED", score = c( 0.4,1 ))
                summarytable <- extract(results)
                write.csv(summarytable, file="SummaryDisGeNET.csv",row.names = FALSE)
                print ("Get summary of query searches, output file: SummaryDisGeNET.csv ")
                dbList <- as.character(unique(summarytable$geneid))
                if (length(dbList)>0){
                   Genes <- AnnotationDbi::select(org.Hs.eg.db, keys = dbList, c("ENTREZID", "SYMBOL","GENETYPE"), keytype = "ENTREZID")
                   write.table(Genes,file=gene_file, quote=FALSE,sep="\t",row.names = FALSE)
                   print ("Get disease-related genes, output file: Genes_DisGeNET.txt ")
                   PCGenes <- dplyr::filter(Genes, grepl("protein",tolower(GENETYPE)))
                   PCGenes <- unique(PCGenes[,1])
                   if (length(PCGenes)>0){
                       MapGenes <- paste(PCGenes, collapse="\n")
                       write(MapGenes,file=pcgene_file)
                       print ("Get protein-coding genes, output file: PCGenes_DisGeNET.txt ")
                       MappedGenesx <- PCGenes
                   } else {
                       print ("No disease-related protein-coding genes from DisGeNET.")
                       MappedGenesx <- list()
                   }
                } else {
                    print ("No disease-related genes from DisGeNET.")
                    MappedGenesx <- list()
                }
            }, error=function(e)
                {
                    print ("Error: Unable to get data from DisGeNET.")
                    MappedGenesx <- list()
                }
            )
        }
        MappedGenesx <- get_disgenetx()
    } else {
        print ("No results from DisGeNET.")
        MappedGenesx <- list()
    }
    return(MappedGenesx)
}

#queryTerms <- c("endometriosis")
#queryTerms <- c("colorectal cancer")
#MeSHTerms <- get_meshterms(queryTerms,"id")
#userEmail <- 'nsag@ukm.edu.my'
#userPassword <- 'ns65ns65'
#get_disgenet(MeSHTerms,userEmail,userPassword)


#D004715
#68005831
