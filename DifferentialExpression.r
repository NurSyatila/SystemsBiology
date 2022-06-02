#============= Differential_Expression ====================
#Function:
# 1) Get GEO datasets for given MeSH terms
# 2) For individual GEO dataset, perform differential expression analysis
# 3) Retrieve differentially expressed genes from individual GSE data set
check_libraries <- function(pack_cran, pack_bioc) {
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
    print("Loading required libraries, this might takes a while..")
    all_packages <- c(pack_cran,pack_bioc)
    for (pack in all_packages){
        suppressPackageStartupMessages({library(pack, character.only = TRUE)})
    }
}

get_meshterms <- function(queryTerms,type) {
    #type == "term","id"
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

get_gse_datasets <- function(queryTerms) {
    MeSHTerms <- get_meshterms(queryTerms,"term")
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
               print ("Processing (GEO DataSets)...")
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
                write.csv(summarytable, file="GSESummary.csv",row.names = FALSE)
                Filteredsummarytable <- dplyr::filter(summarytable, grepl("disease",tolower(summarytable$subsetinfo)))
                Filteredsummarytable <- Filteredsummarytable[!duplicated(Filteredsummarytable[,"gse"]),]
            }, error=function(e)
                {
                    print("Please check the output files.")
                    #GSElist <- list()
                    Filteredsummarytable <- data.frame(matrix(ncol = 5, nrow = 0))
                }
                )
            }
            #GSElist <- getData()
            Filteredsummarytable <- getData()
   } else {
       print ("No results from GEO Datasets.")
       #GSElist <- list()
       Filteredsummarytable <- data.frame(matrix(ncol = 5, nrow = 0))
    }
    return (Filteredsummarytable)
}

analyse_deg <- function(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile){
    gene_annot <- suppressMessages(fData(gset))
    gene_ids <- unlist(strsplit(gene_annot$"Gene ID","///"),recursive=TRUE)
    pc_genes <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = gene_ids, c("ENTREZID", "GENETYPE"), keytype = "ENTREZID"))
    pc_genes <- dplyr::filter(pc_genes, grepl("protein",tolower(GENETYPE)))
    pc_genes <- unique(pc_genes$ENTREZID)
    gene_annot2 <- gene_annot[gene_annot$"Gene ID" %in% pc_genes,]
    gset2 <- gset[rownames(gset) %in% gene_annot2$"ID",]
    # log2 transformation
    ex <- exprs(gset2)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
        exprs(gset2) <- log2(ex) }

    # assign samples to groups and set up design matrix
    gs <- factor(sampleGrouping)
    groups <- make.names(c("groupA","groupB"))
    levels(gs) <- groups
    gset2$group <- gs
    design <- model.matrix(~group + 0, gset2)
    colnames(design) <- levels(gs)
    fit <- lmFit(gset2, design)  # fit linear model
    # set up contrasts of interest and recalculate model coefficients
    cts <- paste(groups[1], groups[2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    full_results <- topTable(fit2, adjust="BH", sort.by="p", number=Inf)
    # save DEG results
    new_results <- filter(full_results, adj.P.Val < 0.05, abs(logFC) > 2.0)
    if (length(new_results)>0){
        new_results <- subset(new_results, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.ID","Gene.symbol","Gene.title"))
        new_results <- subset(new_results, !is.na(Gene.symbol))
        new_results <- subset(new_results, Gene.symbol != "")
        deg_list <- unique(na.omit(new_results$Gene.ID))
        if (length(deg_list)>0){
            dir.create(GSEaccession)
            write.csv(new_results,file=paste(GSEaccession,"/","FilteredDE.csv",sep=""))
            MapGenes <- paste(deg_list, collapse="\n")
            write(MapGenes,file=paste(GSEaccession,"/","DEGenes.txt",sep=""))
            # summarize test results as "up", "down" or "not expressed"
            dT<-decideTests(fit2, lfc = 2.0,p.value=0.05,method="separate",adjust.method="BH")
            write("** Summary **",file=summaryFile,append=TRUE)
            write(groupDescription,file=summaryFile,append=TRUE)
            write(knitr::kable(summary(dT)),file=summaryFile,append=TRUE)
            write("\n",file=summaryFile,append=TRUE)
            pdf(paste(GSEaccession,"/","Plots.pdf",sep=""))
            # Venn diagram of results
            vennDiagram(dT, circle.col=palette())
            # create Q-Q plot for t-statistic
            t.good <- which(!is.na(fit2$F)) # filter out bad probes
            qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
            #get volcano plot
            p_cutoff <- 0.05
            fc_cutoff <- 2.0
            ggplot(new_results,aes(x = logFC, y=B)) + geom_point()
            new_results %>%
              mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
            # Visualize and quality control test results.
            # Build histogram of P-values for all genes. Normal test
            # assumption is that most genes are not differentially expressed.
            tT2 <- topTable(fit2, adjust="BH", sort.by="p", number=Inf)
            hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj", ylab = "Number of genes", main = "P-adj value distribution")
            # volcano plot (log P-value vs log fold change)
            colnames(fit2) # list contrast names
            ct <- 1        # choose contrast of interest
            volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20, highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
            # MD plot (log fold change vs mean log expression)
            # highlight statistically significant (p-adj < 0.05) probes
            plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
            abline(h=0)
            ################################################################
            # General expression data analysis
            ex <- exprs(gset2)
            # box-and-whisker plot
            ord <- order(gs)  # order samples by group
            palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02", "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
            par(mar=c(7,4,2,1))
            title <- paste (GSEaccession, "/", annotation(gset2), sep ="")
            boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
            legend("topleft", groups, fill=palette(), bty="n")
            # expression value distribution
            par(mar=c(4,4,2,1))
            title <- paste (GSEaccession, "/", annotation(gset2), " value distribution", sep ="")
            plotDensities(ex, group=gs, main=title, legend ="topright")
            # mean-variance trend, helps to see if precision weights are needed
            plotSA(fit2, main=paste("Mean variance trend, ",GSEaccession,sep=""))
            if (length(deg_list)>20){
                #get top 20 genes
                topN <- 20
                ##
                ids_of_interest <- mutate(new_results, Rank = 1:n()) %>%
                filter(Rank < topN) %>%
                pull(ID)
                gene_names <- mutate(new_results, Rank = 1:n()) %>%
                filter(Rank < topN) %>%
                pull(Gene.symbol)
                ## Get the rows corresponding to ids_of_interest and all columns
                tryCatch({
                    gene_matrix <- exprs(gset2)[ids_of_interest,]
                    pheatmap(gene_matrix, labels_row = gene_names, scale="row")
                    #write(paste(gene_names,collapse="\n"),file=paste(GSEaccession,"/","Top20DEGenes.txt",sep=""))
                    top20genesx <- new_results[new_results$ID %in% ids_of_interest,]
                    top20genes <- subset(top20genesx, select=c("ID","logFC","Gene.ID","Gene.symbol","Gene.title"))
                    top20genes$DE <- ""
                    top20genes$DE[top20genes$DE > 1.9] <- "Up"
                    top20genes$DE[top20genes$DE < - 1.9] <- "Down"
                    suppressMessages(write.table(top20genes, file=paste(GSEaccession,"/","Top20DEGenes.txt",sep=""), quote=FALSE,sep="\t",row.names = FALSE))
                    #write(MapGenes,file=paste(GSEaccession,"/","DEGenes.txt",sep=""))
                }, error=function(e)
                {
                    print("Error: Unable to create heatmap.")
                })

            }
            dev.off()
            degsx_list <- deg_list
        }
        else {
            write("** No protein-coding genes (DEG) identified **\n",file=summaryFile,append=TRUE)
            write(groupDescription,file=summaryFile,append=TRUE)
            degsx_list <- list()
        }
    } else {
        write("** No differentially expressed genes found **\n",file=summaryFile,append=TRUE)
        degsx_list <- list()
    }
    return (degsx_list)
}

get_deg <- function(queryTerms,GSEaccession,GSEplatform,GSETitle) {
    summaryFile <- "GSESummary.txt"
    write(paste("-- ",GSEaccession,": ",GSETitle," --",collapse=" "),file=summaryFile,append=TRUE)
    parse_gsefile <- function(GSEaccession){
        tryCatch({
            gset <- suppressMessages(getGEO(GSEaccession, GSEMatrix =TRUE, AnnotGPL=TRUE,destdir="."))
            if (length(gset) > 1) idx <- grep(GSEplatform, attr(gset, "names")) else idx <- 1
            gset <- gset[[idx]]
            samples <- pData(gset)$characteristics_ch1
            if (!is.null(samples)==TRUE){
              groups <- unique(samples)
              if (length(groups)==2){
                   sample_groupings <- gsub(groups[1],"0",samples)
                   sampleGrouping <- gsub(groups[2],"1",sample_groupings)
                   groupDescription <- paste("--GroupA: ",groups[1],"\n","--GroupB: ",groups[2],sep="")
                   degsx_list <- suppressMessages(analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile))
              } else {
                   patterns <- tolower(paste(queryTerms, collapse = "|"))
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
                       degsx_list <- suppressMessages(analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile))
                   } else {
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
                           degsx_list <- suppressMessages(analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile))
                       } else {
                           print ("** Error: More/less than two sample types detected. The GSE data set will not be used. **")
                           write("** More than two sample types detected. The GSE data set will not be used. **\n--GSE samples:",file=summaryFile,append=TRUE)
                           write(paste(groups,collapse="\n"),file=summaryFile,append=TRUE)
                           degsx_list <- list()
                       }
                   }
              }
            } else {
              print("Error (2): Sample groupings cannot be determined.")
              write("** Error - Sample groupings cannot be determined.",file=summaryFile,append=TRUE)
              degsx_list <- list()
            }

        }, error=function(e)
                {
                    print("Error (1): Unable to download / parse the GSE file.")
                    write("** Error - Unable to download / parse the GSE file.",file=summaryFile,append=TRUE)
                    degsx_list <- list()
                }
        )
    }
    deg_list <- parse_gsefile(GSEaccession)
    return (deg_list)
}
DifferentialExpressionFramework <- function(queryTerms) {
    pack_cran = c("umap","pheatmap","ggplot2","dplyr","rentrez", "stringr","VennDiagram","ggnewscale")
    pack_bioc = c("GEOquery","limma","AnnotationDbi","org.Hs.eg.db")
    check_libraries(pack_cran, pack_bioc)
    queryTerms <- unlist(strsplit(queryTerms,"\\|"),recursive=TRUE)
    GSESum <- get_gse_datasets(queryTerms)
    if (length(GSESum[,1])>0){
     write.csv(GSESum, file="GSESummaryFiltered.csv",row.names = FALSE)
     for (i in 1:length(GSESum[1,])){
       if (!is.na(GSESum[i,]$gse)==TRUE){
          GSEaccession <- GSESum[i,]$gse
          GSEplatform <- GSESum[i,]$gpl
          GSETitle <- GSESum[i,]$seriestitle
          print(GSEaccession)
          suppressMessages(get_deg(queryTerms,GSEaccession,GSEplatform,GSETitle))
       }
    }
  } else {print('No results from GEO DataSets.')}
}

cat("Please enter a disease or keywords separated by '|' (e.g. colorectal cancer|colorectal carcinoma): ");
query_terms <- readLines("stdin",n=1);
cat( "\n" )
dir.create("DifferentialExpression")
setwd("DifferentialExpression")
DifferentialExpressionFramework(query_terms)
setwd("../")
