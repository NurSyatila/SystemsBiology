#============= Differential_Expression ====================
#Function: Perform differential expression analysis and retrieve differentially expressed protein-coding genes
#Usage: analyse_DEG(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile)
#Arguments:
# - GSEaccession: GSE dataset accession number
# - gset: Matrix data retrieved from GEOQuery (from get_DEG function)
# - sampleGrouping: Sample groupings (from get_DEG function)
# - groupDescription: Sample group description (from get_DEG function)
# - summaryFile: Summary file to record output
# Value: List of differentially expressed genes

analyse_deg <- function(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile){
    library(stringr)
    library(umap)
    library(pheatmap)
    library(ggplot2)
    library(dplyr)
    library(igraph)
    library(VennDiagram)
    library(ggnewscale)
    library(limma)
    summaryFile <- paste(GSEaccession,"/","Summary.txt",sep="")
    gene_annot <- fData(gset)
    pc_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_annot$"Gene ID", c("ENTREZID", "GENETYPE"), keytype = "ENTREZID")
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
            print ("Get DE protein-coding genes")
            MapGenes <- paste(deg_list, collapse="\n")
            write(MapGenes,file=paste(GSEaccession,"/","DEGenes.txt",sep=""))
            # summarize test results as "up", "down" or "not expressed"
            print ("Summarize")
            dT<-decideTests(fit2, lfc = 2.0,p.value=0.05,method="separate",adjust.method="BH")
            write(paste("** ",GSEaccession," **",sep=""),file=summaryFile,append=TRUE)
            write(groupDescription,file=summaryFile,append=TRUE)
            write("Summary:",file=summaryFile,append=TRUE)
            write(knitr::kable(summary(dT)),file=summaryFile,append=TRUE)
            pdf(paste(GSEaccession,"/","Plots.pdf",sep=""))
            # Venn diagram of results
            print ("Venn")
            vennDiagram(dT, circle.col=palette())
            # create Q-Q plot for t-statistic
            print ("Q-Q Plot")
            t.good <- which(!is.na(fit2$F)) # filter out bad probes
            qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
            #get volcano plot
            print ("Volcano plot")
            p_cutoff <- 0.05
            fc_cutoff <- 2.0
            ggplot(new_results,aes(x = logFC, y=B)) + geom_point()
            new_results %>%
              mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
            # Visualize and quality control test results.
            # Build histogram of P-values for all genes. Normal test
            # assumption is that most genes are not differentially expressed.
            print ("Histogram")
            tT2 <- topTable(fit2, adjust="BH", sort.by="p", number=Inf)
            hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj", ylab = "Number of genes", main = "P-adj value distribution")
            # volcano plot (log P-value vs log fold change)
            colnames(fit2) # list contrast names
            ct <- 1        # choose contrast of interest
            print ("Volcano plot")
            volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20, highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
            # MD plot (log fold change vs mean log expression)
            # highlight statistically significant (p-adj < 0.05) probes
            print ("MD Plot")
            plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
            abline(h=0)
            ################################################################
            # General expression data analysis
            ex <- exprs(gset2)
            # box-and-whisker plot
            print ("box-and-whisker plot")
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
                print ("Top 20")
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
                    print ("Heatmap")
                    pheatmap(gene_matrix, labels_row = gene_names, scale="row")
                    #write(paste(gene_names,collapse="\n"),file=paste(GSEaccession,"/","Top20DEGenes.txt",sep=""))
                    top20genesx <- new_results[new_results$ID %in% ids_of_interest,]
                    top20genes <- subset(top20genesx, select=c("ID","logFC","Gene.ID","Gene.symbol","Gene.title"))
                    top20genes$DE <- ""
                    top20genes$DE[top20genes$DE > 1.9] <- "Up"
                    top20genes$DE[top20genes$DE < - 1.9] <- "Down"
                    write.table(top20genes, file=paste(GSEaccession,"/","Top20DEGenes.txt", quote=FALSE,sep="\t",row.names = FALSE)
                }, error=function(e)
                {
                    print("Error: Unable to create heatmap.")
                })

            }
            dev.off()
            degsx_list <- deg_list
        }
        else {
            print (paste(GSEaccession,": No protein-coding genes identified.",collapse=" "))
            write(paste("** ",GSEaccession,": No protein-coding genes identified. **\n",collapse=" "),file=summaryFile,append=TRUE)
            write(groupDescription,file=summaryFile,append=TRUE)
            degsx_list <- list()
        }
    } else {
        print ("No differentially expressed genes found in the GSE data set")
        write(paste("** ",GSEaccession,": No differentially expressed genes found. **\n",collapse=" "),file=summaryFile,append=TRUE)
        degsx_list <- list()
    }
    return (degsx_list)
}

