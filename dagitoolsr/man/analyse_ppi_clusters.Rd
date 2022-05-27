\name{analyse_ppi_clusters}
\alias{analyse_ppi_clusters}
\title{Analyse PPI clusters}
\description{The function analyses gene modules resulting from PPI clustering and perform KEGG enrichment for individualgene modules.}
\usage{analyse_ppi_clusters(dfClusters,summaryFile)}
\arguments{
\item{dfClusters}{Data frame of clusters derived from clustering (from get_gene_modules function)}
\item{summaryFile}{Summary file to record output}
}
\value{List of gene symbols for individual gene modules (gene_clusters/ directory), a summary file}
\references{
https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html
}
\author{Nur Syatila Ab Ghani}
