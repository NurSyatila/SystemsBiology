\name{get_gene_modules}
\alias{get_gene_modules}
\title{Retrieves gene modules}
\description{The function performs network clustering using walktrap community clustering in igraph package (walktrap.community function).}
\usage{get_gene_modules(geneNetwork, summaryFile)}
\arguments{
\item{geneNetwork}{Gene network for clustering}
\item{summaryFile}{Summary file to record output}
}
\value{data frame containing information on gene modules, a summary file}
\references{
https://cran.r-project.org/web/packages/igraph/index.html
}
\author{Nur Syatila Ab Ghani}
