\name{get_kegg_enrichment}
\alias{get_kegg_enrichment}
\title{Retrieves KEGG enrichment results}
\description{The function performs KEGG pathway enrichment using ClusterProfiler package and generate visualization plots for a given gene list.}
\usage{get_kegg_enrichment(geneList,summaryFile)}
\arguments{
\item{geneList}{Input gene list (Entrez identifiers separated by new line))}
\item{summaryFile}{Summary file to record output}
}
\value{Functional Enrichment results in CSV-formatted files and a summary file}
\references{
https://bioconductor.org/packages/clusterProfiler/
}
\author{Nur Syatila Ab Ghani}
