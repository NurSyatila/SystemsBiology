\name{get_enrichment_results}
\alias{get_enrichment_results}
\title{Retrieves enrichment results}
\description{The function performs functional enrichment using ClusterProfiler package and generate visualization plots for a given gene list. Available functional terms for enrichment: BP, MF, CC, KEGG, DO. For KEGG enrichment, top enriched pathways with mapped genes will be automatically retreived and saved in 'Pathview' directory.}
\usage{get_enrichment_results(geneList,ORAterm, summaryFile)}
\arguments{
\item{geneList}{Input gene list (Entrez identifiers separated by new line))}
\item{ORAterm}{Functional enrichment terms, "BP":Biological Process,"MF":Molecular Function,"CC":Cellular Component,"KEGG":KEGG Pathway,"DO":Disease Ontology)}
\item{summaryFile}{Summary file to record output}
}
\value{Functional Enrichment results in CSV-formatted files, visualization plots, Mapped pathway files ['Pathview' directory] and a summary file}
\references{
https://bioconductor.org/packages/clusterProfiler/
}
\author{Nur Syatila Ab Ghani}
