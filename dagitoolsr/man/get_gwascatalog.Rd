\name{get_gwascatalog}
\alias{get_gwascatalog}
\title{Retrieves genes from GWAS Catalog}
\description{Retrieves related genes for a given query disease in GWAS Catalog using gwasrapidd package}
\usage{get_gwascatalog(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [gwas_genes.txt] and a summary file}
\references{
https://github.com/ramiromagno/gwasrapidd
}
\author{Nur Syatila Ab Ghani}
