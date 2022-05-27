\name{get_malacards}
\alias{get_malacards}
\title{Retrieves genes from MalaCards}
\description{Retrieves related genes for a given query disease in MalaCards using mala and AutoSeed package}
\usage{get_malacards(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [malacards_genes.txt] and a summary file}
\references{
https://www.malacards.org/, https://cran.r-project.org/web/packages/Autoseed/index.html
}
\author{Nur Syatila Ab Ghani}
