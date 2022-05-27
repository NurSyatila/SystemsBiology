\name{get_medgen}
\alias{get_medgen}
\title{Retrieves genes from MedGen}
\description{Retrieves related genes for a given query disease in MedGen using rentrez package}
\usage{get_medgen(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [medgen_genes.txt] and a summary file}
\references{
https://cran.r-project.org/web/packages/rentrez/index.html
}
\author{Nur Syatila Ab Ghani}
