\name{get_clinvar}
\alias{get_clinvar}
\title{Retrieves genes from ClinVar}
\description{Retrieves related genes for a given query disease in ClinVar using rentrez package}
\usage{get_clinvar(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [clinvar_genes.txt] and a summary file}
\references{
https://cran.r-project.org/web/packages/rentrez/index.html
}
\author{Nur Syatila Ab Ghani}
