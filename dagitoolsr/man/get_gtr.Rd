\name{get_gtr}
\alias{get_gtr}
\title{Retrieves genes from Gene Testing Registry (GTR)}
\description{Retrieves related genes for a given query disease in GTR through html processing}
\usage{get_gtr(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [gtr_genes.txt] and a summary file}
\references{
https://www.malacards.org/
}
\author{Nur Syatila Ab Ghani}
