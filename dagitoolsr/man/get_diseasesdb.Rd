\name{get_diseasesdb}
\alias{get_diseasesdb}
\title{Retrieves genes from DISEASES database (Jensen)}
\description{Retrieves related genes for a given query disease in MalaCards through text processing}
\usage{get_diseasesdb(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [diseasesdb_genes.txt] and a summary file}
\references{
https://diseases.jensenlab.org/
}
\author{Nur Syatila Ab Ghani}
