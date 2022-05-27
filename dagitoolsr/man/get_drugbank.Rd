\name{get_drugbank}
\alias{get_drugbank}
\title{Retrieves genes from DrugBank}
\description{Retrieves related genes for a given query disease in DrugBank using AutoSeed package}
\usage{get_drugbank(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [drugbank_genes.txt] and a summary file}
\references{
https://cran.r-project.org/web/packages/Autoseed/index.html
}
\author{Nur Syatila Ab Ghani}
