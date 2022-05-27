\name{get_phegeni}
\alias{get_phegeni}
\title{Retrieves genes from PheGenI}
\description{Retrieves related genes for a given query disease in PheGenI through text processing}
\usage{get_phegeni(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [phegeni_genes.txt] and a summary file}
\references{
https://www.ncbi.nlm.nih.gov/gap/phegeni
}
\author{Nur Syatila Ab Ghani}
