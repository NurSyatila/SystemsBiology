\name{get_disgenet}
\alias{get_disgenet}
\title{Retrieves genes from DisGeNET}
\description{Retrieves related genes for a given query disease based on disease-gene associations in DisGENET using disgenet2r package}
\usage{get_disgenet(queryTerms,summaryFile)}
\arguments{
\item{MeSHTerms}{MeSH terms in the form of list (Query terms may include DO terms and MeSH terms)}
\item{userEmail}{DisGeNET account ID / Email}
\item{MeSHTerms}{DisGeNET account password}
\item{summaryFile}{Summary file to record output}
}
\value{List of disease-related genes [disgenet_genes.txt] and a summary file}
\references{
https://www.disgenet.org/static/disgenet2r/disgenet2r.html
}
\author{Nur Syatila Ab Ghani}
