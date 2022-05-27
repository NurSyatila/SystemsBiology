\name{get_mesh_terms}
\alias{get_mesh_terms}
\title{Retrieves MeSH terms from MeSH database}
\description{Retrieves MeSH terms related to a given query disease using rentrez package}
\usage{get_mesh_terms(queryTerms,summaryFile)}
\arguments{
\item{queryTerms}{Query terms in the form of list}
\item{summaryFile}{Summary file to record output}
}
\value{List of MeSH terms and a summary file}
\references{
https://cran.r-project.org/web/packages/rentrez/index.html
}
\author{Nur Syatila Ab Ghani}
