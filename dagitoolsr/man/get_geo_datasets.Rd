\name{get_geo_datasets}
\alias{get_geo_datasets}
\title{Retrieves GEO datasets for given MeSH terms}
\description{The function retrieves GSE accession numbers in NCBi GEO Datasets for given MeSH terms using rentrez package}
\usage{get_geo_datasets(MeSHterms, summaryFile)}
\arguments{
\item{MeSHterms}{MeSHterms in a list}
\item{summaryFile}{Summary file to record output}
}
\value{List of GSE accession numbers and a summary file}
\references{
https://cran.r-project.org/web/packages/rentrez/index.html
}
\author{Nur Syatila Ab Ghani}
