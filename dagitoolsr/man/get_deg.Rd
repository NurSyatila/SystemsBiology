\name{get_deg}
\alias{get_deg}
\title{Retrieves differentially expressed genes for a given GEO dataset}
\description{The automated function retrieves GEO dataset matrix data and analyses sample groupings for further processing using analyse_DEG function. Only GEO dataset with comparable sample groups will be subjected to DEG analysis.}
\usage{get_deg(queryTerms,summaryFile)}
\arguments{
\item{GSEaccession}{GSE dataset accession number)}
\item{summaryFile}{Summary file to record output}
}
\value{List of differentially expressed genes from analyse_DEG function and a summary file}
\references{
https://www.ncbi.nlm.nih.gov/gds
}
\author{Nur Syatila Ab Ghani}
