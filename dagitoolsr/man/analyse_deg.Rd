\name{analyse_deg}
\alias{analyse_deg}
\title{Analyse differential expression data set}
\description{The function performs differential expression analysis and retrieves differentially expressed protein-coding genes. Visualization plots will be automatically retrieved.}
\usage{analyse_deg(GSEaccession,gset,sampleGrouping,groupDescription,summaryFile)}
\arguments{
\item{GSEaccession}{GSE dataset accession number}
\item{gset}{Matrix data retrieved from GEOQuery (from get_DEG function)}
\item{sampleGrouping}{Sample groupings (from get_DEG function)}
\item{groupDescription}{Sample group description (from get_DEG function)}
\item{summaryFile}{Summary file to record output}
}
\value{List of differentially expressed genes (protein-coding) for given GEO dataset [GSE_genes.txt], top 20 genes [top_20_DE_genes.txt], visualization plots [visualization_plots.pdf], and a summary file }
\author{Nur Syatila Ab Ghani}
