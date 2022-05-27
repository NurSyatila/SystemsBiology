\name{analyse_ppi_network}
\alias{analyse_ppi_network}
\title{Analyse PPI network from STRINGdb}
\description{The automated function analyses protein-protein interaction using STRINGdb, where it (1) retrieves PPI network from STRINGdb, (2) filter isolated nodes, (3) identify top genes, (4) perform functional enrichment analysis using get_enrichment_results function, and (5) perform network clustering using get_gene_modules function )}
\usage{analyse_ppi_network(geneList,summaryFile)}
\arguments{
\item{geneList}{Input gene list (Entrez identifiers separated by new line)}
\item{summaryFile}{Summary file to record output}
}
\value{Visualization plots [string_ppi_networks.pdf], filtered genes ['filtered_genes' directory], top 10/20 genes ['top_genes' directory], and a summary file}
\references{
https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html
}
\author{Nur Syatila Ab Ghani}
