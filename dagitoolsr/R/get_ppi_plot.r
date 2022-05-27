#============= PPI_Network ====================
#Function: Plot PPI network
#Usage: get_PPI_plot(geneNetwork,graphTitle,geneAnnotation)
# Arguments:
# - geneNetwork: Gene network for plotting (PPI network derived from STRINGdb)
# - graphTitle: Graph title
# - geneAnnotation: Gene annotation
# Value: graph plot

get_ppi_plot <- function(geneNetwork,graphTitle,geneAnnotation){
  library(ggplot2)
  library(igraph)
  geneAnnotation_sorted <- geneAnnotation[match(V(geneNetwork)$name, geneAnnotation$V1),]
  V(geneNetwork)$name <- geneAnnotation_sorted$V2
  V(geneNetwork)$vertex.frame.color <- "white"
  edgeweights <- E(geneNetwork)$weight * 2.0
  suppressWarnings(plot(
    geneNetwork,
    layout=layout.fruchterman.reingold,
    edge.curved=FALSE,
    vertex.label.color="black",
    asp=FALSE,
    vertex.label.cex=0.6,
    vertex.shape="circle",
    vertex.color=rainbow(betweenness(geneNetwork), start=0, end=2/6),
    edge.width=edgeweights,
    edge.arrow.mode=0,
    main=graphTitle))
}