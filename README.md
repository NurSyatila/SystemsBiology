# SystemsBiology

## Part 1: Gene Search (GeneSearch.r)  
##### Description: A user-guided R script to search for disease-associated genes from publicly available data sets: ClinVar, DISEASES, GTR, GWAS Catalog, MedGen, OMIM, PheGenI, and DisGeNET  
##### User Input: A disease or keywords separated by '|' (e.g. colorectal cancer|colorectal carcinoma)  
##### Options: To include/exclude DisGeNET search (Email and Password to DisGeNET account are required)  
##### Output: A directory named GeneSearch containing:   
  - <Dataset>Summary.csv: A summary for individual data set search  
  - <Dataset>_Genes.txt: A list of genes (symbols, entrez IDs, and gene type) for individual data set search  
  - <Dataset>_PCGenes.txt: A list of Entrez IDs for protein-coding genes (separated by new lines) for individual data set search   
  - GeneSearch.pdf: A graphical summary of protein-coding genes in bar chart  
##### Functions:   
  - check_libraries: Check required libraries and install if not available   
  - get PCGenes: Extract protein-coding genes based on annotation in the 'org.Hs.eg.db' library   
  - get_clinvar: Retrieve disease-associated genes from ClinVar  
  - get_diseasesdb: Retrieve disease-associated genes from DISEASES database    
  - get_meshterms: Retrieve MeSH terms    
  - get_disgenet : Retrieve disease-associated genes from DisGeNET. User email and password to DisGeNET account are required.   
  - get_gtr: Retrieve disease-associated genes from GTR    
  - get_gwascatalog: Retrieve disease-associated genes from GWAS Catalog    
  - get_medgen: Retrieve disease-associated genes from MedGen    
  - get_omim: Retrieve disease-associated genes from OMIM    
  - get_phegeni: Retrieve disease-associated genes from PheGenI    
  
## Part 2: Differential Expression (DifferentialExpression.r)  
##### Description: A user-guided R script to retrieve disease-associated genes (differentially expressed genes) from Gene Expression data sets retrieved from GEO DataSets.  
##### User Input: A disease or keywords separated by '|' (e.g. colorectal cancer|colorectal carcinoma)  
##### Output: A directory named DifferentialExpression containing:  
  - GSESummary.csv: A summary of retrieved GEO data sets - GDS Accession, GPL Platform ID, GSE Accession, Series Title, Subset Info  
  - GSESummaryFiltered.csv: A summary of retrieved GEO data sets that contain samples derived from disease state (samples with annotation under 'disease state')  
  - GSESummary.txt: A summary of gene expression data analysis, sample groupings and gene search for individual data sets.  
  - Individual folder for GSE accession with the following information:  
      - FilteredDE.csv: Gene expresssion data for differentially expressed protein-coding genes (Log2 fold change < 2.0 | > 2.0)  
      - DEGenes.txt: A list of Entrez IDs for differentially expressed protein-coding genes (separated by new lines)  
      - Top20DEGenes.txt: A summary of top 20 differentially expressed protein-coding genes in tabulated form  
      - Plots.pdf: A graphical summary of gene expression analysis;  
          ~ Venn diagram of results  
          ~ Q-Q plot for t-statistics  
          ~ p-value histogram  
          ~ Volcano plot (log P-value vs. log fold change)  
          ~ MD plot (MD plot (log fold change vs. mean log expression)  
          ~ Box-and-whisker plot  
          ~ Expression value distribution (density vs. intensity)  
          ~ Mean variance trend (residual standard deviation vs. average log expression)  
          ~ Heat map for top 20 differentially expressed genes  
  ##### Functions:  
   - check_libraries: Check required libraries and install if not available  
   - get_meshterms: Retrieve MeSH terms  
   - get_gse_datasets: Retrieve GEO Data sets from NCBI GEO DataSet  
   - analyse_deg: Analyse gene expression data and save results for a given GEO Series accession  
   - get_deg: Determine sample groupings from metadata and use analyse_deg function to obtain differentially expressed genes for a given GEO Series accession     
      
## Part 3: Functional Enrichment (FunctionalEnrichment.r)  
##### Description: A user-guided R script to get protein-protein interaction (PPI) network and perform functional enrichment analysis / over-representation analysis for a set of genes.  
##### User Input: A filename containing Entrez identifiers (separated by new line)  
##### Output: A directory named FunctionalEnrichment containing:  
  - PPINetworkAnalysis.txt: A summary of functional enrichment analysis for a given gene list and identified gene clusters  
  - STRING_PPI_Network.pdf: A graphical summary of STRING protein-protein interaction network and functional enrichment analysis 
  - FilteredGenes.txt: A list of genes from STRING protein-protein interaction network in a tabulated form  
  - Top_PPI_genes.txt: A list of top 20 or 10 genes from STRING PPI network (with highest number of interactions)
  - FilteredGenes directory containing:  
      - Functional enrichment results (CSV format) and plots (PDF)
      - Pathview directory containing mapped KEGG pathway for top 5 enriched pathways
  - TopGenes directory containing:  
      - Functional enrichment results (CSV format) and plots (PDF)  
      - Pathview directory containing mapped KEGG pathway for top 5 enriched pathways  
  - GeneClusters directory containing:
      - Individual cluster file with a list of genes
  ##### Functions:  
   - check_libraries: Check required libraries and install if not available  
   - get_ppi_plot: Plot a graphical network using igraph for a given PPI network    
   - get_kegg_enrichment: Perform KEGG enrichment analysis for a given gene list     
   - get_gene_modules: Perform gene clustering for a given STRING PPI network  
   - get_enrichment_results: Perform over-representation analysis for a given gene list and functional term (BP, MF, CC, KEGG or DO)  
   - analyse_ppi_clusters: Get annotation and perform KEGG enrichment analysis for individual gene modules from gene network clustering  
   - analyse_ppi_network: Retrieve PPI network from STRING, filtered genes (genes with isolated nodes are discarded), top 20/10 genes, and identify gene modules   
