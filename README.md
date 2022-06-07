# SystemsBiology

Part 1: Gene Search (GeneSearch.r)
Description: A user-guided R script to search for disease-associated genes from publicly available data sets: ClinVar, DISEASES, GTR, GWAS Catalog, MedGen, OMIM, PheGenI, and DisGeNET
User Input: A disease or keywords separated by '|' (e.g. colorectal cancer|colorectal carcinoma)
Options: To include/exclude DisGeNET search (Email and Password to DisGeNET account are required)
Output: A directory named GeneSearch containing: 
  a) <Dataset>Summary.csv: A summary for individual data set search
  b) <Dataset>_Genes.txt: A list of genes (symbols, entrez IDs, and gene type) for individual data set search
  c) <Dataset>_PCGenes.txt: A list of Entrez IDs for protein-coding genes (separated by new lines) for individual data set search 
  d) GeneSearch.pdf: A graphical summary of protein-coding genes in bar chart
Functions:
  a) check_libraries: Check libraries required for all functions in GeneSearch.r. Install packages if not available. Parameters: 
      ~ pack_cran:  a list of packages stored in CRAN repository
      ~ pack_bioc: a list of packages under Bioconductor
      ~ pack_github:  a list of packages from GitHub 
  b) get PCGenes: Extract protein-coding genes based on annotation in the 'org.Hs.eg.db' library. Parameters: 
      ~ genelist: A list of genes (Entrez IDs or symbols)
      ~ gene_file: A filename to store the list of genes
      ~ pcgene_file: A filename to store protein-coding genes
      ~ key_type: type of key for gene in genelist ('ENTREZID' or 'SYMBOL')
  c) get_clinvar: Retrieve disease-associated genes from ClinVar. Parameter:
      ~ queryTerms: a list of diseases / keywords for search
  d) get_diseasesdb: Retrieve disease-associated genes from DISEASES database. Parameter:
      ~ queryTerms: a list of diseases / keywords for search
  c) get_meshterms: Retrieve MeSH terms. Parameters:
      ~ queryTerms: a list of diseases / keywords for search
      ~ type: Type of output to be returned, "term" (for MeSH terms) or "id" (for MeSH identifier)
  d) get_disgenet: Retrieve disease-associated genes from DisGeNET. Parameters:
      ~ queryTerms: a list of diseases / keywords for search
      ~ userEmail: User email to DisGeNET account
      ~ userPassword: User password to DisGeNET account
  e) get_gtr: Retrieve disease-associated genes from GTR. Parameter:
      ~ queryTerms: a list of diseases / keywords for search
  f) get_gwascatalog: Retrieve disease-associated genes from GWAS Catalog. Parameter:
      ~ queryTerms: a list of diseases / keywords for search
  g) get_medgen: Retrieve disease-associated genes from MedGen. Parameter:
      ~ queryTerms: a list of diseases / keywords for search
  h) get_omim: Retrieve disease-associated genes from OMIM. Parameter:
      ~ queryTerms: a list of diseases / keywords for search
  i) get_phegeni: Retrieve disease-associated genes from PheGenI. Parameter:
      ~ queryTerms: a list of diseases / keywords for search  
  
Part 2: Differential Expression (DifferentialExpression.r)
Description: A user-guided R script to retrieve disease-associated genes (differentially expressed genes) from Gene Expression data sets retrieved from GEO DataSets.
User Input: A disease or keywords separated by '|' (e.g. colorectal cancer|colorectal carcinoma)
Output: A directory named DifferentialExpression containing:
  a) GSESummary.csv: A summary of retrieved GEO data sets - GDS Accession, GPL Platform ID, GSE Accession, Series Title, Subset Info
  b) GSESummaryFiltered.csv: A summary of retrieved GEO data sets that contain samples derived from disease state (samples with annotation under 'disease state')
  c) GSESummary.txt: A summary of gene expression data analysis, sample groupings and gene search for individual data sets.
  d) Individual folder for GSE accession with the following information:
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
  Functions:
   a) check_libraries: A function to check libraries required for all functions in GeneSearch.r. Install packages if not available. Parameters: 
        ~ pack_cran:  a list of packages stored in CRAN repository
        ~ pack_bioc: a list of packages under Bioconductor
        ~ pack_github:  a list of packages from GitHub 
   b) get_meshterms: Retrieve MeSH terms. Parameters:
        ~ queryTerms: a list of diseases / keywords for search
        ~ type: Type of output to be returned, "term" (for MeSH terms) or "id" (for MeSH identifier)
   c) get_gse_datasets: Retrieve GEO Data sets from NCBI GEO DataSet. Parameter:
        ~ queryTerms: a list of diseases / keywords for search
   d) analyse_deg: Analyse gene expression data and save results. Parameters:
        ~ GSEaccession: GSE accession for analysis 
        ~ gset: A GEO object for a query GEO series accession from getGEO function in GEOquery library
        ~ sampleGrouping: Groups assigned for individual samples
        ~ groupDescription: Description of individual samples
        ~ summaryFile: A filename to store the summary of the gene expresion analysis and gene search
   e) get_deg: Determine sample groupings from metadata and use analyse_deg function to obtain differentially expressed genes. Parameters:
        ~ queryTerms: a list of diseases / keywords for search
        ~ GSEaccession: GEO Series accession for analysis
        ~ GSEplatform: Corresponding GPL / Platform identifier
        ~ GSETitle: Title of the GEO data set
    
