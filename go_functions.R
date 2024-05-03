
# Objects needed for GO analysis
GO_analysis_data <-  readRDS("path/to/GO_analysis_data.rds")
term2name <- readRDS("path/to/term2name.rds")

# Function to perform ego, provide as argument a vector containing geneID (characters type)
ego_analysis <- function(vector_geneID){
  
  require(tidyverse)
  require(clusterProfiler)
  require(DOSE)
  require(enrichplot)
  require(xlsx)
  require(readxl)
  require(AnnotationDbi)
  require(GO.db)

  # Keep in each dataframes only the defined ontology (BP, CC, or MF)
  GO_analysis_data_BP <- GO_analysis_data %>% filter(GO %in% term2name$BP$go_id)
  GO_analysis_data_CC <- GO_analysis_data %>% filter(GO %in% term2name$CC$go_id)
  GO_analysis_data_MF <- GO_analysis_data %>% filter(GO %in% term2name$MF$go_id)
  
  # Check if all objects are loaded
  if(!exists("GO_analysis_data")){stop("GO_analysis_data object not loaded")}
  if(!exists("term2name")){stop("term2name object not loaded")}

  # Perform enrichment analysis on each ontology
  ego_BP <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=GO_analysis_data_BP,
    TERM2NAME = term2name$BP
  )
  
  ego_CC <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=GO_analysis_data_CC,
    TERM2NAME = term2name$CC
  )
  
  
  ego_MF <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=GO_analysis_data_MF,
    TERM2NAME = term2name$MF
  )
  
  
  ego_list <- list(ego_BP = ego_BP, ego_CC = ego_CC, ego_MF = ego_MF)
  
  # Create variable Enrichment factor (Count / number of gene in the given GO)
  ego_list2 <- lapply(ego_list, function(x) mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))))
  
  # Create variable Fold enrichment (GeneRatio / BgRatio)
  ego_final <- lapply(ego_list2, function(x) mutate(x, FoldEnrich = parse_ratio(GeneRatio) / parse_ratio(BgRatio)))
  
  return(ego_final)
  
}

# Additional practical functions to explore clusterProfile enrich output
go_search <- function(method="gene2GO",id){
  if(method=="gene2GO"){
    GO <- GO_analysis_data %>% filter(gene==id) %>% pull(GO)
    df <- term2name$ALL %>% filter(go_id %in% GO)
    return(df)
  }
  
  if(method=="GO2gene"){
    return(GO_analysis_data %>% filter(GO==id))
    #geneID <- GO_analysis_data %>% filter(GO==id) %>% pull(gene)
    #print(geneID)
  }
  
  if(method=="GO2term"){
    return(term2name$ALL %>% filter(go_id %in% id))
  }
  
  if(method=="term2GO"){
    return(term2name$ALL %>% filter(Term %in% id))
  }
  
}
