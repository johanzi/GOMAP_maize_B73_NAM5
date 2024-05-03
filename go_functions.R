
# Objects needed for GO analysis
TERM2GENE <-  readRDS("TERM2GENE.rds")
TERM2NAME <- readRDS("TERM2NAME.rds")

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

  # Keep in each data frames only the defined ontology (BP, CC, or MF)
  TERM2GENE_BP <- TERM2GENE %>% filter(GO %in% TERM2NAME$BP$go_id)
  TERM2GENE_CC <- TERM2GENE %>% filter(GO %in% TERM2NAME$CC$go_id)
  TERM2GENE_MF <- TERM2GENE %>% filter(GO %in% TERM2NAME$MF$go_id)
  
  # Check if all objects are loaded
  if(!exists("TERM2GENE")){stop("TERM2GENE object not loaded")}
  if(!exists("TERM2NAME")){stop("TERM2NAME object not loaded")}

  # Perform enrichment analysis on each ontology separately
  ego_BP <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=TERM2GENE_BP,
    TERM2NAME = TERM2NAME$BP
  )
  
  ego_CC <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=TERM2GENE_CC,
    TERM2NAME = TERM2NAME$CC
  )
  
  
  ego_MF <-enricher(
    gene=vector_geneID,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff =0.2,
    TERM2GENE=TERM2GENE_MF,
    TERM2NAME = TERM2NAME$MF
  )
  
  # Create a list of ego results
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
    GO <- TERM2GENE %>% filter(gene==id) %>% pull(GO)
    df <- TERM2NAME$ALL %>% filter(go_id %in% GO)
    return(df)
  }
  
  if(method=="GO2gene"){
    return(TERM2GENE %>% filter(GO==id))
    #geneID <- TERM2GENE %>% filter(GO==id) %>% pull(gene)
    #print(geneID)
  }
  
  if(method=="GO2term"){
    return(TERM2NAME$ALL %>% filter(go_id %in% id))
  }
  
  if(method=="term2GO"){
    return(TERM2NAME$ALL %>% filter(Term %in% id))
  }
  
}
