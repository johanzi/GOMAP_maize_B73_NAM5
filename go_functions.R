
# Objects needed for GO analysis
TERM2NAME <- readRDS("TERM2NAME.rds")
#TERM2GENE <-  readRDS("TERM2GENE.rds")
TERM2GENE <- readRDS("TERM2GENE_Fattel_2024.rds")

# Gather all Terms for the GO_search function
TERM2NAME_ALL <- do.call("rbind", TERM2NAME)

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
  TERM2GENE_BP <- TERM2GENE %>% dplyr::filter(GO %in% TERM2NAME$BP$go_id)
  TERM2GENE_CC <- TERM2GENE %>% dplyr::filter(GO %in% TERM2NAME$CC$go_id)
  TERM2GENE_MF <- TERM2GENE %>% dplyr::filter(GO %in% TERM2NAME$MF$go_id)
  
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

  # Return final list of enrichResult class objects
  return(ego_final)
  
}

# Additional function to explore clusterProfile enrich output
go_search <- function(method=c("gene2GO","GO2gene","GO2term","term2GO","GO2ontology"),id){
  
  if(method=="gene2GO"){
    GO <- TERM2GENE %>% dplyr::filter(gene==id) %>% pull(GO)
    df <- TERM2NAME_ALL %>% dplyr::filter(go_id %in% GO)
    return(df)
  }
  
  if(method=="GO2gene"){
    df <- TERM2GENE %>% dplyr::filter(GO==id)
    return(df)
  }
  
  if(method=="GO2term"){
    df <-TERM2NAME_ALL %>% dplyr::filter(go_id==id)
    return(df)
  }
  
  if(method=="term2GO"){
    df <- TERM2NAME_ALL %>% dplyr::filter(Term==id)
    return(df)
  }
  
  if(method=="GO2ontology"){
    if(sum(TERM2NAME$BP$go_id==id)==1){
      message("BP") } 
    else if (sum(TERM2NAME$CC$go_id==id)==1){
      message("CC") } 
    else if (sum(TERM2NAME$MF$go_id==id)==1){
      message("MF") } else {
      message("GO ID not found in any ontology")
      }
  }
}

# Testing examples for go_search function
#go_search(method="gene2GO", "Zm00001eb248920")
#go_search(method="GO2gene", "GO:0000011")
#go_search(method="GO2term", "GO:0000062")
#go_search(method="term2GO", "reproduction")
#go_search(method="GO2ontology", "GO:0000062")


# Function to convert list of enrichResult objects into one dataframe
enrichResult2dataframe <-  function(output_ego_analysis){
  
  if(!is.null(output_ego_analysis$ego_BP)){
    df_BP <- as.data.frame(output_ego_analysis$ego_BP@result)
    df_BP$ontology <- "BP"
  }
  
  if(!is.null(output_ego_analysis$ego_CC)){
    df_CC <- as.data.frame(output_ego_analysis$ego_CC@result)
    df_CC$ontology <- "CC"
  }
  
  if(!is.null(output_ego_analysis$ego_MF)){
    df_MF <- as.data.frame(output_ego_analysis$ego_MF@result)
    df_MF$ontology <- "MF"
  }
  
  # Bind all. Use get0 so that there is no error when binding non-existing
  # data frames
  df_ego <- rbind(get0("df_BP"), get0("df_CC"), get0("df_MF"))
  
  # Check if df_ego contains anything
  if(!is.null(df_ego)){
    # Put ontology in first column
    df_ego <- df_ego %>% dplyr::select(ontology, everything())
    return(df_ego)
  } else {
    message("No ontology contained a significant hit.")
  }
}
