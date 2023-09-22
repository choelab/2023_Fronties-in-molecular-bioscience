library(data.table);library(dplyr)
library(DEP)
library(ggplot2);library(ggrepel)

# library(pheatmap)
# require(SummarizedExperiment)
################################################################
require(httr)
require(stringr)
httr::set_config(config(ssl_verifypeer = 0L))

#library(jsonlite)
uniprot_mapping <- function(ids) {
  uri <- 'https://rest.uniprot.org/uniprotkb/search?size=1&&query='
  #uri <- 'https://www.uniprot.org/uniprotkb/search?query='
  #idStr <- paste(ids, collapse="+or+")
  format <- '&fields=accession%2Cgene_names&format=tsv'
  #fullUri <- paste0(uri,idStr,format)
  fullUri <- paste0(uri,ids,format)
  #dat <- read.delim(fullUri)
  dat <- GET(fullUri,
             accept_json(),
             add_headers(Accept = 'application/json'))
  return(dat)
}


require(pbapply)
pboptions(type="txt", char=":")

require(reshape2)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


pheatmaps <- function(deps, k = 4){

  data_results <- get_results(deps)
  require(SingleCellExperiment)
  require(SummarizedExperiment)
  # selGene <- data_results %>% 
  #   filter(X100.cell_vs_Bulk_p.adj < 0.1 & abs(X100.cell_vs_Bulk_ratio) > log2(2))
  #SummarizedExperiment::rowData(dep, use.names = TRUE)$'significant' <- data_results$name %in% selGene$name
  # SummarizedExperiment::rowData(dep, use.names = TRUE)$'significant' <- SummarizedExperiment::rowData(dep, use.names = TRUE)$name %in% selGene$name
  p<-plot_heatmap(deps,  type = "centered", show_column_dend = F, row_dend_reorder = F, kmeans = TRUE, k= k)  
  return(p)

}

get.ppiNCBI <- function(g.n) {
  require(XML)
  ppi <- data.frame()
  for(i in 1:length(g.n)){
    o <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", g.n[i], sep=''))
    # check if interaction table exists
    exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
    if(exist){
      p <- getNodeSet(o, "//table")
      ## need to know which table is the good one
      for(j in 1:length(p)){
        int <- readHTMLTable(p[[j]])
        if(colnames(int)[2]=="Interactant"){break}
      }
      ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
    }
    # play nice! and avoid being kicked out from NCBI servers
    Sys.sleep(1)
  }
  if(dim(ppi)[1]>0){
    ppi <- unique(ppi)
    print(paste(dim(ppi)[1], "interactions found"))
    return(ppi)
  } else{
    print("No interaction found")
  }
}

