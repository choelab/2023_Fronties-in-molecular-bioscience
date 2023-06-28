options(stringsAsFactors = FALSE)

###### Load data

brain_cells.subset <- readRDS("/choelab/Fronties-in-molecular-bioscience-2023/data/brain_cells_merged_including_06m_2023-06-16 23-39-03.RDS")
Idents(object = brain_cells.subset) <- brain_cells.subset@meta.data$SCT_snn_res.0.8

brain_cells.markers <- FindAllMarkers(brain_cells.subset,
                                    only.pos = TRUE, logfc.threshold = 1, 
                                    assay = "RNA")

mg.subsetted <- subset(brain_cells.subset, ident = 13)

###### Dimension reduction and clustering
pcs <- 1:50
mg.subset <- mg.subsetted %>%
                 RunPCA(verbose = FALSE, features = VariableFeatures(mg.subsetted)) %>%
                 RunUMAP(reduction = "pca", umap.method = "umap-learn", dims=pcs) %>%
                 FindNeighbors(reduction = "pca", dims=pcs) %>%
                 FindClusters(resolution =1.8, algorithm = 2)

Idents(object = mg.subset) <- sprintf("%02d",mg.subset@meta.data$SCT_snn_res.1.8)
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '01'))] <- "MG-01"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '02'))] <- "MG-02"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '03'))] <- "MG-03"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '04'))] <- "MG-04"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '05'))] <- "MG-05"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '06'))] <- "MG-06"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '07'))] <- "MG-07"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '08'))] <- "MG-08"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '09'))] <- "MG-09"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '10'))] <- "MG-10"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '11'))] <- "MG-11"
mg.subset@meta.data$cell_type_age[colnames(mg.subset) %in% WhichCells(subset(mg.subset, ident = '12'))] <- "MG-12"

###### Gene scoring
mark <- read.xlsx("/choelab/Fronties-in-molecular-bioscience-2023/data/41467_2022_28473_MOESM8_ESM.xlsx", sheetIndex=1)
DAM <- mark$Disease.associated.Microglia..DAM. %>% na.omit()
ATM <- mark$Axon.Tract.Microglia..ATM. %>% na.omit()

markers <- list()
markers$DAM <- intersect(grep(paste0(c("Lpl", "Cst7","Apoe","Spp1","Itgax","Clec7a","Csf1","Igf1","Axl","Ccl6","Ctsz","Cd9","Trem2","Tyrobp","Ctsl","Cd63","Ank","Ctsd"), collapse = "$|^"),rownames(mg.subset), value = TRUE),gsub("\\-","\\.",DAM))
markers$ATM <- gsub("\\-","\\.",ATM)[c(1:20)]
markers$IRM <- grep(paste0(c("Ifit1", "Ifit3","Oas1a","Oas1g","Oas3","Oasl1","Oasl2","Isg15","Isg20","Irf7","Ly6a","Usp18","Ifitm2","Ifitm3","Csprs","Ccrl2","Stat2","Dhx58","Ifi274l2a","Ifit2","Ifit3"), collapse = "$|^"),rownames(mg.subset), value = TRUE)

require(UCell)
mg.subset.mark <- AddModuleScore_UCell(subset(mg.subset,subset = genetype == names(table(mg.subset$genetype))[2], 
                                              idents = c("01","02","03","10")), features = markers)
signature.names <- paste0(names(markers), "_UCell")

df<-data.frame("DAM"=as.numeric(mg.subset.mark@meta.data$DAM_UCell), 
               "IRM"=as.numeric(mg.subset.mark@meta.data$IRM_UCell), 
                "ATM"=as.numeric(mg.subset.mark@meta.data$ATM_UCell), 
                "cell_type_age"=mg.subset.mark@meta.data$cell_type_age)

rownames(df) <- colnames(mg.subset.mark)

dataf <- df %>%
  group_by(cell_type_age)      %>%  
   summarise(across(everything(), mean),
            .groups = 'drop')  %>%
  as.data.frame()


dataf.dam <-  df[,c(2,5)] %>% 
    group_by(cell_type_age)    %>%
    summarise(across(everything(), mean),
              .groups = 'drop')  %>%
    as.data.frame()

dataf.irm <-  df[,c(3,5)] %>% 
    group_by(cell_type_age)    %>%
    summarise(across(everything(), mean),
              .groups = 'drop')  %>%
    as.data.frame()

dataf.atm <-  df[,c(4,5)] %>% 
    group_by(cell_type_age)    %>%
    summarise(across(everything(), mean),
              .groups = 'drop')  %>%
    as.data.frame()

datafs.dam <- dataf.dam[,2]
names(datafs.dam) <- dataf.dam$cell_type_age
datafs.irm <- dataf.irm[,2]
names(datafs.irm) <- dataf.irm$cell_type_age
datafs.atm <- dataf.atm[,2]
names(datafs.atm) <- dataf.atm$cell_type_age

###### Differentially expressed genes between 5xFAD and WT
mg.subset.DEG <- lapply(levels(Idents(mg.subset)), function(id){
            if(length(colnames(subset(mg.subset, idents = id, subset = genetype == names(table(mg.subset$genetype))[1]))) > 3 & length(colnames(subset(mg.subset, idents = id, subset = genetype == names(table(mg.subset$genetype))[2]))) > 3) {
                res <- FindMarkers(
                        subset(mg.subset, idents = id),
                         group.by = "genetype",
                         ident.1 = names(table(mg.subset$genetype))[1],
                         ident.2 = names(table(mg.subset$genetype))[2],
                        logfc.threshold = 0.00001, 
                        test.use = "MAST",
                        only.pos = FALSE)
                # write.xlsx(res, file = paste0(file.path("/home","choelab","working","brain_cells","results"),"/DEGs_",id,"_5xFAD_vs_WT.xlsx"))
                        } else {
                  res <- c()
                        }
                return(res)  
})

names(mg.subset.DEG) <- levels(Idents(mg.subset))

###### Differentially expressed genes between MG-01, MG-03, MG-10 and MG-02 in WT mice
mg.subset.DE.wt <- list()
for(ids in c("03","01","10")) {
  mg.subset.DE.wt[[ids]] <- FindMarkers(
                          subset(mg.subset, subset= genetype == names(table(mg.subset$genetype))[2]),
                          ident.1 = ids,
                          ident.2 = "02", 
                          logfc.threshold = 0.00001, 
                          assay = "RNA",
                          test.use = "MAST",
                          only.pos = FALSE)
}
names(mg.subset.DE.wt) <- c("03","01","04","05","07","10")

###### pahtway analysis of WT
require(DOSE)
# library(plyr)
# library(dplyr)

require(AnnotationDbi)
#library(org.Hs.eg.db)
require(org.Mm.eg.db)
require(clusterProfiler)

Mm <- org.Mm.eg.db
microglia_vs_HM.wt.cc <- list()
microglia_vs_HM.wt.bp <- list()
microglia_vs_HM.wt.mf <- list()
microglia_vs_HM.wt.reactome <- list()
microglia_vs_HM.wt.kegg <- list()
microglia_vs_HM.wt.wiki <- list()
genenlist <- list()

for(id in c("03","01","10")) {
    de<-c()
    selectedGenes <- c()
    gene_list <-c()
    gene_lists <-c()
    selectedGene <- mg.subset.DE.wt[[id]] %>% dplyr::filter(p_val < 0.01  & abs(avg_log2FC) > 0.25)
    selectedGenes <- selectedGene[!grepl("^mt\\.", rownames(selectedGene)),]   
    gene_lists<- selectedGenes %>% dplyr::select(avg_log2FC)
  
    de<-AnnotationDbi::select(Mm,
                              keys = rownames(gene_lists),
                              columns = c("ENTREZID", "SYMBOL"),
                              keytype = "SYMBOL")
    
    geneLST<-data.frame("SYMBOL"=rownames(gene_lists), "avg_log2FC" = gene_lists$avg_log2FC)
    
    gene_df <- de %>% mutate(rank = rank(de$ENTREZID,  ties.method = "random")) # %>% 
      
    gene_dfs <-merge(gene_df, geneLST, by = "SYMBOL", all=TRUE) %>%   arrange(desc(rank))
    
    
    gene_dfs<-gene_dfs[!is.na(gene_dfs$ENTREZID),]
    
    if(length(gene_dfs$ENTREZID)>0) {
      
      gene.list<-gene_dfs$'avg_log2FC'
      names(gene.list)<-gene_dfs$ENTREZID
      
      
      gene.list = sort(gene.list, decreasing = TRUE)
      gene.list <- na.omit(gene.list)
      genenlist[[id]] <- gene.list
      gse.cc <- gseGO(geneList=gene.list, 
                      ont ="CC",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      microglia_vs_HM.wt.cc[[id]] <- setReadable(gse.cc, 'org.Mm.eg.db', 'ENTREZID')

       gse.bp <- gseGO(geneList=gene.list, 
                      ont ="BP",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
       microglia_vs_HM.wt.bp[[id]] <- setReadable(gse.bp, 'org.Mm.eg.db', 'ENTREZID')

        gse.mf <- gseGO(geneList=gene.list, 
                      ont ="MF",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
       microglia_vs_HM.wt.mf[[id]] <- setReadable(gse.mf, 'org.Mm.eg.db', 'ENTREZID')

      require(ReactomePA)
       microglia_vs_HM.wt.reactome[[id]] <- gsePathway(gene.list, 
                organism = "mouse",
                  pvalueCutoff = 0.1,
                  minGSSize = 3,
                  maxGSSize = 1000,
                  pAdjustMethod = "none", 
                  verbose = TRUE)

     microglia_vs_HM.wt.wiki[[id]] <- gseWP(gene.list,
               organism     = 'Mus musculus',
               
                minGSSize = 3,
                  maxGSSize = 1000,
               pvalueCutoff = 0.1,
               verbose      = FALSE)

    } 
}

