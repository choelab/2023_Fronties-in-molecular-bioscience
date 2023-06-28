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

###### single-cell RNA inference analysis using monocle3

require(SeuratWrappers)
require(monocle3)

get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]       
  root_pr_nodes
}

cds <- as.cell_data_set(mg.subset)

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
# recreate.partitions

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- mg.subset@reductions$umap@cell.embeddings

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- mg.subset@active.ident#as.factor(mg.subset@meta.data$cell_type_age)
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds <- learn_graph(cds, 
                    close_loop = FALSE, 
                    learn_graph_control = list(minimal_branch_len = 5),
                    use_partition = F)
DN_node_id = get_correct_root_state(cds, cell_phenotype = 'cell_type_age', "MG-02")

cds <- order_cells(cds, reduction_method = "UMAP",  root_pr_nodes = DN_node_id) #root_cells = colnames(cds[, clusters(cds) == "MG-02"])


require(magrittr)
# Get the closest vertice for every cell
y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

mst <- principal_graph(cds)$UMAP

# Get the root vertices
# It is the same node as above
root <- cds@principal_graph_aux$UMAP$root_pr_nodes

# Get the other endpoints
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]


# For each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  # We find the path between the endpoint and the root
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  # We find the cells that map along that path
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()

rownames(cellWeights) <- colnames(cds)
pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
                     nrow = ncol(cds), byrow = FALSE)

cds_sub <-list()
deg <- list()
for(i in 1:length(endpoints)) {

    cds_sub[[i]]<-choose_graph_segments(
              cds,
              reduction_method = "UMAP",
              starting_pr_node = root,
              ending_pr_nodes = endpoints[i],
              return_list = TRUE)
}

for(i in 1:length(endpoints)) {
    cds_subset <- c()
    cds_subset <- cds[,rownames(colData(cds)) %in% cds_sub[[i]]$cells]
    cds_subset$monocle3_pseudotime <- pseudotime(cds_subset)
    fData(cds_subset)$gene_short_name <- rownames(fData(cds_subset))
    data.pseudo <- as.data.frame(colData(cds_subset))

    # ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()

    deg[[i]] <- graph_test(cds_subset, neighbor_graph = "principal_graph")
}

cds_subset <- lapply(cds_sub, function(cds.cell){
    return(cds[,rownames(colData(cds)) %in% cds.cell$cells])
})

###### Trajectory Differential Expression analysis using TradeSeq

counts <- exprs(cds)

require(tradeSeq)
sces <- fitGAM(counts = counts,
              pseudotime = pseudotime,
              cellWeights = cellWeights)


# dynamically expressed genes

assoRes <- associationTest(sces, contrastType = "start",lineages = TRUE, global = FALSE)
patternRes <- patternTest(sces,  global = FALSE,  pairwise = TRUE)
colnames(patternRes)

rownames(patternRes) <- rownames(sces)

gene.p <- data.frame(patternRes) %>% dplyr::filter(pvalue_4vs5 < 0.01) %>% rownames()
gene.4 <- data.frame(assoRes) %>% dplyr::filter(pvalue_4 < 0.01) %>% rownames()
gene.5 <- data.frame(assoRes) %>% dplyr::filter(pvalue_5 < 0.01) %>% rownames()

insGenes <- intersect(gene.p[!grepl("Rik$|^Gm|^mt.|n.", gene.p)],intersect(gene.4, gene.5))

cds_subset.5 <- cds[,rownames(colData(cds)) %in% cds_sub[[5]]$cells]
cds_subset.4 <- cds[,rownames(colData(cds)) %in% cds_sub[[4]]$cells]

pt.matrix.5 <- exprs(cds_subset.5)[match(insGenes,rownames(rowData(cds_subset.5))),order(pseudotime(cds_subset.5))]
pt.matrix.4 <- exprs(cds_subset.4)[match(insGenes,rownames(rowData(cds_subset.4))),order(pseudotime(cds_subset.4))]

pt.matrix.5 <- t(apply(pt.matrix.5, 1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix.5 <- t(apply(pt.matrix.5, 1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix.5) <- insGenes;

pt.matrix.4 <- t(apply(pt.matrix.4, 1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix.4 <- t(apply(pt.matrix.4, 1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix.4) <- insGenes;

###### Figure 2A
require(scales)
Figure2A_total<-DimPlot(brain_cells.subset, 
              # group.by = "age", 
              # cols = c(rev(brewer.pal(10, "Paired")),brewer.pal(9, "Pastel1"),brewer.pal(6, "Set2")),
              label = F, label.size = 5, pt.size = 0.001,raster=FALSE)+
  theme(axis.title = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
        ggtitle(paste0("Total brain - ",comma(length(colnames(brain_cells.subset)),format ="d")," cells")) + NoAxes()

Figure2A_highlight<-DimPlot(brain_cells.subset, 
              sizes.highlight = 0.001,
              cells.highlight = colnames(subset(brain_cells.subset, idents = "13")), 
              cols.highlight = "red", 
              cols = "gray", 
              label = T, label.size = 5, pt.size = 0.001,raster=FALSE)+
  theme(axis.title = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
        ggtitle(paste0("Total brain - ",comma(length(colnames(brain_cells.subset)),format ="d")," cells")) + NoAxes()

Figure2A_WT<-DimPlot(subset(mg.subset, subset = genetype ==names(table(mg.subset$genetype))[2]), 
  group.by = "cell_type_age", 
  label =T, label.size = 5, pt.size = 1,raster=FALSE)+
  theme(axis.title = element_blank(), 
      legend.position = "none",
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("WT - ",comma(length(colnames(subset(mg.subset, subset = genetype ==names(table(mg.subset$genetype))[2]))),format ="d")," cells")) + NoAxes()

Figure2A_5xFAD <- DimPlot(subset(mg.subset, subset = genetype ==names(table(mg.subset$genetype))[1]), 
  group.by = "cell_type_age", 
  label =T, label.size = 5, pt.size = 1,raster=FALSE)+
  theme(axis.title = element_blank(), 
        legend.position = "none",
        plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("5xFAD - ",comma(length(colnames(subset(mg.subset, subset = genetype ==names(table(mg.subset$genetype))[1]))),format ="d")," cells")) + NoAxes()

###### Figure 2B
require(Nebulosa)
Figure2B <- plot_density(mg.subset, 
              features = c("Tmem119","P2ry12","Cx3cr1",
                          "Apoe","Trem2","Cst7", 
                          "Axl","Dkk2","Hif1a"), reduction = "umap" ,
               pal = "inferno", 
               method = c("wkde"), size = 1, adjust=1)

###### Figure 2C
Figure2C_routeA <-plot_cells(cds_subset[[5]],
               color_cells_by = "pseudotime",
               show_trajectory_graph = TRUE,
               trajectory_graph_color = "grey50",
               trajectory_graph_segment_size = 0.25,
               label_cell_groups=FALSE,
               label_roots = TRUE,
               label_leaves=TRUE,
               label_branch_points=FALSE,
               graph_label_size = 2,
               cell_size = 0.5)
Figure2C_routeB <- plot_cells(cds_subset[[4]],
               color_cells_by = "pseudotime",
               show_trajectory_graph = TRUE,
               trajectory_graph_color = "grey50",
               trajectory_graph_segment_size = 0.25,
               label_cell_groups=FALSE,
               label_roots = TRUE,
               label_leaves=TRUE,
               label_branch_points=FALSE,
               graph_label_size = 2,
               cell_size = 0.5)

###### Figure 2D
require(ggvenn)
list_venn <- list("Track A" = gene.5,
                  "Track B" = gene.4)

Figure2D <- ggvenn(list_venn, c("Track A", "Track B"))  

###### Figure 2E
Figure2E <- plotSmoothers(sces, counts, gene = c("Trem2"), lwd = 1, nPoints = 100, lineagesToPlot = c(4,5))

###### Figure 2F
hthc.5 <- Heatmap(
  pt.matrix.5,
  name                         = "z-score",
  # col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  # row_labels = ifelse(rownames(pt.matrix) %in% as.character(convert_human_to_mouse_symbols(unique(unlist(str_split(enrich_result.sel$G,"\\;"))))), rownames(pt.matrix), ""),
  # row_labels = ifelse(rownames(pt.matrix) %in% c("Trem2",gene.4s), rownames(pt.matrix), ""),
  row_names_gp                 = gpar(fontsize = 8),
  km = 5,
  # clustering_method_rows = "ward.D2",
  # clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

hthc.4 <- Heatmap(
  pt.matrix.4,
  name                         = "z-score",
  # col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  # row_labels = ifelse(rownames(pt.matrix.4) %in% as.character(convert_human_to_mouse_symbols(unique(unlist(str_split(enrich_result.sel4$Genes,"\\;"))))), rownames(pt.matrix.4), ""),
  # row_labels = ifelse(rownames(pt.matrix.4) %in% c("Cyld","Tacr1","Tbc1d24","Lrba"), rownames(pt.matrix.4), ""),
  row_names_gp                 = gpar(fontsize = 8),
  km = 4,
  # clustering_method_rows = "ward.D2",
  # clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

Figure2F <-  hthc.5 + hthc.4
