# Note: Much of the below code has been lifted from the Satija Lab Tutorial version 3 2019-10-17
# https://satijalab.org/seurat/immune_alignment.html
# Which has been archive in wayback machine for posterity

library("Seurat")

download.file("https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1","immune_alignment_expression_matrices.zip")
unzip("immune_alignment_expression_matrices.zip")

ctrl.data <- read.table("immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table("immune_stimulated_expression_matrix.txt.gz",  sep = "\t")

###############################################
# Setup the Seurat objects
###############################################

# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

###############################################
# Perform integration
###############################################
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

###############################################
# Perform an integrated analysis
###############################################

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

###############################################
# Identify conserved cell type markers
###############################################

DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
    "CCL2", "PPBP"), min.cutoff = "q9")

immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

DimPlot(immune.combined, label = TRUE)

Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("Mono/Mk Doublets", "pDC", 
    "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", 
    "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
    split.by = "stim") + RotatedAxis()

t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"

dev.off()

###############################################
# extract for each cell type
###############################################
b_cell <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
mk <- FindMarkers(immune.combined, ident.1 = "Mk_STIM", ident.2 = "Mk_CTRL", verbose = FALSE)
pDC <- FindMarkers(immune.combined, ident.1 = "pDC_STIM", ident.2 = "pDC_CTRL", verbose = FALSE)
eryth <- FindMarkers(immune.combined, ident.1 = "Eryth_STIM", ident.2 = "Eryth_CTRL", verbose = FALSE)
DC <- FindMarkers(immune.combined, ident.1 = "DC_STIM", ident.2 = "DC_CTRL", verbose = FALSE)
CD14_mono <- FindMarkers(immune.combined, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
CD16_mono <- FindMarkers(immune.combined, ident.1 = "CD16 Mono_STIM", ident.2 = "CD16 Mono_CTRL", verbose = FALSE)
B_act <- FindMarkers(immune.combined, ident.1 = "B Activated_STIM", ident.2 = "B Activated_CTRL", verbose = FALSE)
CD8_T <-FindMarkers(immune.combined, ident.1 = "CD8 T_STIM", ident.2 = "CD8 T_CTRL", verbose = FALSE)
NK <-FindMarkers(immune.combined, ident.1 = "NK_STIM", ident.2 = "NK_CTRL", verbose = FALSE)
T_act <-FindMarkers(immune.combined, ident.1 = "T activated_STIM", ident.2 = "T activated_CTRL", verbose = FALSE)
CD4_nai<-FindMarkers(immune.combined, ident.1 = "CD4 Naive T_STIM", ident.2 = "CD4 Naive T_CTRL", verbose = FALSE)
CD4_mem<-FindMarkers(immune.combined, ident.1 = "CD4 Memory T_STIM", ident.2 = "CD4 Memory T_CTRL", verbose = FALSE)

tbl<-list("b_cell"=b_cell,
"mk"=mk,
"pDC"=pDC,
"eryth"=eryth,
"DC"=DC,
"CD14_mono"=CD14_mono,
"CD16_mono"=CD16_mono,
"B_act"=B_act,
"CD8_T"=CD8_T,
"NK"=NK,
"T_act"=T_act,
"CD4_nai"=CD4_nai,
"CD4_mem"=CD4_mem)

save.image("seurat_data.RData")

