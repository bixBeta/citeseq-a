v4t <- readRDS("/workdir/DOCKER_DUMP_LM31/data/TCELLS_sct__harmony__ndim50__res0.4_0.5_0.6_0.8.RDS")
cst <- readRDS("/home/MAIN2/February_03_CITEseq_fixed_AB_annot/res_explorer/adt.tcells.object.allRES.RDS")

ref.anchors = FindIntegrationAnchors(object.list = list(v4t, cst), dims = 1:50)

anchors <- FindTransferAnchors(reference = v4t, query = cst,
                                        dims = 1:30, reference.reduction = "pca")

test.query <- MapQuery(anchorset = cst, reference = v4t, query = cst,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")



DefaultAssay(v4t) <- "RNA"
DefaultAssay(cst) <- "RNA"

DotPlot(object = cst, assay="ADT", features=rownames(cst), col.min=-1.5,
        group.by = paste0("SCT_snn_res.", 1.2), cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + 
  RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle(paste0("CITE-seq ADT DotPlot -- res ", 1.2))


sv4markers = c("AHR","FOXO4","CCR10","CTLA4","FOXP3","IL2RA","LAG3","IL10","NCAM1","FCGR3A","KLRK1","NCR3","LILRB1","TRPM3","CD14","LGALS2","S100A8","PF4","MS4A1","CD1C","HLA-DRA","ITGAX","CD34")
DotPlot(object = sv4, assay="RNA", features=c("CXCR3", "PDCD1"), col.min=-1.5,
        group.by = paste0("SCT_snn_res.", 0.6), cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + 
  RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle(paste0("sv4 scRNA-seq  -- res ", 0.6))
