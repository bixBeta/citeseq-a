adt.obj <- readRDS("/home/MAIN2/February_03_CITEseq_fixed_AB_annot/data/good_samples__doublet__filtered__sobj.sct___Cite__Seq__RNA__sct__harmony__ndim50__res0.4_0.5_0.6_0.8.RDS")

DefaultAssay(adt.obj) <- "RNA"
DimPlot(adt.obj, label = T, label.size = "7")
adt.obj.meta = adt.obj@meta.data
table(adt.obj.meta$seurat_clusters == adt.obj.meta$SCT_snn_res.0.4)

tcell.clusters = as.factor(c("1","2","3","4","5"))

# select t-cell barcodes
adt.tcell.barcodes.orig.idents = adt.obj.meta %>% 
  rownames_to_column("BC") %>% 
  filter(SCT_snn_res.0.4 %in% tcell.clusters ) %>% 
  select(BC, orig.ident, SCT_snn_res.0.4)

cells.use = unname(unlist(adt.tcell.barcodes.orig.idents %>% select(BC)))

adt.Tcells = subset(adt.obj, cells = cells.use)

saveRDS(adt.Tcells, "data/adt__Tcells__nonC9334__subset.Rds")

DefaultAssay(adt.Tcells) <- "RNA"
DimPlot(adt.Tcells, reduction = "umap", label = T)
