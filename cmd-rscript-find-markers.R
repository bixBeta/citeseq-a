.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")

library("Seurat")
library("glmGamPoi")
library("harmony")
library("patchwork")
library("ggplot2")



v4t.obj.added.clusts <- readRDS("/home/MAIN2/February_03_CITEseq_fixed_AB_annot/res_explorer/v4t.object.allRES.RDS")

DefaultAssay(v4t.obj.added.clusts) <- "RNA"
Idents(v4t.obj.added.clusts) <- "SCT_snn_res.0.8"
v4t.obj.added.clusts$seurat_clusters <- v4t.obj.added.clusts$SCT_snn_res.0.8
v4t.obj.added.clusts$condition = paste0(v4t.obj.added.clusts$phenoGroup, "_", v4t.obj.added.clusts$day)
v4t.obj.added.clusts.meta = v4t.obj.added.clusts@meta.data



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD2.controlD2.list = list()

for (i in 1:22) {
  caseD2.controlD2.list[[i]] <- FindMarkers(v4t.obj.added.clusts, ident.1 = "G1_D2", 
                                            ident.2 = "G2_D2", group.by = "condition",
                                            subset.ident = i-1, 
                                            only.pos = F, logfc.threshold = 0)
  
  names(caseD2.controlD2.list)[[i]] <- paste0("CLUSTER__", i-1)
  
}

#save(caseD1.controlD1.list, caseD2.controlD2.list, file= "SEURAT_OBJECT_V3/v4t.obj.added.clusts__DE__case.v.ctrl__top20Clusters.RData")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD1.caseD2.list = list()

for (i in 1:22) {
  caseD1.caseD2.list[[i]] <- FindMarkers(v4t.obj.added.clusts, ident.1 = "G1_D1", 
                                         ident.2 = "G1_D2", group.by = "condition",
                                         subset.ident = i-1, 
                                         only.pos = F, logfc.threshold = 0)
  
  names(caseD1.caseD2.list)[[i]] <- paste0("CLUSTER__", i-1)
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
controlD1.controlD2.list = list()

for (i in 1:22) {
  controlD1.controlD2.list[[i]] <- FindMarkers(v4t.obj.added.clusts, ident.1 = "G2_D1", 
                                               ident.2 = "G2_D2", group.by = "condition",
                                               subset.ident = i-1, 
                                               only.pos = F, logfc.threshold = 0)
  
  names(controlD1.controlD2.list)[[i]] <- paste0("CLUSTER__", i-1)
  
}


save(caseD1.caseD2.list, controlD1.controlD2.list, caseD2.controlD2.list, file = "SEURAT_OBJECT_V4_T_Added_Clusters/v4t.obj.added.clusts__DE__COMPLETE__case.v.ctrl__.RData")