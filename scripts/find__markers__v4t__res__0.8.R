DefaultAssay(v4t.obj.added.clusts) <- "RNA"
Idents(v4t.obj.added.clusts) <- "SCT_snn_res.0.8"
v4t.obj.added.clusts$seurat_clusters <- v4t.obj.added.clusts$SCT_snn_res.0.8
v4t.obj.added.clusts$condition = paste0(v4t.obj.added.clusts$phenoGroup, "_", v4t.obj.added.clusts$day)
v4t.obj.added.clusts.meta = v4t.obj.added.clusts@meta.data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD1.controlD1.list = list()

for (i in 1:22) {
  caseD1.controlD1.list[[i]] <- FindMarkers(v4t.obj.added.clusts, ident.1 = "G1_D1", 
                                            ident.2 = "G2_D1", group.by = "condition",
                                            subset.ident = i-1, 
                                            only.pos = F, logfc.threshold = 0)
  
  names(caseD1.controlD1.list)[[i]] <- paste0("CLUSTER__", i-1)
  
}
saveRDS(caseD1.controlD1.list, "SEURAT_OBJECT_V4_T_Added_Clusters/caseD1.controlD1.markers.RDS")



for (i in 1:length(caseD1.controlD1.list)) {
  nde.1[[i]] = length(which(pluck(caseD1.controlD1.list, i, "p_val_adj") < 1))
  names(nde.1)[[i]] <- names(caseD1.controlD1.list)[[i]]
}

nDE.caseD1.controlD1.list = do.call("rbind",nde.1)
colnames(nDE.caseD1.controlD1.list ) = "nDE.caseD1.controlD1.list"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
caseD2.controlD2.list = list()

for (i in 1:22) {
  caseD2.controlD2.list[[i]] <- FindMarkers(v4t.obj.added.clusts, ident.1 = "G1_D2", 
                                            ident.2 = "G2_D2", group.by = "condition",
                                            subset.ident = i-1, 
                                            only.pos = F, logfc.threshold = 0)
  
  names(caseD2.controlD2.list)[[i]] <- paste0("CLUSTER__", i-1)
  
}

save(caseD1.controlD1.list, caseD2.controlD2.list, file= "SEURAT_OBJECT_V3/v4t.obj.added.clusts__DE__case.v.ctrl__top20Clusters.RData")


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


save(caseD1.caseD2.list, controlD1.controlD2.list, caseD1.controlD1.list , caseD2.controlD2.list, file = "SEURAT_OBJECT_V4_T_Added_Clusters/v4t.obj.added.clusts__DE__COMPLETE__case.v.ctrl__.RData")