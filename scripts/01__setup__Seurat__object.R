source("packages.R")
#----------------------------------------------------------------------------------------------------
samples.paths = list.dirs(path = "/workdir/CITE_SEQ_DATA/good_samples/", full.names = T, recursive = F)

h5.matrices = list()
for (i in 1:length(samples.paths)) {
  h5.matrices[[i]] <- Read10X_h5(list.files(samples.paths, full.names = T)[i], use.names = T)
  names(h5.matrices)[[i]] <- basename(samples.paths)[i]
}

h5.matrices.rna = lapply(h5.matrices, function(x){
  
  x =  x$`Gene Expression`
  
})

#----------------------------------------------------------------------------------------------------
sobjs.list.rna = lapply(h5.matrices.rna, FUN = function(x){
  
  x = CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
})

for (i in 1:length(sobjs.list.rna)) {
  
  sobjs.list.rna[[i]]$orig.ident <- names(sobjs.list.rna)[i]  
}

#----------------------------------------------------------------------------------------------------

for(i in 1:length(sobjs.list.rna)){
  sobjs.list.rna[[i]][["ADT"]] = 
    CreateAssayObject(h5.matrices[[i]][["Antibody Capture"]][, colnames(x = sobjs.list.rna[[i]])])
  
}

# sobj$log10GenesPerUMI <- log10(sobj$nFeature_RNA) / log10(sobj$nCount_RNA)
# sobj.filtered <- subset(sobj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 30  & log10GenesPerUMI > 0.80)

#----------------------------------------------------------------------------------------------------

sobjs.list.rna.filtered = lapply(sobjs.list.rna, function(x){
  AddMetaData(object = x, metadata = log10(x$nFeature_RNA) / log10(x$nCount_RNA), col.name = "log10GenesPerUMI")
})

sobjs.list.rna.filtered = lapply(sobjs.list.rna.filtered, function(x){
  AddMetaData(object = x, metadata = PercentageFeatureSet(x, "^MT-"), col.name = "percent.mt")
})

# sobj.filtered <- subset(sobj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 30  & log10GenesPerUMI > 0.80)

sobjs.list.rna.filtered = lapply(sobjs.list.rna.filtered, function(x){
  x = subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 30  & log10GenesPerUMI > 0.80)
  
})

saveRDS(sobjs.list.rna.filtered, "data/nFeature_RNA_500_nFeature_RNA_5000_percent.mt_30_log10GenesPerUMI_0.80.sobj.list.citeSeq.RDS")

# ---------------------------------- part B ------------------------------------

sobjs.list.rna.filtered <- lapply(X = sobjs.list.rna.filtered, 
                                  FUN = function(x){SCTransform(object = x, method = "glmGamPoi", 
                                 return.only.var.genes = FALSE)})

var.features <- SelectIntegrationFeatures(object.list = sobjs.list.rna.filtered, nfeatures = 3000)

sobj.sct <- merge(x = sobjs.list.rna.filtered[[1]], y = sobjs.list.rna.filtered[2:length(sobjs.list.rna.filtered)], merge.data=TRUE)

VariableFeatures(sobj.sct) <- var.features

sobj.sct.meta = sobj.sct@meta.data
#----------------------------------------------------------------------------------------------------
#save.image(file = "data/dec.15.Rdata")
save.image(file = "data/partB.Rdata")
#----------------------------------------------------------------------------------------------------

write.csv(sobj.sct.meta, file = "results/sct.meta.preIntegration.csv", quote = F)

