library(tibble)
library(dplyr)
library(networkD3)
library(tidyr)
library(reshape2)
library(viridis)
adt.v.gex = as.matrix(table(adt.obj$adt_snn_res.0.2, adt.obj$SCT_snn_res.0.4))
colnames(adt.v.gex) = paste0("GEX_", colnames(adt.v.gex))
row.names(adt.v.gex) = paste0("ADT_", row.names(adt.v.gex))


getSankey <- function(matrix, floor){
  
  x = as.data.frame(matrix) 
  colnames(x) = c("source", "target", "value")
  
  x = x %>% filter(value > floor)
  
  nodes = as.data.frame(c(as.character(x$source), as.character(x$target)) %>% unique())
  
  colnames(nodes) = "name"
  
  x$IDsource=match(x$source, nodes$name)-1 
  x$IDtarget=match(x$target, nodes$name)-1

  y = sankeyNetwork(Links = x, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=12, nodePadding=20, iterations = 0)
  
  return(y)
  
}



getSankey(matrix = adt.v.gex, floor = 10)


# sankey for comparing Tcell clusters

tcell.clusters = as.factor(c("1","2","3","4","5"))

sobj.sct.tcell.barcodes.orig.idents = sobj.sct.meta %>% 
  rownames_to_column("BC") %>% 
  filter(SCT_snn_res.0.4 %in% tcell.clusters) %>% 
  select(BC, SCT_snn_res.0.4) %>% column_to_rownames("BC") 


megaT.v.reT = as.matrix(table(sobj.sct.tcell.barcodes.orig.idents$SCT_snn_res.0.4, adt.sobj$SCT_snn_res.0.4))
colnames(megaT.v.reT) = paste0("SUB_", colnames(megaT.v.reT))
row.names(megaT.v.reT) = paste0("MEGA_", row.names(megaT.v.reT))

getSankey(matrix = megaT.v.reT, floor = 10)
