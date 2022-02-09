library(tibble)
library(dplyr)
library(networkD3)
library(tidyr)
library(reshape2)
library(viridis)
adt.v.gex = as.matrix(table(adt.obj$adt_snn_res.0.2, adt.obj$SCT_snn_res.0.4))
colnames(adt.v.gex) = paste0("GEX_", colnames(adt.v.gex))
row.names(adt.v.gex) = paste0("ADT_", row.names(adt.v.gex))

x = as.data.frame(adt.v.gex)
colnames(x) = c("source", "target", "value")
nodes = as.data.frame(c(as.character(x$source), as.character(x$target)) %>% unique())
colnames(nodes) = "name"

x$IDsource=match(x$source, nodes$name)-1 
x$IDtarget=match(x$target, nodes$name)-1

xf = x %>% filter(value > 10)

sankeyNetwork(Links = xf, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, nodeWidth=40, fontSize=12, nodePadding=20, iterations = 0)



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
