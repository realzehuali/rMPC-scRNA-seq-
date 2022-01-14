#### KRM cel-cell communication in fig2h ####
# Loading Seurat object as AKI and then split the matrix by days
AKI <- SplitObject(AKI, split.by = "Days")
# Get the DEG and significant LRPs of each day
for (i in c("Day0","Day1","Day3")) {
  data <- GetAssayData(AKI[[i]], slot = "counts", assay = "RNA") %>% as.matrix()
  cluster <- AKI[[i]]@meta.data$clusters %>% as.character()
  Idents(AKI[[i]]) <- "clusters" 
  clust.ana <- cluster_analysis(data = data, genes = rownames(data), cluster = cluster)
  cell_signal <- cell_signaling(data = data,
                                genes = rownames(data),
                                cluster = cluster,
                                int.type = "autocrine", 
                                gene_resive = T,
                                species = 'mus musculus')
  saveRDS(cell_signal, paste("./cluster-analysis/KRM",i,"cell_signal.RDS"))
}

# Establish the adjacency matrix
cell_signal_NC = readRDS("./cluster-analysis/KRM DAY0 cell_signal.RDS")
cell_signal_D1 = readRDS("./cluster-analysis/KRM DAY1 cell_signal.RDS")
cell_signal_D3 = readRDS("./cluster-analysis/KRM DAY3 cell_signal.RDS")

data_NC = data.frame(From = unlist(lapply(names(cell_signal_NC), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                     To = unlist(lapply(names(cell_signal_NC), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                     Value = unlist(lapply(cell_signal_NC, FUN = function(x) {nrow(x)})))
adj_matirx_NC = reshape2::acast(data_NC, From ~ To)

data_D1 = data.frame(From = unlist(lapply(names(cell_signal_D1), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                     To = unlist(lapply(names(cell_signal_D1), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                     Value = unlist(lapply(cell_signal_D1, FUN = function(x) {nrow(x)})))
adj_matirx_D1 = reshape2::acast(data_D1, From ~ To)

data_D3 = data.frame(From = unlist(lapply(names(cell_signal_D3), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                     To = unlist(lapply(names(cell_signal_D3), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                     Value = unlist(lapply(cell_signal_D3, FUN = function(x) {nrow(x)})))
adj_matirx_D3 = reshape2::acast(data_D3, From ~ To)

colnames(adj_matirx_NC) = paste(c(1:4), "NC", sep = "_")
colnames(adj_matirx_D1) = paste(c(1:4), "D1", sep = "_")
colnames(adj_matirx_D3) = paste(c(1:4), "D3", sep = "_")

rownames(adj_matirx_NC) = c(1:4)
rownames(adj_matirx_D1) = c(1:4)
rownames(adj_matirx_D3) = c(1:4)

# Circlize visualization
grid.col = NULL
grid.col[c("cluster 1", "cluster 2", "cluster 3", "cluster 4")] = my_color_palette3[c(1,2,3,4)]
cr = colorRampPalette(c("blue", "white","yellow", "red"))(max(data_NC$Value, data_D1$Value, data_D3$Value)-min(data_NC$Value, data_D1$Value, data_D3$Value)+1)
data_NC$Col = cr[data_NC$Value-min(data_NC$Value, data_D1$Value, data_D3$Value)+1]
data_D1$Col = cr[data_D1$Value-min(data_NC$Value, data_D1$Value, data_D3$Value)+1]
data_D3$Col = cr[data_D3$Value-min(data_NC$Value, data_D1$Value, data_D3$Value)+1]
max = max(data_NC$Value, data_D1$Value, data_D3$Value)
min = min(data_NC$Value, data_D1$Value, data_D3$Value)
data_NC$Value = 1
data_D1$Value = 1
data_D3$Value = 1
data_D1$arrow_col = "lightgrey"
data_D3$arrow_col = "lightgrey"
data_NC$arrow_col = "lightgrey"

pdf("./cell-signaling/new_circos_count_adj_matrix.pdf", width = 9, height = 3)
par(mfrow=c(1,3))
chordDiagram(data_NC,
             grid.col = grid.col,
             col=data_NC$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_NC$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"), 
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_D1,
             grid.col = grid.col,
             col=data_D1$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_D1$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_D3,
             grid.col = grid.col,
             col=data_D3$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_D3$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.05, 0.1)
)
bx = par("usr")
coords = c(bx[1] * 0.95, bx[1] * 0.82, bx[4] * 0.94, 
           bx[4] * 0.53)
l = length(cr)
dy = (coords[3] - coords[4])/l
for (i in 1:l) {
  x = c(coords[1], coords[2], coords[2], coords[1])
  y = c(coords[4] + dy * (i - 1), coords[4] + dy * 
          (i - 1), coords[4] + dy * i, coords[4] + dy * 
          i)
  polygon(x, y, col = cr[i], border = cr[i])
}
text(coords[1] * 0.925, coords[3] * 0.97, labels = as.character(max), 
     col = "black")
text(coords[1] * 0.92, coords[4] * 1.05, labels = as.character(min), 
     col = "black")
text((coords[1] + coords[2])/2, coords[3] * 1.03, labels = "counts")
circos.clear() 
dev.off()

#### KRM cel-cell communication in fig2j ####
# First calculated K_KRM to B_IM(B_monocyte) CCI/CCC
# Loading Seurat object as AKI and then prepare the label
Idents(AKI) = "clusters" # The same as harm.75.1.2
AKI = subset(AKI, idents = c(1,2,3,4,6,7,8,9,11,12))
AKI@meta.data$filteration = paste(AKI@meta.data$Organ, AKI@meta.data$clusters, sep = "_")
Idents(AKI) = "filteration"
AKI = subset(AKI, idents = c("Kidney_1","Kidney_2","Kidney_3","Kidney_4","Blood_6","Blood_7","Blood_8","Blood_9","Blood_11","Blood_12"))
AKI@meta.data = AKI@meta.data %>% mutate(typo = ifelse(grepl("^6|^7|^8|^9|^11|^12", clusters),  "1", "2"), cell = rownames(AKI@meta.data))
rownames(AKI@meta.data) = AKI@meta.data$cell
AKI = SplitObject(AKI, split.by = "Days")

for (i in c("Day0","Day1","Day3")) {
  data = GetAssayData(AKI[[i]], slot = "counts", assay = "RNA") %>% as.matrix()
  cluster = AKI[[i]]@meta.data$typo %>% as.character()
  Idents(AKI[[i]]) = "typo" 
  # clust.ana <- cluster_analysis(data = data, genes = rownames(data), cluster = cluster)
  cell_signal <- cell_signaling(data = data,
                                genes = rownames(data),
                                cluster = cluster,
                                int.type = "autocrine", 
                                gene_resive = T,
                                species = 'mus musculus')
  saveRDS(cell_signal, paste("./cluster-analysis/K_KRM B_IM", i, "cell_signal.RDS"))
}

# Then calculated K_KRM to K_IM CCI/CCC
# Loading Seurat object as AKI and then prepare the label
Idents(AKI) = "clusters" # The same as harm.75.1.2
AKI = subset(AKI, idents = c(1,2,3,4,6,7,8,9,11,12))
AKI@meta.data$filteration = paste(AKI@meta.data$Organ, AKI@meta.data$clusters, sep = "_")
Idents(AKI) = "filteration"
AKI = subset(AKI, idents = c("Kidney_1","Kidney_2","Kidney_3","Kidney_4","Kidney_6","Kidney_7","Kidney_8","Kidney_9","Kidney_11","Kidney_12"))
AKI@meta.data = AKI@meta.data %>% mutate(typo = ifelse(grepl("^6|^7|^8|^9|^11|^12", clusters),  "1", "2"), cell = rownames(AKI@meta.data))
rownames(AKI@meta.data) = AKI@meta.data$cell
AKI = SplitObject(AKI, split.by = "Days")

for (i in c("Day1","Day3")) {
  data = GetAssayData(AKI[[i]], slot = "counts", assay = "RNA") %>% as.matrix()
  cluster = AKI[[i]]@meta.data$typo %>% as.character()
  Idents(AKI[[i]]) = "typo" 
  clust.ana <- cluster_analysis(data = data, genes = rownames(data), cluster = cluster)
  cell_signal <- cell_signaling(data = data,
                                genes = rownames(data),
                                cluster = cluster,
                                int.type = "autocrine", 
                                gene_resive = T,
                                species = 'mus musculus')
  saveRDS(cell_signal, paste("./cluster-analysis/K_KRM K_IM", i, "cell_signal.RDS"))
}
# Merge all 5 LRP group into the same adjacency matrix
cell_signal_NC = readRDS("./cluster-analysis/K_KRM B_IM Day0 cell_signal.RDS")
cell_signal_D1 = readRDS("./cluster-analysis/K_KRM B_IM Day1 cell_signal.RDS")
cell_signal_D3 = readRDS("./cluster-analysis/K_KRM B_IM Day3 cell_signal.RDS")
cell_signal_D1_2 = readRDS("./cluster-analysis/K_KRM K_IM Day1 cell_signal.RDS")
cell_signal_D3_2 = readRDS("./cluster-analysis/K_KRM K_IM Day3 cell_signal.RDS")

data_NC = data.frame(From = unlist(lapply(names(cell_signal_NC), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                     To = unlist(lapply(names(cell_signal_NC), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                     Value = unlist(lapply(cell_signal_NC, FUN = function(x) {nrow(x)})))
adj_matirx_NC = reshape2::acast(data_NC, From ~ To)

data_D1 = data.frame(From = unlist(lapply(names(cell_signal_D1), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                     To = unlist(lapply(names(cell_signal_D1), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                     Value = unlist(lapply(cell_signal_D1, FUN = function(x) {nrow(x)})))
adj_matirx_D1 = reshape2::acast(data_D1, From ~ To)

data_D3 = data.frame(From = unlist(lapply(names(cell_signal_D3), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                     To = unlist(lapply(names(cell_signal_D3), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                     Value = unlist(lapply(cell_signal_D3, FUN = function(x) {nrow(x)})))
adj_matirx_D3 = reshape2::acast(data_D3, From ~ To)

data_D1_2 = data.frame(From = unlist(lapply(names(cell_signal_D1_2), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                       To = unlist(lapply(names(cell_signal_D1_2), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                       Value = unlist(lapply(cell_signal_D1_2, FUN = function(x) {nrow(x)})))
adj_matirx_D1_2 = reshape2::acast(data_D1_2, From ~ To)

data_D3_2 = data.frame(From = unlist(lapply(names(cell_signal_D3_2), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][1])})),
                       To = unlist(lapply(names(cell_signal_D3_2), FUN = function(x) {return(strsplit(x, split = "-",fixed=T)[[1]][2])})),
                       Value = unlist(lapply(cell_signal_D3_2, FUN = function(x) {nrow(x)})))
adj_matirx_D3_2 = reshape2::acast(data_D3_2, From ~ To)


colnames(adj_matirx_NC) = paste(c("B_IM","K_KRM"), "NC", sep = "_")
colnames(adj_matirx_D1) = paste(c("B_IM","K_KRM"), "D1", sep = "_")
colnames(adj_matirx_D3) = paste(c("B_IM","K_KRM"), "D3", sep = "_")
colnames(adj_matirx_D1_2) = paste(c("K_IM","K_KRM"), "D1", sep = "_")
colnames(adj_matirx_D3_2) = paste(c("K_IM","K_KRM"), "D3", sep = "_")

rownames(adj_matirx_NC) = c("IM","K_KRM")
rownames(adj_matirx_D1) = c("IM","K_KRM")
rownames(adj_matirx_D3) = c("IM","K_KRM")
rownames(adj_matirx_D1_2) = c("IM","K_KRM")
rownames(adj_matirx_D3_2) = c("IM","K_KRM")

# Circlize visualization
grid.col = NULL
grid.col[c("cluster 1", "cluster 2")] = my_color_palette3[c(6,4)]
cr = colorRampPalette(c("blue", "white","yellow", "red"))(max(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)-min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)+1)
data_NC$Col = cr[data_NC$Value-min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)+1]
data_D1$Col = cr[data_D1$Value-min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)+1]
data_D3$Col = cr[data_D3$Value-min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)+1]
data_D1_2$Col = cr[data_D1_2$Value-min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)+1]
data_D3_2$Col = cr[data_D3_2$Value-min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)+1]
max = max(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)
min = min(data_NC$Value, data_D1$Value, data_D3$Value, data_D1_2$Value, data_D3_2$Value)
data_NC$Value = 1
data_D1$Value = 1
data_D3$Value = 1
data_D1_2$Value = 1
data_D3_2$Value = 1
data_D1$arrow_col = "lightgrey"
data_D3$arrow_col = "lightgrey"
data_NC$arrow_col = "lightgrey"
data_D3_2$arrow_col = "lightgrey"
data_D1_2$arrow_col = "lightgrey"

pdf("./cell-signaling/circos_count_adj_matrix.pdf", width = 9, height = 3)
par(mfrow=c(2,3))
chordDiagram(data_NC,
             grid.col = grid.col,
             col=data_NC$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_NC$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_D1,
             grid.col = grid.col,
             col=data_D1$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_D1$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"), 
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_D3,
             grid.col = grid.col,
             col=data_D3$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_D3$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_NC,
             grid.col = grid.col,
             col=data_NC$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_NC$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"), 
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_D1_2,
             grid.col = grid.col,
             col=data_D1_2$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_D1_2$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"), 
             annotationTrackHeight = c(0.05, 0.1)
)
chordDiagram(data_D3_2,
             grid.col = grid.col,
             col=data_D3_2$Col,
             transparency = 0.2,
             directional = 1,
             direction.type = "arrows",
             link.arr.col = data_D3_2$arrow_col, 
             link.arr.length = 0.2,
             annotationTrack = c("name", "grid"), 
             annotationTrackHeight = c(0.05, 0.1)
)
bx = par("usr")
coords = c(bx[1] * 0.95, bx[1] * 0.82, bx[4] * 0.94, 
           bx[4] * 0.53)
l = length(cr)
dy = (coords[3] - coords[4])/l
for (i in 1:l) {
  x = c(coords[1], coords[2], coords[2], coords[1])
  y = c(coords[4] + dy * (i - 1), coords[4] + dy * 
          (i - 1), coords[4] + dy * i, coords[4] + dy * 
          i)
  polygon(x, y, col = cr[i], border = cr[i])
}
text(coords[1] * 0.925, coords[3] * 0.97, labels = as.character(max), 
     col = "black")
text(coords[1] * 0.92, coords[4] * 1.05, labels = as.character(min), 
     col = "black")
text((coords[1] + coords[2])/2, coords[3] * 1.03, labels = "counts")
circos.clear() 
dev.off()

