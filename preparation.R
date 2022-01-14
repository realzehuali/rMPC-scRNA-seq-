#### general setting and loading packages required ####
my_color_palette = c("grey94","blueviolet","#00A087FF","red1","#F39B7FFF","#8491B4FF","#925E9FFF","#FAFD7CFF","#7E6148FF","#B09C85FF",
                     "#FFDEAD","lightyellow","#42B540FF","#0099B4FF","green2","#FDAF91FF","#E89242FF","#00FFFF","#FA8072","#3C5488FF",
                     "#82491EFF","springgreen","#E762D7FF","#AD002AFF","#917C5DFF","#008000","blue2","thistle1","lightblue","#F0E68C",
                     "mistyrose","beige","#8B008B")

from = c(2,10,14,17,7,3,4,6,21,26,5,1,9,27,22,15,25,24,28,23,16,20,0,12,31,8,11,29,13,18,19,30)
my_color_palette3 = my_color_palette[from+1]
my_color_palette4 = my_color_palette3[c(1,2,3,4,rep(23,28))]
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(reticulate)
library(viridis)
library(scales)
library(plyr)
library(harmony)
library(SingleCellExperiment)
library(scDblFinder)
library(pheatmap)
library(clusterProfiler)
library(clusterProfiler.dplyr)
library(ReactomePA)
library(reactome.db)
library(enrichplot)
library(msigdbr)
library(msigdf)
library(hrbrthemes)
library(monocle)
library(ggsci)
library(ggsignif)
library(Hmisc)
library(data.table)
library(igraph)
library(circlize)
library(scsrctdb)

#### helper function required ####
#### CellOriginScoring, CellDefineScoring and CellTypeScoring####
CheckGC <- function() {
  if (getOption(x = "Seurat.memsafe")) {
    gc(verbose = FALSE)
  }
}

CellOriginScoring <- function (object, spleen.features, blood.features, kidney.features, set.ident = FALSE, ...) 
{
  name <- "Cell.Origin"
  features <- list(Spleen.Score = spleen.features, Blood.Score = blood.features, Kidney.Score = kidney.features)
  object.co <- AddModuleScore(object = object, features = features, name = name, ctrl = min(vapply(X = features, FUN = length, 
                                                                                                   FUN.VALUE = numeric(length = 1))), ...)
  co.columns <- grep(pattern = name, x = colnames(x = object.co[[]]), value = TRUE)
  co.scores <- object.co[[co.columns]]
  rm(object.co)
  CheckGC()
  assignments <- apply(X = co.scores, MARGIN = 1, FUN = function(scores, first = "Spleen", second = "Blood", third = "Kidney", null = "Undecided") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second, third)[which(x = scores == 
                                               max(scores))])
      }
    }
  })
  co.scores <- merge(x = co.scores, y = data.frame(assignments), by = 0)
  colnames(x = co.scores) <- c("rownames", "Spleen.Score", "Blood.Score", "Kidney.Score", "Origin")
  rownames(x = co.scores) <- co.scores$rownames
  co.scores <- co.scores[, c("Spleen.Score", "Blood.Score", "Kidney.Score", "Origin")]
  object[[colnames(x = co.scores)]] <- co.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Origin"
  }
  return(object)
}

CellDefineScoring <- function (object, trm.features, ly6chi.features, ly6clo.features, set.ident = FALSE, ...) 
{
  name <- "Cell.Define"
  features <- list(TRM.Score = trm.features, Ly6chi.Score = ly6chi.features, Ly6clo.Score = ly6clo.features)
  object.co <- AddModuleScore(object = object, features = features, name = name, ctrl = min(vapply(X = features, FUN = length, 
                                                                                                   FUN.VALUE = numeric(length = 1))), ...)
  co.columns <- grep(pattern = name, x = colnames(x = object.co[[]]), value = TRUE)
  co.scores <- object.co[[co.columns]]
  rm(object.co)
  CheckGC()
  assignments <- apply(X = co.scores, MARGIN = 1, FUN = function(scores, first = "TRM", second = "Ly6chi", third = "Ly6clo", null = "Undecided") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second, third)[which(x = scores == 
                                               max(scores))])
      }
    }
  })
  co.scores <- merge(x = co.scores, y = data.frame(assignments), by = 0)
  colnames(x = co.scores) <- c("rownames", "TRM.Score", "Ly6chi.Score", "Ly6clo.Score", "Define")
  rownames(x = co.scores) <- co.scores$rownames
  co.scores <- co.scores[, c("TRM.Score", "Ly6chi.Score", "Ly6clo.Score", "Define")]
  object[[colnames(x = co.scores)]] <- co.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Define"
  }
  return(object)
}

CellTypeScoring <- function (object, m1.features, m2a.features, m2b.features, m2c.features, set.ident = FALSE, ...) 
{
  name <- "Cell.Type"
  features <- list(M1.Score = m1.features, M2a.Score = m2a.features, M2b.Score = m2b.features, M2c.Score = m2c.features)
  object.ct <- AddModuleScore(object = object, features = features, name = name, ctrl = min(vapply(X = features, FUN = length, 
                                                                                                   FUN.VALUE = numeric(length = 1))), ...)
  ct.columns <- grep(pattern = name, x = colnames(x = object.ct[[]]), value = TRUE)
  ct.scores <- object.ct[[ct.columns]]
  rm(object.ct)
  CheckGC()
  assignments <- apply(X = ct.scores, MARGIN = 1, FUN = function(scores, first = "M1", second = "M2a", third = "M2b", fourth = "M2c", null = "Undecided") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second, third, fourth)[which(x = scores == 
                                                       max(scores))])
      }
    }
  })
  ct.scores <- merge(x = ct.scores, y = data.frame(assignments), by = 0)
  colnames(x = ct.scores) <- c("rownames", "M1.Score", "M2a.Score", "M2b.Score", "M2c.Score", "Type")
  rownames(x = ct.scores) <- ct.scores$rownames
  ct.scores <- ct.scores[, c("M1.Score", "M2a.Score", "M2b.Score", "M2c.Score", "Type")]
  object[[colnames(x = ct.scores)]] <- ct.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Type"
  }
  return(object)
}

#### Seurat data processing, visualization and memory saving ####
clean_my_assay <- function(object = object, assay = c("RNA","SCT")) {
  for (i in Seurat::Assays(object)) {
    if (!(i %in% assay)) {
      object[[i]] <- NULL
    }
  }
  return(object)
  gc()
}

clean_my_reduction <- function(object = object) {
  for (i in Reductions(object)){
    object[[i]] <- NULL
  }
  return(object)
}  

StashAllIdent <- function(object = object, ndims = ndims, do.harm = F, title = "harm."){
  a = colnames(object[[]])[grep("SCT_snn_", colnames(object[[]]))]
  if (do.harm == F) {
    b = gsub(pattern = "SCT_snn_res", replacement = paste0("X", ndims), a)
  } else if (do.harm == T) {
    b = gsub(pattern = "SCT_snn_res", replacement = paste0(title, ndims), a)
  }
  for (i in 1:length(a)) {
    Idents(object = object) <- a[i]
    object <- StashIdent(object, save.name = b[i])
  }
  return(object)
}

set_reduction_assay_name <- function(object, set) {
  for (i in Reductions(object)) {
    object[[i]]@assay.used <- set
  }
  return(object)
}

RotatedAxis <- function (angle = 45, ...) 
{
  rotated.theme <- theme(axis.text.x = element_text(angle = angle, 
                                                    hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}

set_active_seurat_assay = function(object, assay) {
  DefaultAssay(object) = assay
  return(object)
}

set_active_seurat_ident = function(object, ident) {
  Idents(object) = ident
  return(object)
}

set_my_reduction_global <- function(object = object) {
  for (i in Reductions(object)){
    object[[i]]@global <- "T" %>% as.logical()
  }
  return(object)
}  

set_my_meta_char = function(object = object) {
  for (i in colnames(object@meta.data)){
    object[[i]][,1] <- object[[i]][,1] %>% as.character()
  }
  return(object)
}  

#### Seurat object V3 converting to Monocle object V2
seuratimportcds2 <- function(object = object, assay = "SCT", slot = "data", expressionFamily = NULL, doestimate = T) 
{
  message("Extracting data.")
  #Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(GetAssayData(object, assay = assay, slot = slot)), "sparseMatrix")
  
  pd <- new('AnnotatedDataFrame', data = object@meta.data)
  
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  #Construct monocle cds
  message("Constructing monocle cds.")
  monocle_cds <- monocle::newCellDataSet(cellData =data,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.1, 
                                         expressionFamily = VGAM::negbinomial.size())
  if (doestimate) {
    message("estimating data")
    monocle_cds <- estimateSizeFactors(monocle_cds)
    
    if(is.null(expressionFamily)) {
      monocle_cds <- estimateDispersions(monocle_cds)
    }
    
    message("detecting Genes") 
    monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  }
  return(monocle_cds)
}
