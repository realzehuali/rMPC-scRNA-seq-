# Load the datasets
S0 <- Read10X(data.dir = "./Spleen/SP-0/filtered_feature_bc_matrix/")
S1 <- Read10X(data.dir = "D:/20190813-AKI-SP/SP-1_new/filtered_feature_bc_matrix/")
S3 <- Read10X(data.dir = "./Spleen/SP-3_new/filtered_feature_bc_matrix/")

B0 <- Read10X(data.dir = "./Blood/B-0rep/filtered_feature_bc_matrix/")
B1 <- Read10X(data.dir = "./Blood/MB-1/filtered_feature_bc_matrix/")
B3 <- Read10X(data.dir = "./Blood/B-3/filtered_feature_bc_matrix")

K0 <- Read10X(data.dir = "./Kidney/KM-0/mm10/")
K1 <- Read10X(data.dir = "./Kidney/KM-1/filtered_feature_bc_matrix/")
K3 <- Read10X(data.dir = "./Kidney/KM-3/filtered_feature_bc_matrix/")

# Create Seurat Objects
S0 <- CreateSeuratObject(counts = S0, project = "S0", min.cells = 3, min.features = 200)
S1 <- CreateSeuratObject(counts = S1, project = "S1", min.cells = 3, min.features = 200)
S3 <- CreateSeuratObject(counts = S3, project = "S3", min.cells = 3, min.features = 200)

B0 <- CreateSeuratObject(counts = B0, project = "B0", min.cells = 3, min.features = 200)
B1 <- CreateSeuratObject(counts = B1, project = "B1", min.cells = 3, min.features = 200)
B3 <- CreateSeuratObject(counts = B3, project = "B3", min.cells = 3, min.features = 200)

K0 <- CreateSeuratObject(counts = K0, project = "K0", min.cells = 3, min.features = 200)
K1 <- CreateSeuratObject(counts = K1, project = "K1", min.cells = 3, min.features = 200)
K3 <- CreateSeuratObject(counts = K3, project = "K3", min.cells = 3, min.features = 200)

# merge into one Seurat Object
AKI <- merge(x = S0, y = c(S1, S3, B0, B1, B3, K0, K1, K3))

#label preparation
AKI@meta.data$Organ <- AKI@meta.data$orig.ident
AKI@meta.data$Organ <- dplyr::recode(AKI@meta.data$Organ,
                                     "S0" = "Spleen",
                                     "S1" = "Spleen",
                                     "S3" = "Spleen",
                                     "B0" = "Blood",
                                     "B1" = "Blood",
                                     "B3" = "Blood",
                                     "K0" = "Kidney",
                                     "K1" = "Kidney",
                                     "K3" = "Kidney")

AKI@meta.data$Days <- AKI@meta.data$orig.ident
AKI@meta.data$Days <- dplyr::recode(AKI@meta.data$Days,
                                    "S0" = "Day0",
                                    "S1" = "Day1",
                                    "S3" = "Day3",
                                    "B0" = "Day0",
                                    "B1" = "Day1",
                                    "B3" = "Day3",
                                    "K0" = "Day0",
                                    "K1" = "Day1",
                                    "K3" = "Day3")

# ordering the label
my_level <- c("B0","B1","B3","K0","K1","K3","S0","S1","S3")
AKI@meta.data$orig.ident <- factor(AKI@meta.data$orig.ident, levels = my_level)

my_level <- c("Day0","Day1","Day3","Day5","Day7","Day10","Day17","Day28")
AKI@meta.data$Days <- factor(AKI@meta.data$Days, levels = my_level)

saveRDS(AKI, "AKI ALL RAW.RDS")

# QC:quality control and filtering
AKI[["percent.mt"]] <- PercentageFeatureSet(AKI, pattern = "^mt-")
my_QC <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
color_palette <- c("#00FFFF","#E64B35FF","#00A087FF","#4DBBD5FF","#F39B7FFF","#8491B4FF","#925E9FFF","#FAFD7CFF","#7E6148FF","#B09C85FF","#FFDEAD",
                   "#ED0000FF","#42B540FF","#0099B4FF","#91D1C2FF","#FDAF91FF","#E89242FF","#ADB6B6FF","#FA8072","#3C5488FF","#82491EFF","#526E2DFF", 
                   "#E762D7FF","#AD002AFF","#917C5DFF","#008000","#FF4500","#8A2BE2","#00FF00","#F0E68C","#FF7F50","#1E90FF","#8B008B","#FFB6C1")
g <- list()
for (i in my_QC)
{
  g[[i]] <- VlnPlot(AKI, features = i, group.by = "orig.ident", pt.size = 0.001)  + 
    scale_fill_manual(values = color_palette) +
    geom_boxplot(width = 0.2) + 
    theme_classic() + 
    theme_minimal() +
    labs(title= i, x = "Batch", y = i) +
    stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +
    NoLegend() +
    RotatedAxis() +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) +
    FontSize(x.text = 40,
             y.text = 40,
             x.title = 0,
             y.title = 0,
             main = 40)
}
pdf("QC.pdf")
plot_grid(plotlist = g, ncol = 3)
dev.off()

# filtering cells based on QC
AKI <- subset(AKI, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 50 & nCount_RNA > 1000 & nCount_RNA < 50000)

# SCT pipeline based on Seurat object
AKI <- SCTransform(AKI, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
AKI <- RunPCA(AKI, npcs = 100, verbose = T, assay = "SCT")
AKI <- RunHarmony(AKI, group.by.vars = c("Days","Organ"), assay.use = "SCT", verbose = TRUE, max.iter.harmony = 10000, plot_convergence = T) ## group.by.vars = c("Days","Organ") got more powerful batch effect remove than group.by.vars = c("orig.ident")
AKI <- RunUMAP(AKI, reduction = "harmony", dims = 1:75, reduction.name = "UMAP_HARM", reduction.key = "UMAP_HARM")
ElbowPlot(AKI, ndims = 100, reduction = "harmony") ## we finally kept 75 into downstream analysis
AKI <- RunUMAP(AKI, reduction = "harmony", dims = 1:75, reduction.name = "UMAP_HARM", reduction.key = "UMAP_HARM")
AKI <- RunTSNE(AKI, reduction = "harmony", dims = 1:75, reduction.name = "TSNE_HARM", reduction.key = "TSNE_HARM")
AKI <- FindNeighbors(AKI, reduction = "harmony", dims = 1:75) %>% FindClusters(resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)) ## we finally kept 0.3 into downstream analysis
AKI <- StashAllIdent(object = AKI, ndims = 75, do.harm = T, title = "harm.") ## we finally choose harm.75.0.3
Idents(AKI) <- "harm.75.0.3"
AKI.markers <- FindAllMarkers(AKI, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
AKI.markers.top200 <- AKI.markers %>% group_by(cluster) %>% top_n(200, avg_logFC) ## choosing top200 to take a close look

# MPC subset 
Idents(AKI) <- 'harm.75.0.3'
AKI <- subset(AKI,idents = c(4,5,6,9,10,12,13,14,15))

# MPC subset scRNA-seq data processing
AKI <- clean_my_reduction(AKI) ## cleaning unused data
AKI <- clean_my_assay(AKI, assay = c("RNA")) ## cleaning unused data
AKI <- SCTransform(AKI, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
AKI <- RunPCA(AKI, npcs = 100, verbose = T, assay = "SCT")
AKI <- RunHarmony(AKI, group.by.vars = c("Days","Organ"), assay.use = "SCT", verbose = TRUE, max.iter.harmony = 10000) ## group.by.vars = c("Days","Organ") got more powerful batch effect remove than group.by.vars = c("orig.ident")
ElbowPlot(AKI, ndims = 100, reduction = "harmony") ## we finally kept 75 into downstream analysis
AKI <- RunUMAP(AKI, reduction = "harmony", dims = 1:75, reduction.name = "UMAP_HARM", reduction.key = "UMAP_HARM")
AKI <- RunTSNE(AKI, reduction = "harmony", dims = 1:75, reduction.name = "TSNE_HARM", reduction.key = "TSNE_HARM")
AKI <- FindNeighbors(AKI, reduction = "harmony", dims = 1:75) %>% FindClusters(resolution = c(0.3, 0.4, 0.45, 0.5, 0.6, 0.75, 0.9, 1.0, 1.2)) ## we finally kept 1.2 into downstream analysis
AKI <- StashAllIdent(object = AKI, ndims = 75, do.harm = T, title = "harm.") ## we finally choose harm.75.1.2
Idents(AKI) <- "harm.75.1.2"
AKI.markers <- FindAllMarkers(AKI, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
AKI.markers.top200 <- AKI.markers %>% group_by(cluster) %>% top_n(200, avg_logFC) ## choosing top200 to take a close look

# Lineage origin, organ origin and M1/M2 subtype identification
TRM.genes <- list(c("C1qa","C1qb","C1qc","Cd81","Apoe","H2-Eb1","Ctss","Ms4a7","Cd74","Mmp12"))
Ly6chiIM.genes <- list(c("Chil3","Plac8","Ms4a4c","Lyz1","Ly6c2","Ms4a6c","Ifitm3","Lyz2","Irf7","Gbp2"))
Ly6chilo.genes <- list(c("Ear2","Fabp4","Eno3","Treml4","Pglyrp1","Cebpb","Msrb1","Nr4a1","Gngt2","Smpdl3a"))
Kidney.genes <- list(c("Cd81","Cd74","Cd72","Cd63","Fxyd2","C5ar1","Adgre1","Scimp","Cx3cr1","Cadm1","C1qa","C1qb","C1qc"))
Blood.genes <- list(c("Ly6c2","Ccr2","Ifitm3","Ccr1","Anxa1","Sirpb1c","Fxyd5","Cx3cr1","Cebpb"))
Spleen.genes <- list(c("Ltb","Traf1","Cd83","Itga4","Kird1","Wdfy4","Adam23","Cytip","Ffar2","Fosb","Atf3","Egr1","Fos","Rel","Ccnd1","Klf4","Nr4a1","Purb","Btg2"))
M2a.genes <- list(c("Ccl17","Ccl18","Ccl22","Ccl24","Il1rn","Tgfb1","Tgfb2","Tgfb3","Arg1","Chil3","Chil4","Mgl1","Cd163","Tgm2","Retnla","Ccr2","Mmp9","Mmp12","Igf1","Fn1","Pdgfa","Pdgfb","Pdgfc","Clec7a","Ccl23","Ccl13","Ccl14","Il4r"))
M2b.genes <- list(c("Ccl1","Tnf","Il1b","Il6","Cd74","Il12","Fcgr1"))
M2c.genes <- list(c("Ccl16","Cxcl13","Il10","Mrc1","Slamf1","Stab1","Il10r","Nr3c1"))
M1.genes <- list(c("Fcgr3","Fcgr2","Cd80","Cd86","Il12","Il23","Tnf","Il1b","Il6","Cxcl1","Cxcl2","Cxcl10","Nos2"))
AKI <- CellOriginScoring(AKI, spleen.features = unlist(Spleen.genes), blood.features = unlist(Blood.genes), kidney.features = unlist(Kidney.genes), set.ident = F)
AKI <- CellTypeScoring(AKI, m1.features = unlist(M1.genes), m2a.features = unlist(M2a.genes), m2b.features = unlist(M2b.genes), m2c.features = unlist(M2c.genes), set.ident = F)
AKI <- CellDefineScoring(AKI, trm.features = unlist(TRM.genes), ly6chi.features = unlist(Ly6chiIM.genes), ly6clo.features = unlist(Ly6chilo.genes), set.ident = F) 

# Cellular function identification
apoptosis.genes <- list(c("Bax","Bcl2l11","Bad","Casp3","Casp6","Casp8","Casp9"))
necroptosis.genes <- list(c("Ripk1","Ripk3","Il1b","Il18","Casp8","Casp9","Cflar","Nfkb2","Myd88","Ikbkb","Nfkb1","Txnip","Pycard","Tnfrsf1a","Tnfrsf1b"))
proliferation.genes <- list(c("Bcl2l1","Pcna","Stmn1","Cdkn1a","Hnrnpd","Hnrnpa2b1","Cct2","Cct3","Cct5","Mki67"))
wouldrepair.genes <- list(c("Vegfa","Vegfb","Vegfc","Igf1","Igf2","Fizz1","Arg1","Wnt7b","Egf","Hbegf","Hgf","Pf4","Mrc1","Lyve1","Mertk","Fcrls","Fgf1","Efna1","Efnb","Fgfr2","Ctgf","Icam1","Saa3","Pdgfa","Pdgfb","Pdgfc","Pdgfd","Timp2","Lgals3"))
Chemokines.genes <- list(c("Ccl2","Cxcl2","Ccl3","Ccl4","Cxcl1","Cxcl3","Cxcl10","Ccl7","Ccl8","Cxcl12","Ccl9","Cxcl5","Ccl6","Ccl12","Xcl1","Ccl1","Ccl24","Ccl25","Ccl6","Cxcl9"))
Chemokinesreceptor.genes <- list(c("Ccr2","Ccr1","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccr10","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Cxcr7","Xcr1","Ackr1","Ackr2","Ackr3","Ackr4","Ccrl2"))
Antiinflammatory.genes <- list(c("Arg1","Il4","Il10","Il6","Mertk","Fcgr1","Smad7","Mrc1","Chil3","Ebi3","Il10ra","Il4r"))
inflammatory.genes <- list(c("Mif","Il1b","Il1a","Il6","Il34","Tnf","Nfkb1","Nos2","Irf5","Socs3","Nfkb2","Myd88","Ikbkb"))
Angiogenesis.genes <- list(c("Fgf1","Efna1","Mmp2","Efnb","Epas1","Fgfr2","Ctgf","Col18a1","Vegfa","Vegfb","Col4a1","Fgfr1","Tnfaip2","Igf1","Vcam1","Hbp1","Angpt1","Agtr2","Icam1"))
Phagocytosis.genes <- list(c("P2ry2","Cx3cr1","Gpr132","Lrp1","Cd14","Cd36","Bai1","Timd4","Stab1","Mertk","Axl","Cd209d","Cd163","Nr1h3","Msr1","Cd169","Gulp1","Elmo1","Itgav","Itgb5","Dock1","Rac","Rhog","Lamp1","Lamp2","Rab7","Rab5","Siglec1","Itgam","Itgax","Itga5","Clec7a","Scarf1","C1qa","C1qb","C1qc","Calr","Cd93","Fcgr1","Fcgr2b","Fcgr3"))
trafficking.genes <- list(c("Ccr2","Cx3cr1","Ccr1","Ccr5","Ccr6","Ccr7","Ccr8","Cxcr2","Sell","Selplg","Itgal","Pecam1","Itgam","Itgb2","Itgb1","Itgb4","Icam1","Icam2","Mcam","F11r","Vcam1","Pvr","Cd99","Cd99l2","Sele","Selp","Sell","Epha1","Epha2","Ephb1","Ephb2","Ephb3","Ephb4"))
mature.genes <- list(c("Ctsa","Ctsb","Ctsc","Ctsd","Manf","Manb")) 
RAFG1 = list(c("Spp1","Trem2","Ccl3","Lgals3","Ccl2","Cxcl8","Pdgfa","Pdgfb","Vegfa","Il1b","Ccr2","Fn1","Fcgr1","Tnfaip8l2","Tmco1","Mgat4a","Stab1","Cldn1","Sparcl1","Clec5a","Tlr7","Mmp14","Tnfsf12","Clec11a","Donson", "Apoc2","Cdh1", "Gpnmb", "Mmp12", "Spp1", "Timp2","Arg1")) ## Repair and fibrosis gene set 1
RAFG2 = list(c("Vegfa","Vegfb","Vegfc","Igf1","Igf2","Fizz1","Arg1","Wnt7b","Chi3l3","Il22","Lcn2","Egf","Hbegf","Egfl7","Nrg4","Ace","Gatm","Mertk","Cxcl10","Oasl1","Axl","Tgfa","Atf3","Nr4a1","Il10","Nr3b1","Mmp2","Mmp9","Mmp13","Hgf","Nrf2","Slc7a11","Slpi","Bmp7","Kl","Rln3","Rln1","Lyve1")) ## Repair and fibrosis gene set 2
RAFG3 = list(c("Fcgr1","Tnfaip8l2","Tmco1","Mgat4a","Fn1","Ccr2","Stab1","Cldn1","Sparcl1","Spp1","Snca","Ltc4s","Trem2","Dse","Cdk6","Gal3st4","Clec5a","Tlr7","Defb1","Fabp4","Mek6","Ascl2","Cst6","Vwf","Mmp14","Mef2a","Tnfsf12","Tmem107","Pigs","C3","Hamp","Fcgbp","Apoc2","Clec11a","Donson","Pdxk","Olfml3")) ## Repair and fibrosis gene set 3

AKI <- AddModuleScore(AKI, features = efferocytosis.genes, name = "Efferocytosis")
AKI <- AddModuleScore(AKI, features = ferroptosis.genes, name = "Ferroptosis")
AKI <- AddModuleScore(AKI, features = apoptosis.genes, name = "Apoptosis")
AKI <- AddModuleScore(AKI, features = necroptosis.genes, name = "Necroptosis")
AKI <- AddModuleScore(AKI, features = pyroptosis.genes, name = "Pyroptosis")
AKI <- AddModuleScore(AKI, features = proliferation.genes, name = "Proliferation")
AKI <- AddModuleScore(AKI, features = wouldrepair.genes, name = "Wouldrepair")
AKI <- AddModuleScore(AKI, features = Chemokines.genes, name = "Chemokines")
AKI <- AddModuleScore(AKI, features = Chemokinesreceptor.genes, name = "Chemokinereceptors")
AKI <- AddModuleScore(AKI, features = Antiinflammatory.genes, name = "Antiinflammation")
AKI <- AddModuleScore(AKI, features = inflammatory.genes, name = "Inflammation")
AKI <- AddModuleScore(AKI, features = Angiogenesis.genes, name = "Angiogenesis")
AKI <- AddModuleScore(AKI, features = Phagocytosis.genes, name = "Phagocytosis")
AKI <- AddModuleScore(AKI, features = trafficking.genes, name = "Trafficking")
AKI <- AddModuleScore(AKI, features = RAFG1, name = "RAFG1")
AKI <- AddModuleScore(AKI, features = RAFG2, name = "RAFG2")
AKI <- AddModuleScore(AKI, features = RAFG4, name = "RAFG3")
