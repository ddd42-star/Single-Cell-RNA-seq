# Set the working directory
setwd("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results")

# run in ther console to speed up computation
#options(future.globals.maxSize = 8000 * 1024^2)

# library needed
library(dplyr)
library(Seurat)
library(patchwork)
library(stringi)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)



# create 12 Seurat's objects
# We already plotted the object without any filter.
cf_mock_33.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/CF_pAEC_Mock_33_count/filtered_feature_bc_matrix.h5")
cf_mock_33 <- CreateSeuratObject(counts = cf_mock_33.data,project = "CF_Mock_33",min.cells = 5, min.features = 1500)
cf_mock_37.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/CF_pAEC_Mock_37_count/filtered_feature_bc_matrix.h5")
cf_mock_37 <- CreateSeuratObject(counts = cf_mock_37.data, project = "CF_Mock_37",min.cells = 5, min.features = 1500)
cf_ns1_33.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/CF_pAEC_NS1_33_count/filtered_feature_bc_matrix.h5")
cf_ns1_33 <- CreateSeuratObject(counts = cf_ns1_33.data, project = "CF_H1N1-R38A_33",min.cells = 5, min.features = 1500)
cf_ns1_37.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/CF_pAEC_NS1_37_count/filtered_feature_bc_matrix.h5")
cf_ns1_37 <- CreateSeuratObject(counts = cf_ns1_37.data, project = "CF_H1N1-R38A_37",min.cells = 5, min.features = 1500)
cf_wt_33.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/CF_pAEC_WT_33_count/filtered_feature_bc_matrix.h5")
cf_wt_33 <- CreateSeuratObject(counts = cf_wt_33.data, project = "CF_H1N1_33",min.cells = 5, min.features = 1500)
cf_wt_37.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/CF_pAEC_WT_37_count/filtered_feature_bc_matrix.h5")
cf_wt_37 <- CreateSeuratObject(counts = cf_wt_37.data, project = "CF_H1N1_37",min.cells = 5, min.features = 1500)

wt_mock_33.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/WT_pAEC_Mock_33_count/filtered_feature_bc_matrix.h5")
wt_mock_33 <- CreateSeuratObject(counts = wt_mock_33.data,project = "WT_Mock_33",min.cells = 5, min.features = 1500)
wt_mock_37.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/WT_pAEC_Mock_37_count/filtered_feature_bc_matrix.h5")
wt_mock_37 <- CreateSeuratObject(counts = wt_mock_37.data, project = "WT_Mock_37",min.cells = 5, min.features = 1500)
wt_ns1_33.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/WT_pAEC_NS1_33_count/filtered_feature_bc_matrix.h5")
wt_ns1_33 <- CreateSeuratObject(counts = wt_ns1_33.data, project = "WT_H1N1-R38A_33",min.cells = 5, min.features = 1500)
wt_ns1_37.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/WT_pAEC_NS1_37_count/filtered_feature_bc_matrix.h5")
wt_ns1_37 <- CreateSeuratObject(counts = wt_ns1_37.data, project = "WT_H1N1-R38A_37",min.cells = 5, min.features = 1500)
wt_wt_33.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/WT_pAEC_WT_33_count/filtered_feature_bc_matrix.h5")
wt_wt_33 <- CreateSeuratObject(counts = wt_wt_33.data, project = "WT_H1N1_33",min.cells = 5, min.features = 1500)
wt_wt_37.data = Read10X_h5("C:/Users/dario/OneDrive/università/MA/research Project/Single-Cell-RNA-seq/results/WT_pAEC_WT_37_count/filtered_feature_bc_matrix.h5")
wt_wt_37 <- CreateSeuratObject(counts = wt_wt_37.data, project = "WT_H1N1_37",min.cells = 5, min.features = 1500)



# concatenate all object together
samples_data <- merge(cf_mock_33,c(cf_mock_37,cf_ns1_33,cf_ns1_37,cf_wt_33,cf_wt_37,wt_mock_33,wt_mock_37, wt_ns1_33, wt_ns1_37, wt_wt_33, wt_wt_37),add.cell.ids=c("CFM33","CFM37","CFN33","CFN37","CFW33","CFW37","WTM33","WTM37","WTN33","WTN37","WTW33","WTW37"), project = "samples")

################################################################################

# Visualize Head of the seurat object
head(colnames(samples_data))

table(samples_data$orig.ident)

# Visualization of the data pre-analysis
# cell counts
# Visualize the number of cell counts per sample
# create dataframe
metadata <- samples_data@meta.data

metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none")+
  ylab("Number of cells")+
  xlab("Samples")


# nCounts_rna (umis) counts per cell
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 10000)

##############################################################################3

# Quality Control before filtering
# Calculate the percentage of RPS,RPL and MT genes
samples_data[["percent.rps"]] <- PercentageFeatureSet(samples_data, pattern = "^Sus-Sscrofa-RPS")

samples_data[["percent.rpl"]] <- PercentageFeatureSet(samples_data, pattern = "^Sus-Sscrofa-RPL")

samples_data[["percent.mt"]] <- PercentageFeatureSet(samples_data, features = c("Sus-Sscrofa-ND1","Sus-Sscrofa-ND2","Sus-Sscrofa-COX1","Sus-Sscrofa-COX2","Sus-Sscrofa-ATP8",
                                                                                "Sus-Sscrofa-ATP6","Sus-Sscrofa-COX3","Sus-Sscrofa-ND3","Sus-Sscrofa-ND4L","Sus-Sscrofa-ND4",
                                                                                "Sus-Sscrofa-ND5","Sus-Sscrofa-ND6","Sus-Sscrofa-CYTB"))
samples_data[["percent.mprl"]] <- PercentageFeatureSet(samples_data, pattern = "^Sus-Sscrofa-MPRL")

# For more information see https://kb.10xgenomics.com/hc/en-us/articles/218169723-What-fraction-of-reads-map-to-ribosomal-proteins

VlnPlot(samples_data, features = "percent.rpl")+NoLegend()
VlnPlot(samples_data, features = "percent.rps")+NoLegend()
VlnPlot(samples_data, features = "percent.mt")+NoLegend()
VlnPlot(samples_data, features = "percent.mprl")+NoLegend() # not present in wt mock in our data

# filter out 2% of ribosomial genes
samples_data <-  subset(samples_data, subset = percent.mt < 20 & percent.rps < 10 & percent.rpl < 10)

# filtered mt
VlnPlot(samples_data, features = "percent.mt")+NoLegend()
VlnPlot(samples_data, features = "percent.rpl")+NoLegend()
VlnPlot(samples_data, features = "percent.rps")+NoLegend()
# nFeature_rna (genes) count per cell
VlnPlot(samples_data, features = "nFeature_RNA")+NoLegend()
# nCounts_rna (umis) counts per cell
VlnPlot(samples_data, features = "nCount_RNA")+NoLegend()


table(samples_data$orig.ident)

#################################################################################

# Secondary Analysis

# normalize
samples_data <- NormalizeData(samples_data)


# feature selection
samples_data <- FindVariableFeatures(samples_data, selection.method = "vst", nfeatures = 3000)

# scaling data regressing out RPL,RPS and MT genes
samples_data <- ScaleData(samples_data, vars.to.regress = c("percent.rpl", "percent.rps", "percent.mt"))

# perform linear dimensionality reduction
samples_data <- RunPCA(samples_data, features = VariableFeatures(object = samples_data))

# Output from PCA

# PC_ 1 
# Positive:  Sus-Sscrofa-IGFBP5, Sus-Sscrofa-F3, Sus-Sscrofa-FNDC3B, Sus-Sscrofa-APOBEC1, Sus-Sscrofa-ARSJ, Sus-Sscrofa-C3, Sus-Sscrofa-RCAN2, Sus-Sscrofa-GASK1B, Sus-Sscrofa-PLAT, Sus-Sscrofa-NECTIN3 
# Sus-Sscrofa-CDH6, Sus-Sscrofa-PADI1, Sus-Sscrofa-EMP1, Sus-Sscrofa-GPRC5A, Sus-Sscrofa-RUNX2, Sus-Sscrofa-KRT4, Sus-Sscrofa-BACH2, Sus-Sscrofa-ADAM28, Sus-Sscrofa-S100A2, Sus-Sscrofa-GPR39 
# Sus-Sscrofa-LAMA3, Sus-Sscrofa-TNIP3, Sus-Sscrofa-ITGB1, Sus-Sscrofa-PAPSS2, Sus-Sscrofa-PALMD, Sus-Sscrofa-CDK6, Sus-Sscrofa-KCNQ5, Sus-Sscrofa-TGFA, Sus-Sscrofa-KRT7, Sus-Sscrofa-BPIFA1 
# Negative:  Sus-Sscrofa-ANKFN1, Sus-Sscrofa-WDR49, Sus-Sscrofa-DNAH9, Sus-Sscrofa-LRRC71, Sus-Sscrofa-DNAH5, Sus-Sscrofa-DTHD1, Sus-Sscrofa-DZIP1L, Sus-Sscrofa-DNAH7, Sus-Sscrofa-ENSSSCG00000040973, Sus-Sscrofa-CFAP45 
# Sus-Sscrofa-ECT2L, Sus-Sscrofa-CFAP44, Sus-Sscrofa-AK7, Sus-Sscrofa-DAW1, Sus-Sscrofa-IQCA1, Sus-Sscrofa-TMEM232, Sus-Sscrofa-CFAP77, Sus-Sscrofa-LRRC36, Sus-Sscrofa-NEK5, Sus-Sscrofa-RIBC2 
# Sus-Sscrofa-DNAI4, Sus-Sscrofa-TPPP3, Sus-Sscrofa-DNAH3, Sus-Sscrofa-FRMPD2, Sus-Sscrofa-VWA3A, Sus-Sscrofa-DNAI1, Sus-Sscrofa-SPAG6, Sus-Sscrofa-CCDC170, Sus-Sscrofa-PACRG, Sus-Sscrofa-CIMIP1 
# PC_ 2 
# Positive:  Sus-Sscrofa-CHST9, Sus-Sscrofa-KRT4, Sus-Sscrofa-CYP26B1, Sus-Sscrofa-DNAH11, Sus-Sscrofa-ELOVL5, Sus-Sscrofa-FMO3, Sus-Sscrofa-KCNMA1, Sus-Sscrofa-NPR3, Sus-Sscrofa-BPIFA1, Sus-Sscrofa-ENPP3 
# Sus-Sscrofa-EPHX1, Sus-Sscrofa-WDPCP, Sus-Sscrofa-ZBTB16, Sus-Sscrofa-FMO4, Sus-Sscrofa-SOX5, Sus-Sscrofa-AKR1C8, Sus-Sscrofa-FMO2, Sus-Sscrofa-PDGFC, Sus-Sscrofa-HSD17B2, Sus-Sscrofa-WFDC2 
# Sus-Sscrofa-ENSSSCG00000015876, Sus-Sscrofa-GALNT18, Sus-Sscrofa-IGFBP5, Sus-Sscrofa-ERBB4, Sus-Sscrofa-PAPSS2, Sus-Sscrofa-ZNF385D, Sus-Sscrofa-PDE4D, Sus-Sscrofa-ENSSSCG00000062186, Sus-Sscrofa-ZNF146, Sus-Sscrofa-TCP11L2 
# Negative:  Sus-Sscrofa-MMP13, Sus-Sscrofa-ALDH1A3, Sus-Sscrofa-SPP1, Sus-Sscrofa-PLAU, Sus-Sscrofa-ERC2, Sus-Sscrofa-FGD6, Sus-Sscrofa-FLNA, Sus-Sscrofa-PMEPA1, Sus-Sscrofa-MAP1B, Sus-Sscrofa-PDPN 
# Sus-Sscrofa-SEMA7A, Sus-Sscrofa-FN1, Sus-Sscrofa-ENSSSCG00000039222, Sus-Sscrofa-FERMT1, Sus-Sscrofa-NWD2, Sus-Sscrofa-TAGLN, Sus-Sscrofa-TUBB2B, Sus-Sscrofa-ODAPH, Sus-Sscrofa-MMP7, Sus-Sscrofa-UPP1 
# Sus-Sscrofa-TNFRSF12A, Sus-Sscrofa-ENSSSCG00000008115, Sus-Sscrofa-SLC12A2, Sus-Sscrofa-SLC1A1, Sus-Sscrofa-ADAMTS17, Sus-Sscrofa-ENSSSCG00000038616, Sus-Sscrofa-TTC9, Sus-Sscrofa-ENSSSCG00000015664, Sus-Sscrofa-NRP2, Sus-Sscrofa-ADAM8 
# PC_ 3 
# Positive:  Sus-Sscrofa-DLK2, Sus-Sscrofa-COL15A1, Sus-Sscrofa-COL17A1, Sus-Sscrofa-AMOTL1, Sus-Sscrofa-TP63, Sus-Sscrofa-KRT5, Sus-Sscrofa-SEMA5A, Sus-Sscrofa-COL18A1, Sus-Sscrofa-SOX6, Sus-Sscrofa-CCDC3 
# Sus-Sscrofa-DENND2C, Sus-Sscrofa-PALLD, Sus-Sscrofa-MCC, Sus-Sscrofa-KLHL29, Sus-Sscrofa-HK2, Sus-Sscrofa-PTPRS, Sus-Sscrofa-ENSSSCG00000012006, Sus-Sscrofa-S100A10, Sus-Sscrofa-GJA1, Sus-Sscrofa-NECTIN1 
# Sus-Sscrofa-NRG1, Sus-Sscrofa-CLDN1, Sus-Sscrofa-VSNL1, Sus-Sscrofa-KRT15, Sus-Sscrofa-SRGAP3, Sus-Sscrofa-BCAM, Sus-Sscrofa-GPC1, Sus-Sscrofa-NFKBIA, Sus-Sscrofa-ENSSSCG00000003521, Sus-Sscrofa-PHYHIP 
# Negative:  Sus-Sscrofa-WFDC2, Sus-Sscrofa-PTI, Sus-Sscrofa-ENSSSCG00000033382, Sus-Sscrofa-RHOB, Sus-Sscrofa-TFF3, Sus-Sscrofa-ENPP3, Sus-Sscrofa-CBR2, Sus-Sscrofa-BPIFA1, Sus-Sscrofa-FMO2, Sus-Sscrofa-CP 
# Sus-Sscrofa-SLC5A8, Sus-Sscrofa-CD24, Sus-Sscrofa-FAM3B, Sus-Sscrofa-PADI1, Sus-Sscrofa-DNAH11, Sus-Sscrofa-ENSSSCG00000058854, Sus-Sscrofa-ENSSSCG00000035053, Sus-Sscrofa-PAM, Sus-Sscrofa-LYPD2, Sus-Sscrofa-SLPI 
# Sus-Sscrofa-ALDH1A1, Sus-Sscrofa-CAPRIN1, Sus-Sscrofa-TMSB4X, Sus-Sscrofa-TMC5, Sus-Sscrofa-ENSSSCG00000006588, Sus-Sscrofa-KRT7, Sus-Sscrofa-ENSSSCG00000017754, Sus-Sscrofa-CYP26B1, Sus-Sscrofa-CAPN13, Sus-Sscrofa-CHST4 
# PC_ 4 
# Positive:  Sus-Sscrofa-CAPN13, Sus-Sscrofa-PRR5L, Sus-Sscrofa-FAP, Sus-Sscrofa-ARHGEF38, Sus-Sscrofa-CFTR, Sus-Sscrofa-ENSSSCG00000033382, Sus-Sscrofa-ATRNL1, Sus-Sscrofa-TNIK, Sus-Sscrofa-SLC5A8, Sus-Sscrofa-CNGA1 
# Sus-Sscrofa-PIP5K1B, Sus-Sscrofa-PROM1, Sus-Sscrofa-KYNU, Sus-Sscrofa-PDE8B, Sus-Sscrofa-DNAH11, Sus-Sscrofa-FNDC3B, Sus-Sscrofa-ITGB3, Sus-Sscrofa-CREB5, Sus-Sscrofa-ADGRA3, Sus-Sscrofa-GOLM1 
# Sus-Sscrofa-TACR1, Sus-Sscrofa-RCAN2, Sus-Sscrofa-SIDT1, Sus-Sscrofa-CYP2C42, Sus-Sscrofa-PLEKHG1, Sus-Sscrofa-ZNF804B, Sus-Sscrofa-DPP4, Sus-Sscrofa-SYT1, Sus-Sscrofa-ADCY2, Sus-Sscrofa-PAM 
# Negative:  Sus-Sscrofa-RPLP1, Sus-Sscrofa-ENSSSCG00000004489, Sus-Sscrofa-RPS15, Sus-Sscrofa-KRT19, Sus-Sscrofa-ENSSSCG00000032003, Sus-Sscrofa-ENSSSCG00000033697, Sus-Sscrofa-RPLP0, Sus-Sscrofa-RPS12, Sus-Sscrofa-RPS29, Sus-Sscrofa-ENSSSCG00000033310 
# Sus-Sscrofa-RPS19, Sus-Sscrofa-S100A6, Sus-Sscrofa-RPL27A, Sus-Sscrofa-RPS27A, Sus-Sscrofa-RPS18, Sus-Sscrofa-ENSSSCG00000035904, Sus-Sscrofa-UBB, Sus-Sscrofa-RPS28, Sus-Sscrofa-RPL23, Sus-Sscrofa-RPL11 
# Sus-Sscrofa-ENSSSCG00000027057, Sus-Sscrofa-UBA52, Sus-Sscrofa-RPL35A, Sus-Sscrofa-ENSSSCG00000013889, Sus-Sscrofa-RACK1, Sus-Sscrofa-RPS14, Sus-Sscrofa-PTMA, Sus-Sscrofa-ENSSSCG00000017509, Sus-Sscrofa-RPL37, Sus-Sscrofa-RPS26 
# PC_ 5 
# Positive:  Sus-Sscrofa-CP, Sus-Sscrofa-CFTR, Sus-Sscrofa-RBMS3, Sus-Sscrofa-SLC39A8, Sus-Sscrofa-STEAP4, Sus-Sscrofa-CALCRL, Sus-Sscrofa-SLC4A4, Sus-Sscrofa-TANGO6, Sus-Sscrofa-SPRY1, Sus-Sscrofa-GLIPR1 
# Sus-Sscrofa-CHST4, Sus-Sscrofa-FNDC3B, Sus-Sscrofa-ENSSSCG00000061874, Sus-Sscrofa-TNIK, Sus-Sscrofa-RGS4, Sus-Sscrofa-PARP8, Sus-Sscrofa-FHIT, Sus-Sscrofa-ETV1, Sus-Sscrofa-SLC5A8, Sus-Sscrofa-PLEKHS1 
# Sus-Sscrofa-KCNJ15, Sus-Sscrofa-SULF1, Sus-Sscrofa-ARHGAP24, Sus-Sscrofa-H1-3, Sus-Sscrofa-MAB21L3, Sus-Sscrofa-FYB2, Sus-Sscrofa-FBXL7, Sus-Sscrofa-KCNC2, Sus-Sscrofa-CAPN13, Sus-Sscrofa-TPT1 
# Negative:  Sus-Sscrofa-IFI6, Sus-Sscrofa-LGALS3, Sus-Sscrofa-C5, Sus-Sscrofa-ENSSSCG00000008115, Sus-Sscrofa-GPRC5A, Sus-Sscrofa-USP43, Sus-Sscrofa-EMP1, Sus-Sscrofa-SCEL, Sus-Sscrofa-CCDC14, Sus-Sscrofa-PHLDA2 
# Sus-Sscrofa-ISG15, Sus-Sscrofa-DSG3, Sus-Sscrofa-ENSSSCG00000036669, Sus-Sscrofa-ACOT7, Sus-Sscrofa-TMPRSS2, Sus-Sscrofa-ENSSSCG00000058854, Sus-Sscrofa-IFFO2, Sus-Sscrofa-EHD1, Sus-Sscrofa-ISG12(A), Sus-Sscrofa-CAVIN3 
# Sus-Sscrofa-SIPA1L2, Sus-Sscrofa-FBXL2, Sus-Sscrofa-CRYBG2, Sus-Sscrofa-CYP26B1, Sus-Sscrofa-ARHGEF28, Sus-Sscrofa-GAS7, Sus-Sscrofa-CDC20B, Sus-Sscrofa-EIF2AK2, Sus-Sscrofa-ENSSSCG00000027118, Sus-Sscrofa-S100A6 

# cluster
samples_data <- FindNeighbors(samples_data, dims = 1:30, reduction = "pca")
samples_data <- FindClusters(samples_data, resolution = 0.5)

# using umap or t-sne
samples_data <- RunTSNE(samples_data, dims = 1:30)

DimPlot(samples_data, reduction = "tsne", split.by = "orig.ident", ncol=3)
DimPlot(samples_data, reduction = "tsne", group.by = "orig.ident")

# integrate CCA (try different one if somenthing wrong with the data)
samples_data <- IntegrateLayers(object = samples_data, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                verbose = TRUE)

# re-join layers after integration
samples_data[["RNA"]] <- JoinLayers(samples_data[["RNA"]])

samples_data <- FindNeighbors(samples_data, reduction = "integrated.cca", dims = 1:30)
samples_data <- FindClusters(samples_data, resolution = 0.5) # add cluster name

# run umap
samples_data <- RunTSNE(samples_data, dims = 1:30, reduction = "integrated.cca")

# Visualization
DimPlot(samples_data, reduction = "tsne" ,split.by = "orig.ident", ncol = 3)
DimPlot(samples_data, reduction = "tsne" ,group.by = "orig.ident")

DimPlot(samples_data, reduction = "tsne")

VlnPlot(samples_data, features = "percent.rpl")
VlnPlot(samples_data, features = "percent.rps")
VlnPlot(samples_data, features = "percent.mt")

# find some cells with manual annotation
# cell annotattion
# 4 - basal (KRT5,TP63,KRT14, NGFR,ITGA6,ITGB4,LAMA3,KRT15,BCAM,DLK2)
# 5 - basal(KRT5,TP63,ITGA6,KRT15,BCAM,DLK2)
# 7 - ciliated (FOXJ1,TPPP3,LRRIQ1,DNAH12,DNAH5,TSPAN1,RIBC2,CCDC170)
# 6 - preciliated(LRRIQ1,DNAH12,TSPAN1)
# 8 - cilliated (FOXJ1,TPPP3,SNTN,LRRIQ1,DNAH12,CAPS,TUBB4B,DNAH5,TSPAN1,RIBC2,LRRC23,CCDC170,CCDC146)
# 1 - NA
# 0 - NA
# 2 - basal(PDPN,ITGA6,ITGB4,BCAM)
# 3 - secretory(SCGB3A1,BPIFB1,PIGR,CXCL8)
# 9 - NA
# 10 - NA
# 11 - NA

# Use this to visualize the Featureplot and Vln plot from our data
FeaturePlot(samples_data, features = "Sus-Sscrofa-BPIFB1", pt.size = 0.2)
VlnPlot(samples_data, features = "Sus-Sscrofa-ENSSSCG00000033452", pt.size = 0.2)

View(samples_data@meta.data)

########################################################################
# After finding the cell types, manually insert the annotation inside our Seurat Object

# add new annotation
# Create a named vector of cluster annotations
cluster_annotations <- c(
  "0" = "NA",
  "1" = "NA",
  "2" = "Basal",
  "3" = "Secretory",
  "4" = "Basal",
  "5" = "Basal",
  "6" = "Preciliated",
  "7" = "Ciliated",
  "8" = "Ciliated",
  "9" = "NA",
  "10" = "NA",
  "11" = "NA"
)

# Assign these annotations to a new metadata column
samples_data@meta.data$cell_type <- cluster_annotations[samples_data@meta.data$seurat_clusters]

# show new clusters
DimPlot(samples_data, label = TRUE)
DimPlot(samples_data, reduction = "tsne", group.by = "cell_type", label = TRUE)
DimPlot(samples_data, reduction = "tsne" ,split.by = "cell_type", ncol = 3)
DimPlot(samples_data, group.by = "cell_type")+ggtitle(NULL)
DimPlot(samples_data, reduction = "tsne", group.by = "cell_type",split.by = "orig.ident", ncol = 3)+ggtitle(NULL)


# find associated gene markers

# add in metadata the information of the differents conditions
# add type
samples_data$type <- NA
samples_data$type[which(stri_detect(samples_data@meta.data$orig.ident, regex="^cf"))] <- "CF"
samples_data$type[which(stri_detect(samples_data@meta.data$orig.ident, regex="^wt"))] <- "WT"

# add infection type
samples_data$infection <- NA
samples_data$infection[which(stri_detect(samples_data@meta.data$orig.ident, regex="_mock_"))] <- "MOCK"
samples_data$infection[which(stri_detect(samples_data@meta.data$orig.ident, regex="_ns1_"))] <- "H1S1"
samples_data$infection[which(stri_detect(samples_data@meta.data$orig.ident, regex="_wt_"))] <- "WT"

# add temperature
samples_data$temperature <- NA
samples_data$temperature[which(stri_detect(samples_data@meta.data$orig.ident, regex="_33"))] <- "33"
samples_data$temperature[which(stri_detect(samples_data@meta.data$orig.ident, regex="_37"))] <- "37"

View(samples_data@meta.data)

# add temperature
samples_data$background <- NA
samples_data$background[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_mock"))] <- "CFMock"
samples_data$background[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_ns1"))] <- "CFh1n1"
samples_data$background[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt"))] <- "CFWt"
samples_data$background[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_mock"))] <- "WTMock"
samples_data$background[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1"))] <- "WTh1n1"
samples_data$background[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_wt"))] <- "WTWt"

View(samples_data@meta.data)

# cf
samples_data$backtemp <- NA
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_mock_33"))] <- "CFMock33"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_ns1_33"))] <- "CFh1n133"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt_33"))] <- "CFWt33"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_mock_33"))] <- "WTMock33"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1_33"))] <- "WTh1n133"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_wt_33"))] <- "WTWt33"

samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_mock_37"))] <- "CFMock37"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_ns1_37"))] <- "CFh1n137"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt_37"))] <- "CFWt37"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_mock_37"))] <- "WTMock37"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1_37"))] <- "WTh1n137"
samples_data$backtemp[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_wt_37"))] <- "WTWt37"

View(samples_data@meta.data)

# cf mock vs wt mock
samples_data$status <- NA
samples_data$status[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_mock"))] <- "CFMock"
samples_data$status[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_mock"))] <- "WTMock"
# cf mock 33 vs wt mock 33
samples_data$statusComb <- NA
samples_data$statusComb[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_mock_33"))] <- "CFMock33"
samples_data$statusComb[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_mock_33"))] <- "WTMock33"
# cf mock 37 vs wt mock 37
samples_data$statusComb[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_mock_37"))] <- "CFMock37"
samples_data$statusComb[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_mock_37"))] <- "WTMock37"
# cf wt+h1n1 vs wt wt h1n1
samples_data$comb <- NA
samples_data$comb[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_ns1"))] <- "CFWtNs1"
samples_data$comb[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1"))] <- "WTWtNs1"
samples_data$comb[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt"))] <- "CFWtNs1"
samples_data$comb[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_wt"))] <- "WTWtNs1"
# cf wt vs wt ns1
samples_data$reco <- NA
samples_data$reco[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt"))] <- "CFWt"
samples_data$reco[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1"))] <- "WTNs1"
# cf wt 33 vs wt ns1 33
samples_data$recomb <- NA
samples_data$recomb[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt_33"))] <- "CFWt33"
samples_data$recomb[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1_33"))] <- "WTNs133"
samples_data$recomb[which(stri_detect(samples_data@meta.data$orig.ident, regex="cf_wt_37"))] <- "CFWt37"
samples_data$recomb[which(stri_detect(samples_data@meta.data$orig.ident, regex="wt_ns1_37"))] <- "WTNs137"

View(samples_data@meta.data)

# save data before continuing, since it takes a lot of time with the amount of data we have
saveRDS(samples_data, file = "FinalAnalysis.rds")


# read data
samples_data <- readRDS(file = "FinalAnalysis.rds")


# Find DE in all clusters independent from cell type

markers <- FindAllMarkers(samples_data)
# change label of the gene name
markers$gene <- gsub("^Sus-Sscrofa-", "", rownames(markers))

# save table
write.csv(markers, "markersInAllGene.csv")

# open in case of further analysis

# filter out all genes with p value adjusted smaller than 0.01
filter_markers <- filter(markers, markers$p_val_adj < 0.01)

# to be continued


###################################################
###################################################

# INFECTED VS UNINFECTED CELLS

samples_data[["percent.virus"]] <- PercentageFeatureSet(samples_data, pattern = "^h1n1")

View(samples_data@meta.data)

# plot distribution of virus percentage

VlnPlot(samples_data,features = "percent.virus", group.by = "orig.ident")

metadata <- samples_data@meta.data

# filter metadata
metadata <- filter(metadata, metadata$orig.ident != "WT_Mock_33" & metadata$orig.ident != "WT_Mock_37" & metadata$orig.ident != "CF_Mock_33" & metadata$orig.ident != "CF_Mock_37")

metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.virus, fill= orig.ident)) + 
  geom_density(alpha=0.2)+
  theme_classic() +
  scale_x_sqrt()+
  ylab("Sqrt-scale number of cells")+
  xlab("Sqrt-scale percentage infected cells")+
  #geom_vline(xintercept = 0.001)+ # median kernel estimation
  theme(legend.title=element_blank())

metadata %>% 
  ggplot(aes(color=orig.ident, x=percent.virus, fill= orig.ident)) + 
  geom_density(alpha=0.2)+
  theme_classic() +
  scale_x_sqrt(limits=c(0,0.1))+
  ylab("Sqrt-scale number of cells")+
  xlab("Sqrt-scale percentage infected cells")+
  geom_vline(xintercept = 0.001)+ # median kernel estimation
  theme(legend.title=element_blank())

 
  
  
# after looking at the distribution of percentage of infected cells
# we will set as threashold 0.001: under this value all the cell will be bystanders
# and above this value all the cells will be infected
samples_data$infection_status <- NA
samples_data$infection_status[samples_data$percent.virus < 0.001] <- "Bystanders"
samples_data$infection_status[samples_data$percent.virus >= 0.001] <- "Infected"

# Count the number of infected and not infected cells for each sample
infection_counts_by_sample <- samples_data@meta.data %>%
  group_by(orig.ident, infection_status) %>%
  summarise(cell_count = n())

# Display the result
print(infection_counts_by_sample)

# Create the barplot
ggplot(infection_counts_by_sample, aes(x = orig.ident, y = cell_count, fill = infection_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Infected vs Not Infected Cells per Sample",
       x = "Sample",
       y = "Number of Cells") +
  scale_fill_manual(values = c("Infected" = "red", "Bystanders" = "blue")) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# Calculate the total number of cells per sample
infection_counts_by_sample <- infection_counts_by_sample %>%
  group_by(orig.ident) %>%
  mutate(total_cells = sum(cell_count)) %>%
  ungroup()

# Calculate the percentage of infected and not infected cells
infection_counts_by_sample <- infection_counts_by_sample %>%
  mutate(percentage = (cell_count / total_cells) * 100)

# Create the normalized barplot showing percentages
ggplot(infection_counts_by_sample, aes(x = orig.ident, y = percentage, fill = infection_status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3)+
  labs(title = "Percentage of Infected vs Not Infected Cells per Sample",
       x = "Sample",
       y = "Percentage of Cells (%)") +
  scale_fill_manual(values = c("Infected" = "red", "Bystanders" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# plot which cells type are infected and which are not

DimPlot(samples_data, reduction = "tsne", group.by = "infection_status", split.by = "cell_type", ncol = 3)+ggtitle(NULL)
DimPlot(samples_data, reduction = "tsne", group.by = "infection_status")+ggtitle(NULL)
DimPlot(samples_data, reduction = "tsne", group.by = c("infection_status", "cell_type"))
DimPlot(samples_data, reduction = "tsne", group.by = "cell_type", split.by = "infection_status", ncol = 2)

## just look

samples_data$count_status <- NA

samples_data$count_status <- (samples_data$percent.virus/100)*samples_data$nCount_RNA
View(samples_data@meta.data)

FeatureScatter(samples_data, feature1 = "nCount_RNA", feature2 = "count_status", group.by = "orig.ident")

metadata <- samples_data@meta.data

#metadata_filter <- metadata[metadata$orig.ident=="cf_wt_37",]
metadata_filter <- metadata %>%
  filter(orig.ident %in% c("cf_wt_37", "cf_wt_33","wt_wt_33", "wt_wt_37"))

plot_multi <- function(metadata){
  # Calculate density
  density_data <- data.frame(
    nCount_RNA = density(metadata$nCount_RNA)$x,
    density_scaled = scales::rescale(density(metadata$nCount_RNA)$y, to = c(0, 6000))
  )
  
  metadata %>%
    ggplot(aes(x=nCount_RNA, y=count_status))+
    geom_point()+
    theme_minimal()+
    ylim(0,6000)+
    geom_line(data = density_data, aes(x = nCount_RNA, y = density_scaled), color = "blue")+
    geom_rug(data = metadata, aes(y=count_status), color="orange", sides = "l")
}


samples_data$cellcondition <- NA

samples_data$cellcondition[which(stri_detect(samples_data@meta.data$orig.ident, regex="_mock_"))] <- "Mock"
samples_data$cellcondition[which(stri_detect(samples_data@meta.data$orig.ident, regex="_ns1_"))] <- "H1N1-R38A"
samples_data$cellcondition[which(stri_detect(samples_data@meta.data$orig.ident, regex="_wt_"))] <- "H1N1"


count_infected_cells <- samples_data@meta.data %>%
  group_by(type,cellcondition,cell_type,temperature,infection_status) %>%
  summarise(cell_count = n())
  ungroup()

count_infected_cells <- count_infected_cells %>%
  group_by(type,cellcondition,cell_type, temperature) %>%
  mutate(proportion = (cell_count / sum(cell_count))*100)
  ungroup()

View(count_infected_cells)


ggplot(count_infected_cells, aes(x = cell_type, y = proportion, fill=infection_status)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(type~cellcondition+temperature)+
  theme_minimal() +
  labs(x = "Cell type", y = "Fraction of cells", fill="Infection status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = sprintf("%.0f", proportion)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 2.5)

##################################################################################
##################################################################################

# CELLULAR BACKGROUND

# lets see how the background influence the infection

# we compare CF_mock vs wt_mock
# we compare CF_mock_33/37 vs wt_mock_33/37
# we compare CF_wt_33/37 vs wt_ns1_33/37s
# we compare CF_mock vs CF_WT/CF_NS1
# we compare WT_mock vs CF_WT/NS1
Idents(samples_data) <- "backtemp"
tmp_markers <- FindMarkers(samples_data, ident.1 = "WTWt37", ident.2 = "WTMock37", assay = "RNA", test.use = "wilcox")

head(tmp_markers)
View(tmp_markers)
# For each comparison write a csv table
#write.csv(tmp_markers, "table22.csv")

tmp_markers$diffexpression <- "NO"
tmp_markers$diffexpression[tmp_markers$avg_log2FC > 2.0 & tmp_markers$p_val_adj < 0.01] <- "UP"
tmp_markers$diffexpression[tmp_markers$avg_log2FC < -2.0 & tmp_markers$p_val_adj < 0.01] <- "DOWN"

tmp_markers$names <- rownames(tmp_markers)
tmp_markers$names <- gsub("^Sus-Sscrofa-", "", tmp_markers$names)

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
tmp_markers$delabel <- ifelse(tmp_markers$diffexpression %in% c("UP", "DOWN"), tmp_markers$names, NA)
# Look at the table and choose the gene to show. Not all gene are related to the immune system
genes_to_label <- c("h1n1--------M","h1n1--------NA","h1n1--------NP",
                    "h1n1--------NS1", "h1n1--------PB2","h1n1--------PB1",
                    "h1n1--------PA","h1n1--------HA","HERC5","USP18","ISG15",
                    "IFIT3","CMPK2","RSAD2","CXCL10","CHI3L2","AMCF-II")

mark <- read.csv("table18.csv")


mark$diffexpression <- "NO"
mark$diffexpression[mark$avg_log2FC > 2.0 & mark$p_val_adj < 0.01] <- "UP"
mark$diffexpression[mark$avg_log2FC < -2.0 & mark$p_val_adj < 0.01] <- "DOWN"


data_lab = subset(tmp_markers,names %in% genes_to_label)
custom_col <- ifelse(tmp_markers$diffexpression == "UP", "#bb0c00", 
                            ifelse(tmp_markers$diffexpression == "DOWN", "#00AFBB", "grey"))
custom_col[is.na(custom_col)] <- 'grey'
names(custom_col)[custom_col == '#bb0c00'] <- 'Upregulated'
names(custom_col)[custom_col == '#00AFBB'] <- 'Downregulated'
names(custom_col)[custom_col == 'grey'] <- 'Not significant'

# plot the DE using EnhancedVolcano
EnhancedVolcano(tmp_markers,
                lab = tmp_markers$names,  # Gene labels
                x = 'avg_log2FC',  # Log fold change
                y = 'p_val_adj',  # Adjusted p-value
                xlim = c(-15, 15),  # Adjust x-axis limits
                ylim = c(0, 350),  # Adjust y-axis limits
                title = 'WT H1N1 37 vs WT Mock 37',  # Plot title
                subtitle = '',
                pCutoff = 0.01,  # p-value cutoff
                FCcutoff = 2.0,  # Fold change cutoff
                pointSize = 2.5,  # Point size
                labSize = 4,  # Label size
                selectLab = genes_to_label,  # Genes to highlight with labels
                boxedLabels = T,  # Draw boxes around labels
                drawConnectors = TRUE,  # Draw lines connecting labels
                widthConnectors = 0.5,  # Width of connector lines
                colConnectors = 'black',  # Color of connectors
                colCustom = custom_col,  # Colors for downregulated, non-significant, upregulated
                colAlpha = 1,  # Transparency of points
                legendLabels = c("Downregulated", "Not significant", "Upregulated"),  # Custom legend labels
                legendPosition = "right",  # Position the legend at the bottom
                legendLabSize = 10,
                axisLabSize = 11,
                caption = NULL,
                gridlines.major = F,
                gridlines.minor = F)


#############################
# You can also do it using ggplot
# Uncomment this section
# ggplot(tmp_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
#   geom_vline(xintercept = c(-2.0, 2.0), col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
#   geom_point(size = 2) +
#   scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
#                      labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
#   coord_cartesian(ylim = c(0, 350), xlim = c(-15, 15)) + # since some genes can have minuslog10padj of inf, we set these limits
#   labs(color = 'Severe', #legend_title, 
#        x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
#   scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
#   ggtitle('title')+ # Plot title
#   geom_label_repel(max.overlaps = Inf) # To show all labels 


