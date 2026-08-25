#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
library(Seurat)
library(dplyr)
library(sctransform)
library(ggplot2)
library(Matrix)
library(patchwork)
library("Seurat")
library('dplyr')
library('ggplot2')
library('MAST')
library('clustree')
library('patchwork')
library('AUCell')
library(readxl)
library(DESeq2)
library(tibble)
library(apeglm)
library(glmGamPoi)
library(abdiv)
# library(DoubletFinder)
library(scRepertoire)
library(stringr)

### 1. Setting working directory and loading data----------------------------------------------------------------
setwd("~/PhD/Fetal 10X/sample_filtered_feature_bc_matrix_LN/sample_filtered_feature_bc_matrix")


# #### not necessary, but to get an idea of the file content
# # Read in matrix.mtx
# counts <- readMM("matrix.mtx.gz")
# 
# # Read in genes.tsv
# genes <- read.table("features.tsv.gz")
# gene_ids <- genes$V2
# 
# # Read in barcodes.tsv
# cells <- read.table("barcodes.tsv.gz")
# cell_ids <- cells$V1
# 
# # Make the column names as the cell IDs and the row names as the gene IDs
# rownames(counts) <- gene_ids
# colnames(counts) <- cell_ids
# count_matrix <- as.data.frame(counts)
# 
# View(count_matrix[1:10,1:10])
# 
# rm(list = ls())
# gc()

#### actual command
FLN_data <- Read10X(data.dir = getwd()) 

### 2. Seurat object  -------------------------------------------------------------------------------------
## Initialize the Seurat object with the raw (non-normalized) data
FLN <- CreateSeuratObject(counts = FLN_data$`Gene Expression`, project = "FLN")
FLN

##add antibody data to new assay (called ADT) https://satijalab.org/seurat/articles/multimodal_vignette
View(FLN_seurat_G@assays[["ADT"]]@data)

FLN_ADT <- CreateAssayObject(counts=FLN_data$`Antibody Capture`)
FLN[['ADT']] <- FLN_ADT

rownames(FLN[['ADT']])
rownames(FLN[['RNA']])

rm(FLN_data)
rm(FLN_ADT)

### 3. Add info to metadata ----------------------------------------------------------------------------------
setwd("~/PhD/Fetal 10X/UsedObjects")
seurat_G <- readRDS("all_ab_old_harmony_standarized.rds")

colnames(seurat_G@meta.data)
rownames(seurat_G@meta.data)

## extract only the FLN data and remove FLN from rownames
FLN_seurat_G <- subset(seurat_G, subset=orig.ident=='FLN_ab') 
rownames(FLN_seurat_G@meta.data) <- str_remove(rownames(FLN_seurat_G@meta.data), 'FLN_')
rownames(FLN_seurat_G@meta.data)

# donor info = 'donor_id_test' column
FLN@meta.data <- cbind(FLN@meta.data, 
                       FLN_seurat_G@meta.data[match(rownames(FLN@meta.data),
                                                    rownames(FLN_seurat_G@meta.data)),
                                              c(12)])
colnames(FLN@meta.data)[6] <- 'Donor'


### 4. Quality control in Seurat--------------------------------------------------------------------------------------
## Change working directory where you want to save the QC images per plate
setwd("~/")
DefaultAssay(FLN) <- 'RNA'

## Pre-processing workflow
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
FLN <- PercentageFeatureSet(FLN, pattern = "^MT-", col.name = 'percent.mt')
FLN <- PercentageFeatureSet(FLN, pattern = '^RP', col.name = 'percent.ribo')

length(rownames(FLN@assays$RNA@counts)[grep('MT', rownames(FLN@assays$RNA@counts))]) #somewhere in name
length(rownames(FLN@assays$RNA@counts)[grep('^MT\\.', rownames(FLN@assays$RNA@counts))]) #start of name followed by .
length(rownames(FLN@assays$RNA@counts)[grep('^MT-', rownames(FLN@assays$RNA@counts))]) #start of name followed by -
length(rownames(FLN@assays$RNA@counts)[grep('^RP', rownames(FLN@assays$RNA@counts))]) #start of name followed by -


# Show QC metrics for the first 5 cells
head(FLN@meta.data, 5)


# Visualize QC metrics as a violin plot
pdf(file="FLN_beforeQC_Vlns.pdf", width = 20, height = 8)
VlnPlot(FLN, features = c("nFeature_RNA", "nCount_RNA", "nCount_ADT","percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf(file="FLN_correlationsQC.pdf", width = 24, height = 8)
plot1 <- FeatureScatter(FLN, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)
plot1.z <- FeatureScatter(FLN, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)+
  scale_x_continuous(limits=c(0,5000))
plot2 <- FeatureScatter(FLN, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.00001)
plot3 <- FeatureScatter(FLN, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.00001)
plot1 + plot1.z+  plot2 + plot3 + plot_layout(ncol=4)
dev.off()

#With filtering cut-offs
pdf(file="FLN_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
plot1 <- ggplot(FLN@meta.data, aes(FLN$nCount_RNA, FLN$nFeature_RNA))+geom_point(size=.5)+theme_bw()+
  geom_vline(xintercept = 32000)+ annotate("text",x=33000,y=1,label=c("32000"),hjust=0, size=2.8)+
  geom_vline(xintercept = 1200)+ annotate("text",x=1100,y=1,label=c("1200"),hjust=0, size=2.8)+
  geom_hline(yintercept = 750)+ annotate("text",x=600,y=1,label=c("750"),hjust=0, size=2.8)
plot2 <- ggplot(FLN@meta.data, aes(FLN$nCount_RNA, FLN$percent.mt))+geom_point(size=.5)+
  theme(axis.text.x = element_text(size=0.2))+theme_bw()+
  geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=15,label=c("1200"),hjust=0, size=2.8)+
  geom_hline(yintercept = 7)+ annotate("text",x=3000,y=8,label=c("7"),vjust=0, size=2.8)+
  xlim(0,5000)
plot3 <- ggplot(FLN@meta.data, aes(FLN$nCount_RNA, FLN$percent.mt))+geom_point(size=.5)+
  theme(axis.text.x = element_text(size=0.2))+theme_bw()+
  geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=20,label=c("1200"),hjust=0, size=2.8)+
  geom_hline(yintercept = 7)+ annotate("text",x=10000,y=9,label=c("7"),vjust=0, size=2.8)
plot4 <- ggplot(FLN@meta.data, aes(FLN$nFeature_RNA, FLN$percent.mt))+geom_point(size=.5)+
  theme(axis.text.x = element_text(size=0.2))+theme_bw()+
  geom_vline(xintercept = 750)+ annotate("text",x=900,y=20,label=c("750"),hjust=0, size=2.8)+
  geom_hline(yintercept = 7)+ annotate("text",x=10000,y=9,label=c("7"),vjust=0, size=2.8)+
  xlim(0,1000)
plot5 <- ggplot(FLN@meta.data, aes(FLN$nCount_RNA, FLN$percent.ribo))+geom_point(size=.5)+theme_bw()+
  geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=50,label=c("1000"),hjust=0, size=2.8)+
  geom_hline(yintercept = 7)+ annotate("text",x=50,y=10,label=c("7"),hjust=0, size=2.8)+xlim(0,5000)
plot6 <- ggplot(FLN@meta.data, aes(FLN$percent.mt, FLN$percent.ribo))+geom_point(size=.5)+theme_bw()
plot1 + plot2 + plot3+plot4+plot5+plot6
dev.off()

## Filtering based on QC metrics
# Check filtering to be set
selected <- WhichCells(FLN, expression = nCount_RNA < 32000 & nFeature_RNA > 750 & 
                         nCount_RNA > 1200 &percent.mt < 7 &percent.ribo >7)
length(selected) #How many cells are left

# Filter cells
FLN_filtered <- subset(FLN, subset = nCount_RNA < 32000 & nFeature_RNA > 750 & 
                         nCount_RNA > 1200 &percent.mt < 7 &percent.ribo >7)


# Visualize QC metrics post-filtering
pdf(file="FLN_afterQC_Vlns.pdf", width = 20, height = 8)
VlnPlot(FLN_filtered, features = c("nFeature_RNA", "nCount_RNA", "nCount_ADT","percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
dev.off()

## remove doublets
FLN_filtered$doublets <- ifelse(FLN_filtered$Donor=='doublet', T, F)
FLN_filtered@meta.data[is.na(FLN_filtered$doublets),'doublets'] <- F
FLN_filtered <- subset(FLN_filtered, subset= doublets==F)

## unassigned donor = NA
FLN_filtered@meta.data[grep('unassigned', FLN_filtered$Donor), 'Donor'] <- NA

## make d004 NA as it is only 1 cell of that donor
FLN_filtered@meta.data[grep('d004', FLN_filtered$Donor), 'Donor'] <- NA

FLN_filtered$Donor <- droplevels(FLN_filtered$Donor)

### 5. Normalization and cell cycle scoring---------------------------------------------------------------------------------------------------------

## Full regression for cell cycle
# Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
FLN_filtered <- NormalizeData(FLN_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign cell cycle scores to the cells 
FLN_filtered <- CellCycleScoring(object = FLN_filtered, s.features = s.genes, g2m.features = g2m.genes, 
                                 set.ident = FALSE)
head(x = FLN_filtered@meta.data)

FLN_filtered <- ScaleData(FLN_filtered)
FLN_filtered <- FindVariableFeatures(FLN_filtered)

## remove TCR genes from variable features 
VariableFeatures(FLN_filtered) <- VariableFeatures(FLN_filtered)[!grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FLN_filtered))]

FLN_filtered <- RunPCA(FLN_filtered)

### 6. ADT data  ------------------------------------------------------------------------------
DefaultAssay(FLN_filtered) <- 'ADT'
FLN_filtered <- NormalizeData(FLN_filtered, normalization.method = 'CLR', margin = 2)

rownames(FLN_filtered@assays$ADT$data)
FLN_filtered$CD4_ADT <- FLN_filtered@assays$ADT$data[6,]
FLN_filtered$CD8_ADT <- FLN_filtered@assays$ADT$data[7,]
FLN_filtered$CD45RA_ADT <- FLN_filtered@assays$ADT$data[3,]

DefaultAssay(FLN_filtered) <- 'RNA'

#### Save and load -----------------------------------------------------------------------------------------------------------

setwd("~/")
save(FLN_filtered, file = 'FLN_filtered_simple.Robj')

## load
setwd("~/")
load('FLN_filtered_simple.Robj')

### 7. Clustering ----------------------------------------------------------------------------------------------------------------------------
DefaultAssay(FLN_filtered) <- 'RNA'
FLN_filtered <- RunUMAP(FLN_filtered, dims = 1:30)

FLN_filtered <- FindNeighbors(object = FLN_filtered, reduction = "pca", dims = 1:30)
FLN_filtered <- FindClusters(object = FLN_filtered, resolution = seq(0.1:1, by=0.1))

DimPlot(object = FLN_filtered, reduction = "umap", pt.size=.2, group.by = 'RNA_snn_res.0.2')

clustree(FLN_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')

FeaturePlot(FLN_filtered, 'CD4.1',reduction='umap', pt.size=.2)
FeaturePlot(FLN_filtered, 'CD8a',reduction='umap', pt.size=.2)
FeaturePlot(FLN_filtered, 'CD45RA',reduction='umap', pt.size=.2)
FeaturePlot(FLN_filtered, 'CD27.1',reduction='umap', pt.size=.2)
DimPlot(FLN_filtered, reduction='umap', group.by = 'Phase')

### 8. Figures ---------------------------------------------------------------------------------------------------
setwd("")
donor_colors <- c(palette.colors(palette = "R4")[2:4],palette.colors(palette = "R4")[6:7])

FLN_filtered@meta.data$FLN_res.0.3_clusters <- as.factor(FLN_filtered@meta.data$RNA_snn_res.0.3)
levels(FLN_filtered@meta.data$FLN_res.0.3_clusters) <- c( '0: CD4 Naive/TCM 1', '1: CD8 Naive-like',
                                                          '2: CD4 Naive/TCM 2','3: Treg',
                                                          '4: CD4 Naive/TCM 3 - stressed',
                                                          '5: Proliferating - S/G2M Phase', '6: CD8 effector/memory')
FLN_filtered@meta.data$FLN_res.0.3_clusters <- factor(FLN_filtered@meta.data$FLN_res.0.3_clusters,
                                                      levels=c(  '5: Proliferating - S/G2M Phase',
                                                                 '1: CD8 Naive-like', '6: CD8 effector/memory',
                                                                 '2: CD4 Naive/TCM 2',
                                                                 '0: CD4 Naive/TCM 1', '4: CD4 Naive/TCM 3 - stressed',
                                                                 '3: Treg'))


FigS1A<- ((FeaturePlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), reduction = "umap", pt.size=.5,
              label.size = 12,
              features = c('CD4.1'), cols = c('gold', 'purple'),
              label = T)+ggtitle('CD4 surface expression')) + 
            (FeaturePlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), reduction = "umap", pt.size=.5,
                        label.size = 12,
                        features = c('CD8a'), cols = c('gold', 'purple'),
                        label = T)+ggtitle('CD8A surface expression'))+
          DimPlot(FLN_filtered, reduction='umap', pt.size=.5, group.by ='Phase', cols=c('grey20', 'violet', 'seagreen1')) + 
          DimPlot(FLN_filtered, reduction='umap', pt.size=.5, group.by ='Donor', cols=donor_colors)  &
          theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
             legend.text = element_text(size=26),legend.title = element_text(size=28),
             legend.key.size = unit(20,'points'),
             plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=2),plot.subtitle = element_text(size=28, hjust=0.5),
             plot.margin = margin(20,20,20,20),
             text=element_text(size=8))) +
          (VlnPlot(object = SetIdent(FLN_filtered, value = "FLN_res.0.3_clusters"),  
                   features = c('CD45RA'), pt.size=0,
                   cols=c('red4','mistyrose', 'plum3',
                          'deepskyblue3','blue', '#669999','turquoise')) +
             ggtitle('CD45RA surface expression')+
             theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=22), 
                   axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
                   legend.text = element_text(size=14),
                   plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),
                   plot.margin = margin(10,10,10,10),
                   text=element_text(size=14))+
             scale_x_discrete(labels=c('5','1','6','2','0','4','3'))+
             xlab('Cluster ID')+
             guides(fill='none'))+
          (ggplot(FLN_filtered@meta.data, aes(x=orig.ident, fill=FLN_res.0.3_clusters)) + theme_classic() +
                geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
                geom_text(stat = 'count', size=6,
                          position = position_fill(vjust=.5), aes(color=FLN_res.0.3_clusters,
                                                                  label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
                labs(fill = "ClusterID", title='Proportion of T cell clusters')+
             theme(axis.text.x = element_blank(), axis.title.x = element_text(size=22), 
                   axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
                   legend.text = element_text(size=24),legend.title = element_text(size=26),
                   plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),plot.subtitle = element_text(size=28, hjust=0.5),
                   plot.margin = margin(40,40,40,40),
                   text=element_text(size=8))+
             scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
                scale_fill_manual(values=c('red4','mistyrose', 'plum3',
                                           'deepskyblue3','blue', '#669999','turquoise'))+
                scale_color_manual(values=c('grey80', rep('grey20',3), 'grey60','grey10', 'grey10'))+
                guides(color='none'))+ 
             plot_layout(ncol=2)+
             plot_annotation(title = "Fetal MLN", 
                             theme=theme(plot.title=element_text(size=50, hjust=0.35, vjust=4), 
                                         plot.margin=margin(60,30,10,30)))
FigS1A 
ggsave(file='FLN_overview_res.0.3_SupplFig1A.pdf', FigS1A, width=20, height=18)

pdf('FLN_SurfaceData_res.0.3.pdf', width = 18, height=18)
FeaturePlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), reduction = "umap", pt.size=.2,
            features =FLN_filtered@assays[["ADT"]]@counts@Dimnames[[1]], 
            label = T, label.size = 14) + plot_layout(ncol=3) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=10))
dev.off()

pdf('FLN_SurfaceData_Vln_res.0.3.pdf', width = 18, height=18)
VlnPlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"),  
        features = FLN_filtered@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
        cols=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3'),
        group.by = 'RNA_snn_res.0.3') + plot_layout(ncol=3)& 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FLN_Clustering_res.0.3.pdf', height=12, width=12)
DimPlot(object = FLN_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.3', 
        label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
        cols=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Lymph Node', color='Cluster ID')+
  scale_color_manual(values=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3'))
dev.off()

pdf('FLN_Clustering_res.0.3_clusternames.pdf', height=12, width=18)
DimPlot(object = FLN_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.3', 
        label = T, label.size = 22,label.color = 'grey50', repel=T,label.box = T,
        cols=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Lymph Node', color='Cluster ID')+
  scale_color_manual(values=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3'),
                     labels= c( '0: CD4 Naive/TCM 1', '0: CD8 Naive-like',
                                '2: CD4 Naive/TCM 2 ','3: Treg',
                                '4: CD4 Naive/TCM 3 - stressed',
                                '5: Proliferating - S/G2M Phase', '6: CD8 effector/memory'))
dev.off()

pdf('FLN_Donor.pdf', height=12, width=12)
DimPlot(object = FLN_filtered, reduction = "umap", pt.size=1, group.by = 'Donor') +
  ggtitle('Donor ID') +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  scale_color_manual(values=palette.colors(palette = "R4")[2:8])
dev.off()

pdf('FLN_QCData_res.0.3.pdf', width = 18, height=10)
FeaturePlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), reduction = "umap", pt.size=.4,
            features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
            label = T, label.size = 10) + DimPlot(FLN_filtered, reduction='umap', pt.size=.4, group.by ='Phase')+
  plot_layout(ncol=3) & 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))
dev.off()

pdf('FLN_QCData_Vln_res.0.3.pdf', width = 14, height=5)
VlnPlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"),  
        features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
        cols=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3')) +
  plot_layout(nrow=1, ncol=4) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FLN_overview_res.0.3.pdf', width = 18, height=16)
FeaturePlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), reduction = "umap", pt.size=.5,
            label.size = 12,
            features = c('CD4.1', 'CD8a', 'CD45RA'),
            label = T) + DimPlot(FLN_filtered, reduction='umap', pt.size=.5, group.by ='Phase') + plot_layout(ncol=2) &
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))
dev.off()

pdf('FLN_CD45RA_res.0.4.pdf', width = 6, height=6)
FeaturePlot(object = SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), reduction = "umap", pt.size=.5,
            label.size = 8,
            features = c('CD45RA'),
            label = T,label.color = c('blue','blue','black','black','blue','turquoise3','red4'))+
  scale_color_gradient(limits=c(0,4),low="lightgrey", high= "blue")+
  ggtitle('CD45RA surface expression')+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20))
dev.off()

pdf('FLN_Clustree_colouredbyCD4ADT.pdf')
clustree(FLN_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')
dev.off()

rownames(FLN_filtered@assays$RNA)[grep('^TRGV|^TRGJ|^TRGD|^TRGC|^TRDV|^TRDJ|^TRDD|^TRDC', rownames(FLN_filtered@assays$RNA))]
pdf('FLN_Vln_TCRgdgenes_res.0.3.pdf', width=14, height=18)
VlnPlot(SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), rownames(FLN_filtered@assays$RNA)[grep('^TRGV|^TRGD|^TRGC|^TRDV|^TRDD|^TRDJ|^TRDC', rownames(FLN_filtered@assays$RNA))], 
        pt.size=0, assay='RNA', slot='data',
        cols=c('black','mistyrose','lightslateblue', 'turquoise','grey','red4','plum3'))&
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FLN_TCRgdgenes_res.0.3.pdf', width=14, height=18)
FeaturePlot(SetIdent(FLN_filtered, value = "RNA_snn_res.0.3"), 
            rownames(FLN_filtered@assays$RNA)[grep('^TRGV|^TRGD|^TRGC|^TRDV|^TRDD|^TRDJ|^TRDC', rownames(FLN_filtered@assays$RNA))], 
            label=T, pt.size=.3, label.size = 4)&
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()


FLN_filtered@meta.data$FLN_res.0.3_clusters <- as.factor(FLN_filtered@meta.data$RNA_snn_res.0.3)
levels(FLN_filtered@meta.data$FLN_res.0.3_clusters) <- c( '0: CD4 Naive/TCM 1', '1: CD8 Naive-like',
                                                          '2: CD4 Naive/TCM 2','3: Treg',
                                                          '4: CD4 Naive/TCM 3 - stressed',
                                                          '5: Proliferating - S/G2M Phase', '6: CD8 effector/memory')
FLN_filtered@meta.data$FLN_res.0.3_clusters <- factor(FLN_filtered@meta.data$FLN_res.0.3_clusters,
                                                      levels=c(  '5: Proliferating - S/G2M Phase',
                                                                 '1: CD8 Naive-like', '6: CD8 effector/memory',
                                                                 '2: CD4 Naive/TCM 2',
                                                                 '0: CD4 Naive/TCM 1', '4: CD4 Naive/TCM 3 - stressed',
                                                                 '3: Treg'))
pdf('FLN_Proportionplot_total_res.0.3.pdf')
ggplot(FLN_filtered@meta.data, aes(x=orig.ident, fill=FLN_res.0.3_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  geom_text(stat = 'count', col='grey40',size=4,
            position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
  labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal Lymph Node')+
  theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('red4','mistyrose','plum3','lightslateblue', 'black','grey','turquoise'))
dev.off()

pdf('FLN_Proportionplot_total_res.0.5_simplecolors.pdf', height=9, width=10)
ggplot(FLN_filtered@meta.data, aes(x=orig.ident, fill=FLN_res.0.3_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  geom_text(stat = 'count', size=5,
            position = position_fill(vjust=.5), aes(color=FLN_res.0.3_clusters,
                                                    label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
  labs(fill = "ClusterID", title='Fetal MLN', subtitle='Proportion of T cell clusters')+
  theme(axis.text.y = element_text(size=16), axis.title.y = element_text(size=18,vjust=3), 
        axis.text.x = element_text(size=20),
        legend.text = element_text(size=16), legend.title = element_text(size=18),
        plot.title = element_text(size=28, hjust=0.5, , vjust=4), plot.subtitle = element_text(size=24, hjust=0.5, , vjust=4),
        plot.margin = margin(30,30,30,30),
        text=element_text(size=16))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('red4','mistyrose', 'plum3',
                             'deepskyblue3','blue', '#669999','turquoise'))+
  scale_color_manual(values=c('grey80', rep('grey20',3), 'grey60','grey10', 'grey10'))+
  guides(color='none')
dev.off()


pdf('FLN_Proportionplot_ClustersperDonor_res.0.3.pdf', width=8,height=5)
ggplot(FLN_filtered@meta.data, aes(x=Donor, fill=FLN_res.0.3_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='', subtitle='')+
  theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('red4','mistyrose','plum3','lightslateblue', 'black','grey','turquoise'))+
  coord_flip()
dev.off()

pdf('FLN_Proportionplot_DonorsperCluster_res.0.3.pdf')
ggplot(FLN_filtered@meta.data, aes(fill=Donor, x=RNA_snn_res.0.3)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  labs(fill = "Donor ID", title='', subtitle='')+
  theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=palette.colors(palette = "R4")[2:8])
dev.off()

pdf('FLN_CD45RA_Vln_res.0.3.pdf', width = 12, height=6)
VlnPlot(object = SetIdent(FLN_filtered, value = "FLN_res.0.3_clusters"),  
        features = c('CD45RA'), pt.size=0,
        cols=c('red4','mistyrose', 'plum3',
               'deepskyblue3','blue', '#669999','turquoise')) +
  ggtitle('CD45RA surface expression')+
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20), legend.text = element_text(size=18),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10))+
  scale_x_discrete(labels=c('5','1','6','2','0','4','3'))
dev.off()

### 9. DEG -----------------------------------------------------------------------------------------------------------
#### ROC/MAST for cluster defining markers 
DefaultAssay(object = FLN_filtered) <- "RNA"
FLN_filtered <- SetIdent(FLN_filtered, value = "RNA_snn_res.0.3")

##Find Markers that are specific for each cluster
Batchedmarkers.mast_i_RNA_data_0.3=FindAllMarkers(FLN_filtered, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_0.3$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.3$pct.1/Batchedmarkers.mast_i_RNA_data_0.3$pct.2

## Create list
listDEgenes_i_RNA_MAST_0.3 <- split(Batchedmarkers.mast_i_RNA_data_0.3, f=Batchedmarkers.mast_i_RNA_data_0.3$cluster)

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
listDEgenes_i_RNA_MAST_0.3 <-lapply(listDEgenes_i_RNA_MAST_0.3, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
listDEgenes_i_RNA_MAST_0.3 <-lapply(listDEgenes_i_RNA_MAST_0.3,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})


##save as Robj
setwd("C:/Users/elise/Documents/PhD/Fetal 10X/DEG")
save(listDEgenes_i_RNA_MAST_0.3, file='SupplData4_listDEG_fetal_FLN_scRNA_0.3.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(listDEgenes_i_RNA_MAST_0.3, file='SupplData4_listDEG_fetal_FLN_scRNA_0.3.xlsx')
detach("package:openxlsx", unload=TRUE)

#### 10. Subset naive and memory CD4 T cells - with Treg ----------------------------------------------------------------------------------
FLN_filtered <- SetIdent(FLN_filtered, value = "RNA_snn_res.0.3")
FLN_CD4 <- subset(FLN_filtered,  idents= c('0','2','3','4'), subset=CD8_ADT<1)
rm(FLN_filtered)

## recluster
FLN_CD4 <- FindVariableFeatures(FLN_CD4)
VariableFeatures(FLN_CD4) <- VariableFeatures(FLN_CD4)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FLN_CD4))]
FLN_CD4 <- ScaleData(FLN_CD4)
FLN_CD4 <- RunPCA(FLN_CD4)
FLN_CD4 <- RunUMAP(FLN_CD4, dims = 1:30)
FLN_CD4 <- FindNeighbors(object = FLN_CD4, reduction = "pca", dims = 1:30)
FLN_CD4 <- FindClusters(object = FLN_CD4, resolution = seq(0.1:1, by=0.1))

DimPlot(object = FLN_CD4, reduction = "umap", pt.size=.5, group.by = 'RNA_snn_res.0.2')+
  DimPlot(object = FLN_CD4, reduction = "umap", pt.size=.5, group.by = 'FLN_res.0.3_clusters')
DimPlot(object = FLN_CD4, reduction = "umap", pt.size=.5, group.by = 'Donor')

clustree(FLN_CD4)

FeaturePlot(FLN_CD4, 'CD4.1',reduction='umap', pt.size=1)
FeaturePlot(FLN_CD4, 'CD8a',reduction='umap', pt.size=1)
FeaturePlot(FLN_CD4, 'CD45RA',reduction='umap', pt.size=1)
FeaturePlot(FLN_CD4, 'CD27.1',reduction='umap', pt.size=.2)
DimPlot(FLN_CD4, reduction='umap', group.by = 'Phase', pt.size=1)

#### Save and load -----------------------------------------------------------------------------------------------------------

setwd("~/")
save(FLN_filtered, file = 'FLN_filtered_simple.Robj')
save(FLN_CD4, file = 'FLN_CD4subset_wTreg.Robj')

## load
setwd("~/")
load('FLN_filtered_simple.Robj')
load('FLN_CD4subset_wTreg.Robj')

#### 10. CD4 - Figures ####
setwd("")

pdf('FLN_CD4_SurfaceData_res.0.2.pdf', width = 18, height=18)
FeaturePlot(object = SetIdent(FLN_CD4, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.2,
            features =FLN_CD4@assays[["ADT"]]@counts@Dimnames[[1]], 
            label = T, label.size = 14) + plot_layout(ncol=3) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=10))
dev.off()

pdf('FLN_CD4_SurfaceData_Vln_res.0.2.pdf', width = 18, height=18)
VlnPlot(object = SetIdent(FLN_CD4, value = "RNA_snn_res.0.2"),  
        features = FLN_CD4@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
        cols=c( 'black','lightslateblue', 'grey','steelblue2'),
        group.by = 'RNA_snn_res.0.2') + plot_layout(ncol=3)& 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FLN_CD4_Clustering_res.0.2.pdf', height=12, width=12)
DimPlot(object = FLN_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.2',
        label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
        cols=c( 'black','lightslateblue', 'grey','steelblue2')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Lymph Node - CD4 subset', color='Cluster ID')
dev.off()

pdf('FLN_CD4_Clustering_res.0.2_clusternames.pdf', height=12, width=18)
DimPlot(object = FLN_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.2', 
        label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
        cols=c( 'black','lightslateblue', 'grey','steelblue2')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Lymph Node - CD4 subset', color='Cluster ID')+
  scale_color_manual(values=c( 'black','lightslateblue', 'grey','steelblue2'),
                     labels=c('0: CD4 Naive/TCM 1', 
                              '1: CD4 Naive/TCM 2 ',
                              '2: CD4 Naive/TCM 3 - stressed',
                              '3: CD4 Naive/TCM 4'))
dev.off()


pdf('FLN_CD4_Oldclusters_totalFLNres0.3.pdf', height=12, width=18)
DimPlot(object = FLN_CD4, reduction = "umap", pt.size=1, group.by = 'FLN_res.0.3_clusters', 
        cols=c('lightslateblue', 'black','grey')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Lymph Node - CD4', color='Cluster ID - total FLN clusters')
dev.off()

pdf('FLN_CD4_Donor.pdf', height=12, width=12)
DimPlot(object = FLN_CD4, reduction = "umap", pt.size=1, group.by = 'Donor') +
  ggtitle('Donor ID') +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  scale_color_manual(values=palette.colors(palette = "R4")[2:8])
dev.off()

pdf('FLN_CD4_QCData_res.0.2.pdf', width = 18, height=10)
FeaturePlot(object = SetIdent(FLN_CD4, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.4,
            features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
            label = T, label.size = 10) + DimPlot(FLN_CD4, reduction='umap', pt.size=.4, group.by ='Phase')+
  plot_layout(ncol=3) & 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))
dev.off()

pdf('FLN_CD4_QCData_Vln_res.0.2.pdf', width = 14, height=5)
VlnPlot(object = SetIdent(FLN_CD4, value = "RNA_snn_res.0.2"),  
        features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
        cols=c('black','lightslateblue', 'grey','steelblue2')) +
  plot_layout(nrow=1, ncol=4) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()


FLN_CD4$CD45RA_ADT <- FLN_CD4@assays$ADT$data[3,]

pdf('FLN_CD4_Clustree_colouredbyCD45RA_ADT.pdf')
clustree(FLN_CD4, node_colour ='CD45RA_ADT', node_colour_aggr = 'median')
dev.off()

FLN_CD4@meta.data$FLN_CD4_res.0.2_clusters <- as.factor(FLN_CD4@meta.data$RNA_snn_res.0.2)
levels(FLN_CD4@meta.data$FLN_CD4_res.0.2_clusters) <- c('0: CD4 Naive/TCM 1', 
                                                        '1: CD4 Naive/TCM 2',
                                                        '2: CD4 Naive/TCM 3 - stressed',
                                                        '3: CD4 Naive/TCM 4')
FLN_CD4@meta.data$FLN_CD4_res.0.2_clusters <- factor(FLN_CD4@meta.data$FLN_CD4_res.0.2_clusters,
                                                     levels=c('0: CD4 Naive/TCM 1', 
                                                              '2: CD4 Naive/TCM 3 - stressed',
                                                              '3: CD4 Naive/TCM 4',
                                                              '1: CD4 Naive/TCM 2'))

pdf('FLN_CD4_Proportionplot_total_res.0.2.pdf')
ggplot(FLN_CD4@meta.data, aes(x=orig.ident, fill=FLN_CD4_res.0.2_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  geom_text(stat = 'count', color='grey30',
            position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
  labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal Lymph Node - CD4 subset')+
  theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('black','grey','steelblue2', 'lightslateblue'))
dev.off()


pdf('FLN_CD4_Proportionplot_ClustersperDonor_res.0.2.pdf', width=8, height=5)
ggplot(FLN_CD4@meta.data, aes(x=Donor, fill=FLN_CD4_res.0.2_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='', subtitle='')+
  theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('black','grey','steelblue2', 'lightslateblue'))+
  coord_flip()
dev.off()

pdf('FLN_CD4_Proportionplot_DonorsperCluster_res.0.2.pdf')
ggplot(FLN_CD4@meta.data, aes(fill=Donor, x=RNA_snn_res.0.2)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  labs(fill = "Donor ID", title='', subtitle='')+
  theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=palette.colors(palette = "R4")[2:8])
dev.off()


#### 11. CD4 - DEG ####
DefaultAssay(object = FLN_CD4) <- "RNA"
FLN_CD4 <- SetIdent(FLN_CD4, value = "RNA_snn_res.0.2")

##Find Markers that are specific for each cluster
Batchedmarkers.mast_i_RNA_data_0.2=FindAllMarkers(FLN_CD4, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_0.2$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.2$pct.1/Batchedmarkers.mast_i_RNA_data_0.2$pct.2

## Create list
listDEgenes_i_RNA_MAST_0.2<- split(Batchedmarkers.mast_i_RNA_data_0.2, f=Batchedmarkers.mast_i_RNA_data_0.2$cluster)

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
listDEgenes_i_RNA_MAST_0.2 <-lapply(listDEgenes_i_RNA_MAST_0.2, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
listDEgenes_i_RNA_MAST_0.2<-lapply(listDEgenes_i_RNA_MAST_0.2,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})




