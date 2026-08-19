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
library(VDJdive)
library(S4Vectors)
library(purrr)
library(harmony)

### 1. Setting working directory and loading data----------------------------------------------------------------
setwd("~/PhD/Fetal 10X/FetalBlood_files_Guillem/")

#### not necessary, but to get an idea of the file content
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
FBL_data <- Read10X(data.dir = getwd()) 

### 2. Seurat object  -------------------------------------------------------------------------------------
## Initialize the Seurat object with the raw (non-normalized) data
FBL <- CreateSeuratObject(counts = FBL_data$`Gene Expression`, project = "FBL")
FBL

FBL_ADT <- CreateAssayObject(counts=FBL_data$`Antibody Capture`)
FBL[['ADT']] <- FBL_ADT

rownames(FBL[['ADT']])
rownames(FBL[['RNA']])

rm(FBL_data)
rm(FBL_ADT)


### 4. Quality control in Seurat--------------------------------------------------------------------------------------
## Change working directory where you want to save the QC images per plate
setwd("~/")
DefaultAssay(FBL) <- 'RNA'

## Pre-processing workflow
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
FBL <- PercentageFeatureSet(FBL, pattern = "^MT-", col.name = 'percent.mt')
FBL <- PercentageFeatureSet(FBL, pattern = '^RP', col.name = 'percent.ribo')

length(rownames(FBL@assays$RNA)[grep('MT', rownames(FBL@assays$RNA))]) #somewhere in name
length(rownames(FBL@assays$RNA)[grep('^MT\\.', rownames(FBL@assays$RNA))]) #start of name followed by .
length(rownames(FBL@assays$RNA)[grep('^MT-', rownames(FBL@assays$RNA))]) #start of name followed by -
length(rownames(FBL@assays$RNA)[grep('^RP', rownames(FBL@assays$RNA))]) #start of name followed by -


# Show QC metrics for the first 5 cells
head(FBL@meta.data, 5)


# Visualize QC metrics as a violin plot
pdf(file="FBL_beforeQC_Vlns.pdf", width = 20, height = 8)
VlnPlot(FBL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf(file="FBL_correlationsQC.pdf", width = 24, height = 8)
plot1 <- FeatureScatter(FBL, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)
plot1.z <- FeatureScatter(FBL, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)+
  scale_x_continuous(limits=c(0,5000))
plot2 <- FeatureScatter(FBL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.00001)
plot3 <- FeatureScatter(FBL, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.00001)
plot1 + plot1.z+  plot2 + plot3 + plot_layout(ncol=4)
dev.off()

#With filtering cut-offs
pdf(file="FBL_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
plot1 <- ggplot(FBL@meta.data, aes(FBL$nCount_RNA, FBL$nFeature_RNA))+geom_point(size=.5)+theme_bw()+
  geom_vline(xintercept = 30000)+ annotate("text",x=31000,y=1,label=c("30000"),hjust=0, size=2.8)+
  geom_vline(xintercept = 2000)+ annotate("text",x=2100,y=1,label=c("2000"),hjust=0, size=2.8)+
  geom_hline(yintercept = 750)+ annotate("text",x=600,y=1,label=c("750"),hjust=0, size=2.8)
plot2 <- ggplot(FBL@meta.data, aes(FBL$nCount_RNA, FBL$percent.mt))+geom_point(size=.5)+
  theme(axis.text.x = element_text(size=0.3))+theme_bw()+
  geom_vline(xintercept = 2000)+ annotate("text",x=2200,y=15,label=c("2000"),hjust=0, size=2.8)+
  geom_hline(yintercept = 4)+ annotate("text",x=3000,y=8,label=c("4"),vjust=0, size=2.8)+
  xlim(0,5000)
plot3 <- ggplot(FBL@meta.data, aes(FBL$nCount_RNA, FBL$percent.mt))+geom_point(size=.5)+
  theme(axis.text.x = element_text(size=0.3))+theme_bw()+
  geom_vline(xintercept = 2000)+ annotate("text",x=2200,y=20,label=c("2000"),hjust=0, size=2.8)+
  geom_hline(yintercept = 4)+ annotate("text",x=10000,y=9,label=c("4"),vjust=0, size=2.8)
plot4 <- ggplot(FBL@meta.data, aes(FBL$nFeature_RNA, FBL$percent.mt))+geom_point(size=.5)+
  theme(axis.text.x = element_text(size=0.3))+theme_bw()+
  geom_vline(xintercept = 750)+ annotate("text",x=900,y=20,label=c("750"),hjust=0, size=2.8)+
  geom_hline(yintercept =4)+ annotate("text",x=10000,y=9,label=c("4"),vjust=0, size=2.8)+
  xlim(0,1000)
plot5 <- ggplot(FBL@meta.data, aes(FBL$nCount_RNA, FBL$percent.ribo))+geom_point(size=.5)+theme_bw()+
  geom_vline(xintercept = 2000)+ annotate("text",x=2200,y=50,label=c("2000"),hjust=0, size=2.8)
plot6 <- ggplot(FBL@meta.data, aes(FBL$percent.mt, FBL$percent.ribo))+geom_point(size=.5)+theme_bw()
plot1 + plot2 + plot3+plot4+plot5+plot6
dev.off()

## Filtering based on QC metrics
# Check filtering to be set
selected <- WhichCells(FBL, expression = nCount_RNA < 30000 & nFeature_RNA > 750 & 
                         nCount_RNA > 2000 &percent.mt < 4)
length(selected) #How many cells are left

# Filter cells
FBL_filtered <- subset(FBL, subset = nCount_RNA < 30000 & nFeature_RNA > 750 & 
                         nCount_RNA > 2000 &percent.mt < 4)

# Visualize QC metrics post-filtering
pdf(file="FBL_afterQC_Vlns.pdf", width = 20, height = 8)
VlnPlot(FBL_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
dev.off()


### 5. Normalization and cell cycle scoring---------------------------------------------------------------------------------------------------------
DefaultAssay(FBL_filtered) <- 'RNA'
## Full regression for cell cycle
# Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
FBL_filtered <- NormalizeData(FBL_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
# Segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign cell cycle scores to the cells 
FBL_filtered <- CellCycleScoring(object = FBL_filtered, s.features = s.genes, g2m.features = g2m.genes, 
                                 set.ident = FALSE)
head(x = FBL_filtered@meta.data)

FBL_filtered <- ScaleData(FBL_filtered)
FBL_filtered <- FindVariableFeatures(FBL_filtered)

## remove TCR genes from variable features 
VariableFeatures(FBL_filtered) <- VariableFeatures(FBL_filtered)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FBL_filtered))]

FBL_filtered <- RunPCA(FBL_filtered)

### 6. ADT data ------------------------------------------------------------------------------
DefaultAssay(FBL_filtered) <- 'ADT'
FBL_filtered <- NormalizeData(FBL_filtered, normalization.method = 'CLR', margin = 2)

rownames(FBL_filtered@assays$ADT$data)
FBL_filtered$CD4_ADT <- FBL_filtered@assays$ADT$data[6,]
FBL_filtered$CD8_ADT <- FBL_filtered@assays$ADT$data[7,]
FBL_filtered$CD45RA_ADT <- FBL_filtered@assays$ADT$data[3,]

VlnPlot(FBL_filtered, 'Hashtag2')
FBL_filtered$Donor <- rep('GD88', length(FBL_filtered$orig.ident))
FBL_filtered@meta.data[FBL_filtered@assays$ADT$data[18,]>4, 'Donor'] <- 'GD89'
FBL_filtered@meta.data[FBL_filtered@assays$ADT$data[19,]>4, 'Donor'] <-'GD108'
FBL_filtered@meta.data[FBL_filtered@assays$ADT$data[20,]>4, 'Donor'] <- 'GD129'
FBL_filtered@meta.data[FBL_filtered@assays$ADT$data[21,]>4, 'Donor'] <- 'GD136'
FBL_filtered@meta.data[FBL_filtered@assays$ADT$data[22,]>4, 'Donor'] <-'GD144'

DefaultAssay(FBL_filtered) <- 'RNA'

#### Save and load -----------------------------------------------------------------------------------------------------------

setwd("~/PhD/Fetal 10X/UsedObjects")
save(FBL_filtered, file = 'FBL_filtered_simple.Robj')
save(FBL_CD4, file = 'FBL_CD4subset_wTreg.Robj')

## load
setwd("~/PhD/Fetal 10X/UsedObjects")
load('FBL_filtered_simple.Robj')
load('FBL_CD4subset_wTreg.Robj')

### 7. Clustering ----------------------------------------------------------------------------------------------------------------------------
DefaultAssay(FBL_filtered) <- 'RNA'
FBL_filtered <- RunUMAP(FBL_filtered, dims = 1:30)

FBL_filtered <- FindNeighbors(object = FBL_filtered, reduction = "pca", dims = 1:30)
FBL_filtered <- FindClusters(object = FBL_filtered, resolution = seq(0.1:1, by=0.1))

DimPlot(object = FBL_filtered, reduction = "umap", pt.size=.2, group.by = 'RNA_snn_res.0.2')

clustree(FBL_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')

FeaturePlot(FBL_filtered, 'CD4.1',reduction='umap', pt.size=.2)
FeaturePlot(FBL_filtered, 'CD8a',reduction='umap', pt.size=.2)
FeaturePlot(FBL_filtered, 'CD45RA',reduction='umap', pt.size=.2)
FeaturePlot(FBL_filtered, 'Hashtag2',reduction='umap', pt.size=.2)
DimPlot(FBL_filtered, group.by = 'Donor',reduction='umap', pt.size=.2)
DimPlot(FBL_filtered, reduction='umap', group.by = 'Phase')
FeaturePlot(FBL_filtered, 'nFeature_RNA', reduction='umap')
FeaturePlot(FBL_filtered, 'percent.mt', reduction='umap')

### 8. Figures ---------------------------------------------------------------------------------------------------
setwd("")

FBL_filtered@meta.data$FBL_res.0.2_clusters <- as.factor(FBL_filtered@meta.data$RNA_snn_res.0.2)
levels(FBL_filtered@meta.data$FBL_res.0.2_clusters) <- c('0: CD4 Naive/TCM', '1: CD8 Naive-like', 
                                                         '2: Donor 088',
                                                         '3: Treg', '4: CD4/CD8 early memory/Th17', '5: Proliferating', '6: CD4/CD8 low quality',
                                                         '7: CD8 Effector/memory', 
                                                         '8: B cell contamination')
FBL_filtered@meta.data$FBL_res.0.2_clusters <- factor(FBL_filtered@meta.data$FBL_res.0.2_clusters,
                                                      levels=c( '8: B cell contamination','5: Proliferating',
                                                                '2: Donor 088', '6: CD4/CD8 low quality',
                                                                '1: CD8 Naive-like',  '7: CD8 Effector/memory',
                                                                '4: CD4/CD8 early memory/Th17', 
                                                                '0: CD4 Naive/TCM', '3: Treg'))


FigS3A<- ((FeaturePlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.5,
                       label.size = 12,
                       features = c('CD4.1'), cols = c('gold', 'purple'),
                       label = T)+ggtitle('CD4 surface expression')) + 
            (FeaturePlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.5,
                         label.size = 12,
                         features = c('CD8a'), cols = c('gold', 'purple'),
                         label = T)+ggtitle('CD8A surface expression'))+
            DimPlot(FBL_filtered, reduction='umap', pt.size=.5, group.by ='Phase', cols=c('grey20', 'violet', 'seagreen1')) + 
            DimPlot(FBL_filtered, reduction='umap', pt.size=.5, group.by ='Donor', cols=palette.colors(palette = "Accent"))  &
            theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                  legend.text = element_text(size=26),legend.title = element_text(size=28),
                  legend.key.size = unit(20,'points'),
                  plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=2),plot.subtitle = element_text(size=28, hjust=0.5),
                  plot.margin = margin(20,20,20,20),
                  text=element_text(size=8))) +
  (VlnPlot(object = SetIdent(FBL_filtered, value = "FBL_res.0.2_clusters"),  
           features = c('CD45RA'), pt.size=0,
           cols=c('grey','red4','grey30','#669999','mistyrose', 'plum3',
                  'deepskyblue3','blue', 'turquoise')) +
     ggtitle('CD45RA surface expression')+
     theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=22), 
           axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
           legend.text = element_text(size=14),
           plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),
           plot.margin = margin(10,10,10,10),
           text=element_text(size=14))+
     scale_x_discrete(labels=c('8','5','2','6','1','7','4', '0','3'))+
     xlab('Cluster ID')+
     guides(fill='none'))+
  (ggplot(FBL_filtered@meta.data, aes(x=orig.ident, fill=FBL_res.0.2_clusters)) + theme_classic() +
     geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
     geom_text(stat = 'count', size=6,
               position = position_fill(vjust=.5), aes(color=FBL_res.0.2_clusters,
                                                       label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
     labs(fill = "ClusterID", title='Proportion of T cell clusters')+
     theme(axis.text.x = element_blank(), axis.title.x = element_text(size=22), 
           axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
           legend.text = element_text(size=24),legend.title = element_text(size=26),
           plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),plot.subtitle = element_text(size=28, hjust=0.5),
           plot.margin = margin(40,40,40,40),
           text=element_text(size=8))+
     scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
     scale_fill_manual(values=c('grey','red4','grey30','#669999','mistyrose', 'plum3',
                                'deepskyblue3','blue', 'turquoise'))+
     scale_color_manual(values=c('grey10','grey80', 'grey80',rep('grey20',4), 'grey60', 'grey20'))+
     guides(color='none'))+ 
  plot_layout(ncol=2)+
  plot_annotation(title = "Fetal blood", 
                  theme=theme(plot.title=element_text(size=50, hjust=0.35, vjust=4), 
                              plot.margin=margin(60,30,10,30)))
FigS3A 
ggsave(file='FBL_overview_res.0.2_SupplFig3A.pdf', FigS3A, width=20, height=18)



pdf('FBL_SurfaceData_res.0.2.pdf', width = 18, height=32)
FeaturePlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.2,
            features =FBL_filtered@assays[["ADT"]]@counts@Dimnames[[1]], 
            label = T, label.size = 14) + plot_layout(ncol=3) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=10))
dev.off()

pdf('FBL_SurfaceData_Vln_res.0.2.pdf', width = 18, height=36)
VlnPlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"),  
        features = FBL_filtered@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
        cols=c('black','mistyrose','grey','turquoise',"gold",'red4','#669999','plum3','grey40'),
        group.by = 'RNA_snn_res.0.2') + plot_layout(ncol=3)& 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_Clustering_res.0.2.pdf', height=12, width=12)
DimPlot(object = FBL_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.2', 
        label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
        cols=c( 'black','mistyrose','grey','turquoise',"gold",'red4','#669999','plum3','grey40')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Blood', color='Cluster ID')+
  scale_color_manual(values=c('black','mistyrose','grey','turquoise',"gold",'red4','#669999','plum3','grey40'))
dev.off()

pdf('FBL_Clustering_res.0.2_clusternames.pdf', height=12, width=18)
DimPlot(object = FBL_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.2', 
        label = T, label.size = 22,label.color = 'grey50', repel=T,label.box = T,
        cols=c('black','mistyrose','grey','turquoise',"gold",'red4','#669999','plum3','grey40')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Blood', color='Cluster ID')+
  scale_color_manual(values=c('black','mistyrose','grey','turquoise',"gold",'red4','#669999','plum3','grey40'),
                     labels= c('0: CD4 Naive/TCM', '1: CD8 Naive-like', 
                               '2: Donor 088',
                               '3: Treg', '4: CD4 early development/Th17', '5: Proliferating', '6: CD4/CD8 low quality',
                               '7: CD8 Effector/memory', 
                               '8: B cell contamination'))
dev.off()

pdf('FBL_Donor.pdf', height=12, width=12)
DimPlot(object = FBL_filtered, reduction = "umap", pt.size=1, group.by = 'Donor') +
  ggtitle('Donor ID') +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  scale_color_manual(values=palette.colors(palette = "R4"))
dev.off()

pdf('FBL_QCData_res.0.2.pdf', width = 18, height=10)
FeaturePlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.4,
            features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
            label = T, label.size = 10) + DimPlot(FBL_filtered, reduction='umap', pt.size=.4, group.by ='Phase')+
  plot_layout(ncol=3) & 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))
dev.off()

pdf('FBL_QCData_Vln_res.0.2.pdf', width = 14, height=5)
VlnPlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"),  
        features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
        cols=c('black','mistyrose','grey','turquoise',"gold",'red4','#669999','plum3','grey40')) +
  plot_layout(nrow=1, ncol=4) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_overview_res.0.2.pdf', width = 18, height=16)
FeaturePlot(object = SetIdent(FBL_filtered, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.5,
            label.size = 12,
            features = c('CD4.1', 'CD8a', 'CD45RA'),
            label = T) + DimPlot(FBL_filtered, reduction='umap', pt.size=.5, group.by ='Phase') + plot_layout(ncol=2) &
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))
dev.off()


pdf('FBL_Clustree_colouredbyCD4ADT.pdf')
clustree(FBL_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')
dev.off()


FBL_filtered@meta.data$FBL_res.0.2_clusters <- as.factor(FBL_filtered@meta.data$RNA_snn_res.0.2)
levels(FBL_filtered@meta.data$FBL_res.0.2_clusters) <- c('0: CD4 Naive/TCM', '1: CD8 Naive-like', 
                                                         '2: Donor 088',
                                                         '3: Treg', '4: CD4 early development/Th17', '5: Proliferating', '6: CD4/CD8 low quality',
                                                         '7: CD8 Effector/memory', 
                                                         '8: B cell contamination')
FBL_filtered@meta.data$FBL_res.0.2_clusters <- factor(FBL_filtered@meta.data$FBL_res.0.2_clusters,
                                                      levels=c( '8: B cell contamination','5: Proliferating',
                                                                '2: Donor 088', '6: CD4/CD8 low quality',
                                                                '1: CD8 Naive-like',  '7: CD8 Effector/memory',
                                                                '4: CD4 early development/Th17', 
                                                                '0: CD4 Naive/TCM', '3: Treg'))
pdf('FBL_Proportionplot_total_res.0.2.pdf')
ggplot(FBL_filtered@meta.data, aes(x=orig.ident, fill=FBL_res.0.2_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  geom_text(stat = 'count', col='grey40',size=4,
            position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
  labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal Blood')+
  theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('grey40','red4','grey','#669999','mistyrose','plum3',"gold", 'black','turquoise'))
dev.off()


pdf('FBL_Proportionplot_ClustersperDonor_res.0.2.pdf', width=8,height=5)
ggplot(FBL_filtered@meta.data, aes(x=Donor, fill=FBL_res.0.2_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='', subtitle='')+
  theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('grey40','red4','grey','#669999','mistyrose','plum3',"gold", 'black','turquoise'))+
  coord_flip()
dev.off()

pdf('FBL_Proportionplot_DonorsperCluster_res.0.2.pdf')
ggplot(FBL_filtered@meta.data, aes(fill=Donor, x=RNA_snn_res.0.2)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  labs(fill = "Donor ID", title='', subtitle='')+
  theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=palette.colors(palette = "R4"))
dev.off()

#### 10. subset naive and memory CD4 T cells - with Treg ----------------------------------------------------------------------------------
FBL_filtered <- SetIdent(FBL_filtered, value = "RNA_snn_res.0.2")
FBL_CD4 <- subset(FBL_filtered,  idents= c('0','3','4'), subset=CD8_ADT<1&CD4_ADT>1)
## recluster
FBL_CD4 <- FindVariableFeatures(FBL_CD4)
VariableFeatures(FBL_CD4) <- VariableFeatures(FBL_CD4)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FBL_CD4))]
FBL_CD4 <- ScaleData(FBL_CD4)
FBL_CD4 <- RunPCA(FBL_CD4)
FBL_CD4 <- RunUMAP(FBL_CD4, dims = 1:30)
FBL_CD4 <- FindNeighbors(object = FBL_CD4, reduction = "pca", dims = 1:30)
FBL_CD4 <- FindClusters(object = FBL_CD4, resolution = seq(0.1:1, by=0.1))

DimPlot(object = FBL_CD4, reduction = "umap", pt.size=.5, group.by = 'RNA_snn_res.0.3')
DimPlot(object = FBL_CD4, reduction = "umap", pt.size=.5, group.by = 'FBL_res.0.2_clusters')
DimPlot(object = FBL_CD4, reduction = "umap", pt.size=.5, group.by = 'donor_id')

clustree(FBL_CD4wTreg)

FeaturePlot(FBL_CD4wTreg, 'CD4.1',reduction='umap', pt.size=.2)
FeaturePlot(FBL_CD4wTreg, 'CD8a',reduction='umap', pt.size=.2)
FeaturePlot(FBL_CD4wTreg, 'CD45RA',reduction='umap', pt.size=1)
FeaturePlot(FBL_CD4wTreg, 'CD27.1',reduction='umap', pt.size=.2)
DimPlot(FBL_CD4wTreg, reduction='umap', group.by = 'Phase')

setwd("")
pdf('FBL_CD4wTreg_SurfaceData_res0.6.pdf', width = 18, height=18)
FeaturePlot(object = SetIdent(FBL_CD4wTreg, value = "RNA_snn_res.0.6"), reduction = "umap", pt.size=.5,
            features = c('CD4.1', 'CD8a', 'CD45RA', 'CD31', 'CD27.1','CD28.1',  'CXCR3', 'CD69', 'CD103', 'CD62L'), 
            label = T, label.size = 10) + plot_layout(ncol=3) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=10))
dev.off()

pdf('FBL_CD4wTreg_SurfaceData_Vln_res0.6.pdf', width = 18, height=18)
VlnPlot(object = SetIdent(FBL_CD4wTreg, value = "RNA_snn_res.0.6"),  
        features = c('CD4.1', 'CD8a', 'CD45RA', 'CD31', 'CD27.1', 'CD28.1', 'CXCR3', 'CD69', 'CD103', 'CD62L'), pt.size=0,
        group.by = 'RNA_snn_res.0.6') + plot_layout(ncol=3)& 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()


pdf('FBL_CD4wTreg_Clustering_res0.6.pdf', height=12, width=14)
DimPlot(object = FBL_CD4wTreg, reduction = "umap", pt.size=.5, group.by = 'RNA_snn_res.0.6', label = T, label.size = 10) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_CD4wTreg_OldClusters_totalTcells.pdf', height=12, width=14)
DimPlot(object = FBL_CD4wTreg, reduction = "umap", pt.size=.5, group.by = 'FBL_res.0.4_clusters', label = F, label.size = 10) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_CD4wTreg_Donor_res0.6.pdf', height=12, width=14)
DimPlot(object = FBL_CD4wTreg, reduction = "umap", pt.size=.5, group.by = 'donor_id') +
  ggtitle('Donor ID') & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=10),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_CD4wTreg_doubles_res0.6.pdf', height=12, width=14)
DimPlot(object = FBL_CD4wTreg, reduction = "umap", pt.size=.4, group.by = 'DF.classifications_0.25_0.22_446') & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=10),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_CD4wTreg_OverviewSurfaceQC_res0.6.pdf', width = 18, height=12)
FeaturePlot(object = SetIdent(FBL_CD4wTreg, value = "RNA_snn_res.0.6"), reduction = "umap", pt.size=.5,
            features = c('CD4.1', 'CD8a', 'CD45RA', 'nFeature_RNA', 'percent.mt'),
            label = T, label.size = 6) + DimPlot(FBL_CD4wTreg, reduction='umap', pt.size=.5, group.by ='Phase')+
  plot_layout(ncol=3) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_CD4wTreg_QCData_Vln_res0.6.pdf', width = 14, height=5)
VlnPlot(object = SetIdent(FBL_CD4wTreg, value = "RNA_snn_res.0.6"),  
        features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0) +
  plot_layout(nrow=1, ncol=4) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('FBL_CD4wTreg_Clustree.pdf')
clustree(FBL_CD4wTreg)
dev.off()

#### 11. DEG CD4 --------------------------------------------------
#### ROC/MAST for cluster defining markers 
DefaultAssay(object = FBL_CD4wTreg) <- "RNA"
FBL_CD4wTreg <- SetIdent(FBL_CD4wTreg, value = "RNA_snn_res.0.6")

##Find Markers that are specific for each cluster
Batchedmarkers.mast_i_RNA_data_0.6=FindAllMarkers(FBL_CD4wTreg, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_0.6$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.6$pct.1/Batchedmarkers.mast_i_RNA_data_0.6$pct.2

## Create list
listDEgenes_i_RNA_MAST_0.6<- split(Batchedmarkers.mast_i_RNA_data_0.6, f=Batchedmarkers.mast_i_RNA_data_0.6$cluster)

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
listDEgenes_i_RNA_MAST_0.6 <-lapply(listDEgenes_i_RNA_MAST_0.6, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
listDEgenes_i_RNA_MAST_0.6<-lapply(listDEgenes_i_RNA_MAST_0.6,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

