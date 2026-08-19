#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
library(Seurat)
library(dplyr)
library(sctransform)
library(ggplot2)
library(Matrix)
library(patchwork)
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
#library(DoubletFinder)
library(scRepertoire)
library(stringr)
library(ggthemes)
library(VDJdive)
library(S4Vectors)
library(purrr)

### 1. Setting working directory and loading data----------------------------------------------------------------
setwd("~/")

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
FSI_data <- Read10X(data.dir = getwd()) 

### 2. Seurat object  -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data
FSI <- CreateSeuratObject(counts = FSI_data$`Gene Expression`, project = "FSI")
FSI

##add antibody data to new assay (called ADT) https://satijalab.org/seurat/articles/multimodal_vignette
FSI_ADT <- CreateAssayObject(counts=FSI_data$`Antibody Capture`)
FSI[['ADT']] <- FSI_ADT

rownames(FSI[['ADT']])
rownames(FSI[['RNA']])

rm(FSI_data)
rm(FSI_ADT)

### 3. Add info to metadata ----------------------------------------------------------------------------------
seurat_G <- readRDS("~/all_ab_old_harmony_standarized.rds")

colnames(seurat_G@meta.data)
rownames(seurat_G@meta.data)

## extract only the FSI data and remove FSI from rownames
FSI_seurat_G <- subset(seurat_G, subset=orig.ident=='FSI_ab') 
rownames(FSI_seurat_G@meta.data) <- str_remove(rownames(FSI_seurat_G@meta.data), 'FSI_')
rownames(FSI_seurat_G@meta.data)

# donor info = 'donor_id_test' column
FSI@meta.data <- cbind(FSI@meta.data, 
                                FSI_seurat_G@meta.data[match(rownames(FSI@meta.data),
                                                             rownames(FSI_seurat_G@meta.data)),
                                                       c(12)])
colnames(FSI@meta.data)[6] <- 'Donor'


### 4. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
setwd("~/")
DefaultAssay(FSI) <- 'RNA'
  
## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  FSI <- PercentageFeatureSet(FSI, pattern = "^MT-", col.name = 'percent.mt')
  FSI <- PercentageFeatureSet(FSI, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(FSI@assays$RNA@counts)[grep('MT', rownames(FSI@assays$RNA@counts))]) #somewhere in name
  length(rownames(FSI@assays$RNA@counts)[grep('^MT\\.', rownames(FSI@assays$RNA@counts))]) #start of name followed by .
  length(rownames(FSI@assays$RNA@counts)[grep('^MT-', rownames(FSI@assays$RNA@counts))]) #start of name followed by -
  length(rownames(FSI@assays$RNA@counts)[grep('^RP', rownames(FSI@assays$RNA@counts))]) #start of name followed by -
  
  
  # Show QC metrics for the first 5 cells
  head(FSI@meta.data, 5)
  
  
  # Visualize QC metrics as a violin plot
  pdf(file="FSI_beforeQC_Vlns.pdf", width = 20, height = 8)
  VlnPlot(FSI, features = c("nFeature_RNA", "nCount_RNA", "nCount_ADT","percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="FSI_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(FSI, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)
  plot1.z <- FeatureScatter(FSI, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)+
    scale_x_continuous(limits=c(0,5000))
  plot2 <- FeatureScatter(FSI, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.00001)
  plot3 <- FeatureScatter(FSI, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.00001)
  plot1 + plot1.z+  plot2 + plot3 + plot_layout(ncol=4)
  dev.off()
  
  #With filtering cut-offs
  pdf(file="FSI_correlations QC_ggplot_withfilterlingablines_nFeature.pdf", width = 24, height = 12)
  plot1 <- ggplot(FSI@meta.data, aes(FSI$nCount_RNA, FSI$nFeature_RNA))+geom_point(size=.5)+theme_bw()+
    geom_vline(xintercept = 22000)+ annotate("text",x=24000,y=1,label=c("22000"),hjust=0, size=2.8)+
    geom_vline(xintercept = 1200)+ annotate("text",x=1500,y=1,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 750)+ annotate("text",x=100,y=900,label=c("750"),hjust=0, size=2.8)
  plot2 <- ggplot(FSI@meta.data, aes(FSI$nCount_RNA, FSI$percent.mt))+geom_point(size=.5)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=15,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=3000,y=8,label=c("7"),vjust=0, size=2.8)+
    xlim(0,5000)
  plot3 <- ggplot(FSI@meta.data, aes(FSI$nCount_RNA, FSI$percent.mt))+geom_point(size=.5)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=20,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=10000,y=9,label=c("7"),vjust=0, size=2.8)
  plot4 <- ggplot(FSI@meta.data, aes(FSI$nFeature_RNA, FSI$percent.mt))+geom_point(size=.5)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 750)+ annotate("text",x=900,y=20,label=c("750"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=10000,y=9,label=c("7"),vjust=0, size=2.8)+
    xlim(0,1000)
  plot5 <- ggplot(FSI@meta.data, aes(FSI$nCount_RNA, FSI$percent.ribo))+geom_point(size=.5)+theme_bw()+
    geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=50,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=50,y=10,label=c("7"),hjust=0, size=2.8)+xlim(0,5000)
  plot6 <- ggplot(FSI@meta.data, aes(FSI$percent.mt, FSI$percent.ribo))+geom_point(size=.5)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5+plot6
  dev.off()
  
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(FSI, expression = nCount_RNA < 22000 & nFeature_RNA > 750 
                         & nCount_RNA > 1200 & percent.mt < 7 & percent.ribo > 7)
  length(selected) #How many cells are left
  
  # Filter cells
  FSI_filtered <- subset(FSI, subset = nCount_RNA < 22000 & nFeature_RNA > 750 
                         & nCount_RNA > 1200 &  percent.mt < 7 & percent.ribo > 7)
  
  # Visualize QC metrics post-filtering
  pdf(file="FSI_afterQC_Vlns.pdf", width = 20, height = 8)
  VlnPlot(FSI_filtered, features = c("nFeature_RNA", "nCount_RNA", "nCount_ADT","percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
  dev.off()
  
  ## remove doublets
  FSI_filtered$doublets <- ifelse(FSI_filtered$Donor=='doublet', T, F)
  FSI_filtered@meta.data[is.na(FSI_filtered$doublets),'doublets'] <- F
  FSI_filtered <- subset(FSI_filtered, subset= doublets==F)

### 5. Normalization and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  FSI_filtered <- NormalizeData(FSI_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  FSI_filtered <- CellCycleScoring(object = FSI_filtered, s.features = s.genes, g2m.features = g2m.genes, 
                                   set.ident = FALSE)
  head(x = FSI_filtered@meta.data)
  
  FSI_filtered <- ScaleData(FSI_filtered)
  FSI_filtered <- FindVariableFeatures(FSI_filtered)
  
  ## remove TCR genes from variable features 
  VariableFeatures(FSI_filtered) <- VariableFeatures(FSI_filtered)[!grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FSI_filtered))]
  
  FSI_filtered <- RunPCA(FSI_filtered)
  
  
### 6. ADT data  ------------------------------------------------------------------------------
  DefaultAssay(FSI_filtered) <- 'ADT'
  FSI_filtered <- NormalizeData(FSI_filtered, normalization.method = 'CLR', margin = 2)
  
  rownames(FSI_filtered@assays$ADT$data)
  FSI_filtered$CD4_ADT <- FSI_filtered@assays$ADT$data[6,]
  FSI_filtered$CD8_ADT <- FSI_filtered@assays$ADT$data[7,]
  
  DefaultAssay(FSI_filtered) <- 'RNA'
  
 
#### Save and load -----------------------------------------------------------------------------------------------------------
  
  setwd("~/")
  save(FSI_filtered, file = 'FSI_filtered_simple.Robj')
  save(FSI_CD4, file = 'FSI_CD4subset_wTreg_TCR.Robj')

  ## load
  setwd("~/")
  load('FSI_filtered_simple.Robj')
  load('FSI_CD4subset_wTreg.Robj')
  
  
### 7. Clustering ----------------------------------------------------------------------------------------------------------------------------
  DefaultAssay(FSI_filtered) <- 'RNA'
  FSI_filtered <- RunUMAP(FSI_filtered, dims = 1:30)
  
  FSI_filtered <- FindNeighbors(object = FSI_filtered, reduction = "pca", dims = 1:30)
  FSI_filtered <- FindClusters(object = FSI_filtered, resolution = seq(0.1:1, by=0.1))
  
  DimPlot(object = FSI_filtered, reduction = "umap", pt.size=.2, group.by = 'RNA_snn_res.0.5')
  
  
  clustree(FSI_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')
  
  FeaturePlot(FSI_filtered, 'CD4.1',reduction='umap', pt.size=.2)
  FeaturePlot(FSI_filtered, 'CD8a',reduction='umap', pt.size=.2)
  FeaturePlot(FSI_filtered, 'CD45RA',reduction='umap', pt.size=.2)
  FeaturePlot(FSI_filtered, 'CD27.1',reduction='umap', pt.size=.2)
  DimPlot(FSI_filtered, reduction='umap', group.by = 'Phase')
  
### 8. Figures ---------------------------------------------------------------------------------------------------
  setwd("")
  donor_colors <- c(palette.colors(palette = "R4")[1:4],palette.colors(palette = "R4")[6:7])
  
  FSI_filtered@meta.data$FSI_res.0.5_clusters <- as.factor(FSI_filtered@meta.data$RNA_snn_res.0.5)
  levels(FSI_filtered@meta.data$FSI_res.0.5_clusters) <- c('0: CD4 Th1/TRM', '1: CD4 Naive-like/T-CM', 
                                                           '2: CD8 Naive-like','3: Treg', 
                                                           '4: CD4 Th2/early development', '5: CD4 Activated', '6: CD4 Th17', '7: CD4 - low quality', 
                                                           '8: CD8 effector/TRM', '9: Proliferating - S/G2M phase')
  FSI_filtered@meta.data$FSI_res.0.5_clusters <- factor(FSI_filtered@meta.data$FSI_res.0.5_clusters,
                                                        levels=c('9: Proliferating - S/G2M phase', 
                                                                 '2: CD8 Naive-like', '8: CD8 effector/TRM','7: CD4 - low quality',
                                                                 '0: CD4 Th1/TRM', '6: CD4 Th17','4: CD4 Th2/early development','5: CD4 Activated',
                                                                 '1: CD4 Naive-like/T-CM',
                                                                 '3: Treg'))
  
  
  FigS1B <- ((FeaturePlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"), reduction = "umap", pt.size=.5,
                         label.size = 12,
                         features = c('CD4.1'), cols = c('gold', 'purple'),
                         label = T)+ggtitle('CD4 surface expression')) + 
              (FeaturePlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"), reduction = "umap", pt.size=.5,
                           label.size = 12,
                           features = c('CD8a'), cols = c('gold', 'purple'),
                           label = T)+ggtitle('CD8A surface expression'))+
              DimPlot(FSI_filtered, reduction='umap', pt.size=.5, group.by ='Phase', cols=c('grey20', 'violet', 'seagreen1')) + 
              DimPlot(FSI_filtered, reduction='umap', pt.size=.5, group.by ='Donor', cols=donor_colors)  &
              theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                    legend.text = element_text(size=26),legend.title = element_text(size=28),
                    legend.key.size = unit(20,'points'),
                    plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=2),plot.subtitle = element_text(size=28, hjust=0.5),
                    plot.margin = margin(20,20,20,20),
                    text=element_text(size=8))) +
    (VlnPlot(object = SetIdent(FSI_filtered, value = "FSI_res.0.5_clusters"),  
             features = c('CD45RA'), pt.size=0,
             cols=c('red4','mistyrose', 'plum3',
                    '#669999',"deepskyblue","cadetblue2","cadetblue1", 'deepskyblue3','blue', 'turquoise')) +
              scale_x_discrete(labels=c('9','2','8','7','0','6','4','5','1','3'))+
       ggtitle('CD45RA surface expression')+
       theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=22), 
             axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
             legend.text = element_text(size=14),
             plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),
             plot.margin = margin(10,10,10,10),
             text=element_text(size=14))+
       xlab('Cluster ID')+
       guides(fill='none'))+
    (ggplot(FSI_filtered@meta.data, aes(x=orig.ident, fill=FSI_res.0.5_clusters)) + theme_classic() +
       geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
       geom_text(stat = 'count', size=6,
                 position = position_fill(vjust=.5), aes(color=FSI_res.0.5_clusters,
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
                                  '#669999',"deepskyblue","cadetblue2","cadetblue1", 'deepskyblue3','blue', 'turquoise'))+
       scale_color_manual(values=c('grey80', rep('grey20',7), 'grey60', 'grey10'))+
       guides(color='none'))+ 
    plot_layout(ncol=2)+
    plot_annotation(title = "Fetal small intestine", 
                    theme=theme(plot.title=element_text(size=50, hjust=0.4, vjust=4), 
                                plot.margin=margin(60,30,10,30)))
  FigS1B
  ggsave(file='FSI_overview_res.0.5_SupplFig1B.pdf', FigS1B, width=20, height=18)

  pdf('FSI_SurfaceData_res.0.5.pdf', width = 18, height=18)
  FeaturePlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"), reduction = "umap", pt.size=.2,
              features =FSI_filtered@assays[["ADT"]]@counts@Dimnames[[1]], 
              label = T, label.size = 14) + plot_layout(ncol=3) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
          text=element_text(size=10))
  dev.off()
  
  pdf('FSI_SurfaceData_Vln_res.0.5.pdf', width = 18, height=18)
  VlnPlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"),  
          features = FSI_filtered@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
          cols=c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4'),
          group.by = 'RNA_snn_res.0.5') + plot_layout(ncol=3)& 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  
  pdf('FSI_Clustering_res.0.5.pdf', height=12, width=12)
  DimPlot(object = FSI_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.5', 
          label = T, label.size = 22, label.color = 'grey60', repel=T,label.box = T,
          cols=c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4')) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Intestine', color='Cluster ID')
 dev.off()
  
  pdf('FSI_Clustering_res.0.5_clusternames.pdf', height=12, width=18)
  DimPlot(object = FSI_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.5', 
          label = T, label.size = 22, label.color = 'grey60', repel=T,label.box = T,
          cols=c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4')) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Intestine', color='Cluster ID')+
    scale_color_manual(values=c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4'),
                       labels= c('0: CD4 Th1/TRM', '1: CD4 Naive-like/T-CM', 
                                 '2: CD8 Naive-like','3: Treg', 
                                 '4: CD4 Th2/early development', '5: CD4 Activated', '6: CD4 Th17', '7: CD4 - low quality', 
                                 '8: CD8 effector/TRM', '9: Proliferating - S/G2M phase'))
  dev.off()
  
  pdf('FSI_Donor.pdf', height=12, width=12)
  DimPlot(object = FSI_filtered, reduction = "umap", pt.size=1, group.by = 'Donor') +
    ggtitle('Donor ID') +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    scale_color_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  pdf('FSI_QCData_res.0.5.pdf', width = 18, height=10)
  FeaturePlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"), reduction = "umap", pt.size=.4,
              features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
              label = T, label.size = 10) + DimPlot(FSI_filtered, reduction='umap', pt.size=.4, group.by ='Phase')+
    plot_layout(ncol=3) & 
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()
  
  pdf('FSI_QCData_Vln_res.0.5.pdf', width = 14, height=5)
  VlnPlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"),  
          features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
          cols=c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4')) +
    plot_layout(nrow=1, ncol=4) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  pdf('FSI_overview_res.0.5.pdf', width = 18, height=16)
  FeaturePlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"), reduction = "umap", pt.size=.5,
              label.size = 12,
              features = c('CD4.1', 'CD8a', 'CD45RA'),
              label = T) + DimPlot(FSI_filtered, reduction='umap', pt.size=.5, group.by ='Phase') + plot_layout(ncol=2) &
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()
  
  pdf('FSI_CD45RA_res.0.5.pdf', width = 6, height=6)
  FeaturePlot(object = SetIdent(FSI_filtered, value = "RNA_snn_res.0.5"), reduction = "umap", pt.size=.5,
              label.size = 8,
              features = c('CD45RA'),
              label = T,label.color = c('blue','blue','black','blue','blue','blue','turquoise3','blue','red4','black'))+
    scale_color_gradient(limits=c(0,4),low="lightgrey", high= "blue")+
    ggtitle('CD45RA surface expression')+
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()

  
  pdf('FSI_Clustree_colouredbyCD4ADT.pdf')
  clustree(FSI_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')
  dev.off()
  
 
  FSI_filtered@meta.data$FSI_res.0.5_clusters <- as.factor(FSI_filtered@meta.data$RNA_snn_res.0.5)
  levels(FSI_filtered@meta.data$FSI_res.0.5_clusters) <- c('0: CD4 Th1/TRM', '1: CD4 Naive-like/T-CM', 
                                                           '2: CD8 Naive-like','3: Treg', 
                                                           '4: CD4 Th2/early development', '5: CD4 Activated', '6: CD4 Th17', '7: CD4 - low quality', 
                                                           '8: CD8 effector/TRM', '9: Proliferating - S/G2M phase')
  FSI_filtered@meta.data$FSI_res.0.5_clusters <- factor(FSI_filtered@meta.data$FSI_res.0.5_clusters,
                                                        levels=c('9: Proliferating - S/G2M phase', 
                                                                 '2: CD8 Naive-like', '8: CD8 effector/TRM','7: CD4 - low quality',
                                                                 '0: CD4 Th1/TRM', '6: CD4 Th17','4: CD4 Th2/early development','5: CD4 Activated',
                                                                 '1: CD4 Naive-like/T-CM',
                                                                 '3: Treg'))
  pdf('FSI_Proportionplot_total_res.0.5.pdf')
  ggplot(FSI_filtered@meta.data, aes(x=orig.ident, fill=FSI_res.0.5_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    geom_text(stat = 'count', col='grey40',size=4,
              position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
    labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal intestine')+
    theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('red4','mistyrose', 'plum3','#669999',
                               "#E69F00","#F0E442","#009E73", 'purple','black', 'turquoise'))
  dev.off()
  
  pdf('FSI_Proportionplot_total_res.0.5_simplecolors.pdf', height=9, width=10)
  ggplot(FSI_filtered@meta.data, aes(x=orig.ident, fill=FSI_res.0.5_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    geom_text(stat = 'count', size=5,
              position = position_fill(vjust=.5), aes(color=FSI_res.0.5_clusters,
                                                      label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
    labs(fill = "ClusterID", title='Fetal small intestine', subtitle='Proportion of T cell clusters')+
    theme(axis.text.y = element_text(size=16), axis.title.y = element_text(size=18,vjust=3), 
          axis.text.x = element_text(size=20),
          legend.text = element_text(size=16), legend.title = element_text(size=18),
          plot.title = element_text(size=28, hjust=0.5, , vjust=4), plot.subtitle = element_text(size=24, hjust=0.5, , vjust=4),
          plot.margin = margin(30,30,30,30),
          text=element_text(size=16))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('red4','mistyrose', 'plum3',
                               '#669999',"deepskyblue","cadetblue2","cadetblue1", 'deepskyblue3','blue', 'turquoise'))+
    scale_color_manual(values=c('grey80', rep('grey20',7), 'grey60', 'grey10'))+
    guides(color='none')
  dev.off()
   
  pdf('FSI_Proportionplot_ClustersperDonor_res.0.5.pdf', width=8,height=5)
  ggplot(FSI_filtered@meta.data, aes(x=Donor, fill=FSI_res.0.5_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    labs(fill = "Cluster", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('red4','#669999','mistyrose', 'plum3',
                               "#E69F00","#F0E442","#009E73", 'purple','black', 'turquoise'))+
    coord_flip()
  dev.off()
  
  pdf('FSI_Proportionplot_DonorsperCluster_res.0.5.pdf')
  ggplot(FSI_filtered@meta.data, aes(fill=Donor, x=RNA_snn_res.0.5)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "Donor ID", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  pdf('FSI_CD45RA_Vln_res.0.5.pdf', width = 12, height=6)
  VlnPlot(object = SetIdent(FSI_filtered, value = "FSI_res.0.5_clusters"),  
          features = c('CD45RA'), pt.size=0,
          cols=c('red4','mistyrose', 'plum3',
                 '#669999',"deepskyblue","cadetblue2","cadetblue1", 'deepskyblue3','blue', 'turquoise')) +
    ggtitle('CD45RA surface expression')+
    theme(axis.text = element_text(size=18), axis.title = element_text(size=20), legend.text = element_text(size=18),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10))+
    scale_x_discrete(labels=c('9','2','8','7','0','6','4','5','1','3'))
  dev.off()

  load("C:/Users/elise/Documents/PhD/Fetal 10X/DEG/FSI/listDEG_RNA_FSI_RNA_MAST_0.5.Robj")
  FSI_filtered <- SetIdent(FSI_filtered, value = "RNA_snn_res.0.5")
  top15 <- lapply(listDEgenes_i_RNA_MAST_0.5, function(df) {df %>% subset(avg_log2FC>0) %>% top_n(15, avg_log2FC)})
  top15 <- bind_rows(top15)
  
  pdf('FSI_heatmap_res.0.5.pdf', height=40, width=40)
  DoHeatmap(
    FSI_filtered,
    features = top15$gene,
    cells = NULL,
    group.by = "RNA_snn_res.0.5",
    group.bar = TRUE,
    group.colors = c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4'),
    disp.min = -2.5,
    disp.max = NULL,
    size = 6,
    hjust = 0,
    angle = 0,
    raster = TRUE,
    draw.lines = TRUE,
    lines.width = NULL,
    group.bar.height = 0.01,
    combine = TRUE,
    label=T) & 
    labs(colour='Cluster identity') &
    theme(text=element_text(size=28), legend.text = element_text(size=26), 
          legend.title = element_text(size=28), 
          plot.margin = margin(20,20,20,20))&
    scale_color_manual(labels=c('0: CD4 Th1/TRM', '1: CD4 Naive-like/T-CM', 
                             '2: CD8 Naive-like','3: Treg', 
                             '4: CD4 Th2/early development', '5: CD4 Activated', '6: CD4 Th17', '7: CD4 - low quality', 
                             '8: CD8 effector/TRM', '9: Proliferating - S/G2M phase'),
                       values=c("#E69F00",'black', 'mistyrose', 'turquoise',"#009E73", 'purple',"#F0E442",'#669999','plum3','red4'))
  dev.off()
  
  FSI_filtered_TCR$cloneSize <- as.factor(FSI_filtered_TCR$cloneSize)
  FSI_filtered_TCR$cloneSize_beta <- as.factor(FSI_filtered_TCR$cloneSize_beta)
  
  pdf('FSI_proportionplot_TCR_screp_both_res.0.5.pdf')
  ggplot(FSI_filtered_TCR@meta.data, aes(x=RNA_snn_res.0.5, fill=cloneSize)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "TCRab clone frequency", title='TCRab clonality - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(4,6,8,10)]),
                      labels=c('Small (3-5)', 'Rare (2)', 'Single (1)', 'NA'))
  dev.off()
  
  pdf('FSI_proportionplot_TCR_screp_both_NAremoved_res.0.5.pdf')
  ggplot(FSI_filtered_TCR@meta.data[!is.na(FSI_filtered_TCR@meta.data$cloneSize),], aes(x=RNA_snn_res.0.5, fill=cloneSize)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "TCRab clone frequency", title='TCRab clonality - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(4,6,8,10)]),
                      labels=c('Small (3-5)', 'Rare (2)', 'Single (1)', 'NA'))
  dev.off()
  
  pdf('FSI_countplot_TCR_screp_both_res.0.5.pdf')
  ggplot(FSI_filtered_TCR@meta.data, aes(x=RNA_snn_res.0.5, fill=cloneSize)) + theme_classic() +
    geom_bar(position = "stack") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "TCRab clone frequency", title='TCRab clonality - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(4,6,8,10)]),
                      labels=c('Small (3-5)', 'Rare (2)', 'Single (1)', 'NA'))
  dev.off()
  
  pdf('FSI_proportionplot_TCR_screp_beta_res.0.5.pdf')
  ggplot(FSI_filtered_TCR@meta.data, aes(x=RNA_snn_res.0.5, fill=cloneSize_beta)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "TCRb clone frequency", title='TCRb clonality - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(2,4,6,8,10)]),
                      labels=c('Medium (6-10)','Small (3-5)', 'Rare (2)', 'Single (1)', 'NA'))
  dev.off()
  
  pdf('FSI_proportionplot_TCR_screp_beta_NAremoved_res.0.5.pdf')
  ggplot(FSI_filtered_TCR@meta.data[!is.na(FSI_filtered_TCR@meta.data$cloneSize_beta),], 
         aes(x=RNA_snn_res.0.5, fill=cloneSize_beta)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "TCRb clone frequency", title='TCRb clonality - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(2,4,6,8,10)]),
                      labels=c('Medium (6-10)','Small (3-5)', 'Rare (2)', 'Single (1)', 'NA'))
  dev.off()
  
  pdf('FSI_countplot_TCR_screp_beta_res.0.5.pdf')
  ggplot(FSI_filtered_TCR@meta.data, aes(x=RNA_snn_res.0.5, fill=cloneSize_beta)) + theme_classic() +
    geom_bar(position = "stack") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "TCRb clone frequency", title='TCRb clonality - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(2,4,6,8)]),
                      labels=c('Medium (6-10)','Small (3-5)', 'Rare (2)', 'Single (1)', 'NA'))
  dev.off()
  
  # pdf('FSI_proportionplot_TCR_VDJdiveCustom_res.0.5.pdf')
  # ggplot(FSI_filtered_TCR@meta.data, aes(x=RNA_snn_res.0.5, fill=cloneSize_vdjd)) + theme_classic() +
  #   geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  #   labs(fill = "TCR clone frequency", title='TCR clonality - Fetal Intestine', subtitle='')+
  #   theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
  #         axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
  #   scale_y_continuous(expand = c(0,0))+
  #   scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(4,6,8,10)]))
  # dev.off()
  # 
  # pdf('FSI_proportionplot_TCR_VDJdiveCustom_NAremoved_res.0.5.pdf')
  # ggplot(FSI_filtered_TCR@meta.data[!is.na(FSI_filtered_TCR@meta.data$cloneSize_vdjd),], 
  #        aes(x=RNA_snn_res.0.5, fill=cloneSize_vdjd)) + theme_classic() +
  #   geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  #   labs(fill = "TCR clone frequency", title='TCR clonality - Fetal Intestine', subtitle='')+
  #   theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
  #         axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
  #   scale_y_continuous(expand = c(0,0))+
  #   scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(4,6,8,10)]))
  # dev.off()
  # 
  # pdf('FSI_countplot_TCR_VDJdiveCustom_res.0.5.pdf')
  # ggplot(FSI_filtered_TCR@meta.data, aes(x=RNA_snn_res.0.5, fill=cloneSize_vdjd)) + theme_classic() +
  #   geom_bar(position = "stack") + xlab("Cluster ID") + ylab("Fraction") + 
  #   labs(fill = "TCR clone frequency", title='TCR clonality - Fetal Intestine', subtitle='')+
  #   theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
  #         axis.title.y = element_text(vjust=2.5),plot.margin = margin(20,20,20,20))+
  #   scale_y_continuous(expand = c(0,0))+
  #   scale_fill_manual(values=c(palette.colors(palette = "Okabe-Ito")[c(4,6,8)]))
  # dev.off()
  # 
  pdf('FSI_clonaloverlap_jaccard_TCR_screp_both_res.0.5.pdf', width=8, height=7)
  clonalOverlap(FSI_filtered_TCR, cloneCall = 'aa', chain = 'both', method = 'jaccard', group.by = 'RNA_snn_res.0.5') +
    scale_x_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) +
    scale_y_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) + 
    xlab("Cluster ID") + ylab("Cluster ID") + 
    labs(fill = "Jaccard index", title='Shared TCRab-clones between clusters - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, vjust=-1,size=18),
          axis.title.y = element_text(vjust=3, size=16),axis.title.x = element_text(vjust=.5, size=14),
          axis.text = element_text(size=14),
          legend.title = element_text(vjust=2, size=14),legend.text = element_text(size=12),
          plot.margin = margin(20,20,20,20))
  dev.off()  

  pdf('FSI_clonaloverlap_jaccard_TCR_screp_beta_res.0.5.pdf', width=8, height=7)
  clonalOverlap(FSI_filtered_TCR, cloneCall = 'aa', chain = 'TRB', method = 'jaccard', group.by = 'RNA_snn_res.0.5') +
    scale_x_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) +
    scale_y_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) + 
    xlab("Cluster ID") + ylab("Cluster ID") + 
    labs(fill = "Jaccard index", title='Shared TCRb-clones between clusters - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, vjust=-1,size=18),
          axis.title.y = element_text(vjust=3, size=16),axis.title.x = element_text(vjust=.5, size=14),
          axis.text = element_text(size=14),
          legend.title = element_text(vjust=2, size=14),legend.text = element_text(size=12),
          plot.margin = margin(20,20,20,20))
  dev.off() 
  
  pdf('FSI_clonaloverlap_raw_TCR_screp_both_res.0.5.pdf', width=8, height=7)
  clonalOverlap(FSI_filtered_TCR, cloneCall = 'aa', chain = 'both', method = 'raw', group.by = 'RNA_snn_res.0.5') +
    scale_x_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) +
    scale_y_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) + 
    xlab("Cluster ID") + ylab("Cluster ID") + 
    labs(fill = "Raw count overlap", title='Shared TCRab-clones between clusters - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, vjust=-1,size=18),
          axis.title.y = element_text(vjust=3, size=16),axis.title.x = element_text(vjust=.5, size=14),
          axis.text = element_text(size=14),
          legend.title = element_text(vjust=2, size=14),legend.text = element_text(size=12),
          plot.margin = margin(20,20,20,20))
  dev.off() 
  
  pdf('FSI_clonaloverlap_raw_TCR_screp_beta_res.0.5.pdf', width=8, height=7)
  clonalOverlap(FSI_filtered_TCR, cloneCall = 'aa', chain = 'TRB', method = 'raw', group.by = 'RNA_snn_res.0.5') +
    scale_x_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) +
    scale_y_continuous(expand = c(0, 0), breaks = 0:9, labels = 0:9) + 
    xlab("Cluster ID") + ylab("Cluster ID") + 
    labs(fill = "Raw count overlap", title='Shared TCRb-clones between clusters - Fetal Intestine', subtitle='')+
    theme(plot.title = element_text(hjust=0.5, vjust=-1,size=18),
          axis.title.y = element_text(vjust=3, size=16),axis.title.x = element_text(vjust=.5, size=14),
          axis.text = element_text(size=14),
          legend.title = element_text(vjust=2, size=14),legend.text = element_text(size=12),
          plot.margin = margin(20,20,20,20))
  dev.off() 
  
### 9. DEG -----------------------------------------------------------------------------------------------------------
  #### ROC/MAST for cluster defining markers ####
  DefaultAssay(object = FSI_filtered) <- "RNA"
  FSI_filtered <- SetIdent(FSI_filtered, value = "RNA_snn_res.0.5")
  
  ##Find Markers that are specific for each cluster
  Batchedmarkers.mast_i_RNA_data_0.5=FindAllMarkers(FSI_filtered, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                    min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
  ##optional: add pct.fold = how large is the absolute difference in percentage?
  Batchedmarkers.mast_i_RNA_data_0.5$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.5$pct.1/Batchedmarkers.mast_i_RNA_data_0.5$pct.2
  
  ## Create list
  listDEgenes_i_RNA_MAST_0.5<- split(Batchedmarkers.mast_i_RNA_data_0.5, f=Batchedmarkers.mast_i_RNA_data_0.5$cluster)
  
  ## Filter on adj.P-value
  ##change name according to test used (MAST, roc, negbinom, et.c)
  listDEgenes_i_RNA_MAST_0.5 <-lapply(listDEgenes_i_RNA_MAST_0.5, function(x){dplyr::filter(x, p_val_adj<0.05)})
  ## Sort on logFC
  listDEgenes_i_RNA_MAST_0.5<-lapply(listDEgenes_i_RNA_MAST_0.5,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})
  
  
  ##save as Robj
  setwd("C:/Users/elise/Documents/PhD/Fetal 10X/DEG")
  save(listDEgenes_i_RNA_MAST_0.5, file='SupplData3_listDEG_fetal_FSI_scRNA_0.5.Robj')
  
  ## Write to Excel
  library('openxlsx')
  write.xlsx(listDEgenes_i_RNA_MAST_0.5, file='SupplData3_listDEG_fetal_FSI_scRNA_0.5.xlsx')
  detach("package:openxlsx", unload=TRUE)
  
  
  
#### 10. subset naive and memory CD4 T cells - with Treg ----------------------------------------------------------------------------------
FSI_filtered <- SetIdent(FSI_filtered, value = "RNA_snn_res.0.5")
FSI_CD4 <- subset(FSI_filtered,  idents= c('0','1','3','4','5','6'), subset=CD8_ADT<1)
rm(FSI_filtered)

## recluster
FSI_CD4 <- FindVariableFeatures(FSI_CD4)
VariableFeatures(FSI_CD4) <- VariableFeatures(FSI_CD4)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FSI_CD4))]
FSI_CD4 <- ScaleData(FSI_CD4)
FSI_CD4 <- RunPCA(FSI_CD4)
FSI_CD4 <- RunUMAP(FSI_CD4, dims = 1:30)
FSI_CD4 <- FindNeighbors(object = FSI_CD4, reduction = "pca", dims = 1:30)
FSI_CD4 <- FindClusters(object = FSI_CD4, resolution = seq(0.1:1, by=0.1))

DimPlot(object = FSI_CD4, reduction = "umap", pt.size=.5, group.by = 'RNA_snn_res.0.6')
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=.5, group.by = 'FSI_res.0.5_clusters')
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=.5, group.by = 'Donor')

clustree(FSI_CD4)

FeaturePlot(FSI_CD4, 'CD4.1',reduction='umap', pt.size=.2)
FeaturePlot(FSI_CD4, 'CD8a',reduction='umap', pt.size=.2)
FeaturePlot(FSI_CD4, 'CD45RA',reduction='umap', pt.size=1)
FeaturePlot(FSI_CD4, 'CD27.1',reduction='umap', pt.size=.2)
DimPlot(FSI_CD4, reduction='umap', group.by = 'Phase')

#### 10. CD4 - Figures ####
setwd("")

pdf('FSI_CD4_SurfaceData_res.0.6.pdf', width = 18, height=18)
FeaturePlot(object = SetIdent(FSI_CD4, value = "RNA_snn_res.0.6"), reduction = "umap", pt.size=.2,
            features =FSI_CD4@assays[["ADT"]]@counts@Dimnames[[1]], 
            label = T, label.size = 14) + plot_layout(ncol=3) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
        text=element_text(size=10))
dev.off()

pdf('FSI_CD4_SurfaceData_Vln_res.0.6.pdf', width = 18, height=18)
VlnPlot(object = SetIdent(FSI_CD4, value = "RNA_snn_res.0.6"),  
        features = FSI_CD4@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
        cols=c('purple','lightskyblue',"#CC79A7", "#E69F00","#009E73" ,"#F0E442" ,"black"),
        group.by = 'RNA_snn_res.0.6') + plot_layout(ncol=3)& 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()


pdf('FSI_CD4_Clustering_res.0.6.pdf', height=12, width=14)
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.6', label = T, label.size = 22,label.color = 'grey50',
        cols=c('purple','lightskyblue',"#CC79A7", "#E69F00","#009E73" ,"#F0E442" ,"black")) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Intestine', color='Cluster ID')+
  scale_color_manual(values=c('purple','lightskyblue',"#CC79A7", "#E69F00","#009E73" ,"#F0E442" ,"black"),
                     labels=c('0: Activated', '1: Naive-like/T-CM - activation/tissue', 
                              '2: Th1 - 1', '3: Th1 - 2', 
                              '4: CD4 Th2/early development', '5: Th17', '6: Naive-like/T-CM - quiescence'))
dev.off()

pdf('FSI_CD4_Clustering_res.0.5.pdf', height=12, width=14)
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.5', label = T, label.size = 22,label.color = 'grey50',
        cols=c( "#009E73" ,'black','purple',"#CC79A7", "#E69F00","#F0E442" )) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Intestine', color='Cluster ID')
dev.off()

pdf('FSI_CD4_Clustering_res.0.8.pdf', height=12, width=14)
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.8', label = T, label.size = 22,label.color = 'grey50',
        cols=c("#CC79A7", "#E69F00",'lightskyblue',"#009E73" ,"#F0E442","purple", 'purple4', "black", "#D55E00")) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Intestine', color='Cluster ID')+
  scale_color_manual(values=c("#CC79A7", "#E69F00",'lightskyblue',"#009E73" ,"#F0E442","purple", 'purple4', "black", "#D55E00"),
                     labels=c('0: Th1 - 1', '1: Th1 - 2',
                              '2: Naive-like/T-CM - activation/tissue', '3: CD4 Th2/early development',
                              '4: Th17','5: Activated - 1', '6: Activated - 2',
                              '7: Naive-like/T-CM - quiescence', '8: Antiviral'))
dev.off()

pdf('FSI_CD4_Oldclusters_totalFSIres0.5.pdf', height=12, width=14)
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=1, group.by = 'FSI_res.0.5_clusters', 
        cols=c("#E69F00","#F0E442","#009E73" ,'purple','black')) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  labs(title='Fetal Intestine', color='Cluster ID')
dev.off()

pdf('FSI_CD4_Donor.pdf', height=12, width=14)
DimPlot(object = FSI_CD4, reduction = "umap", pt.size=1, group.by = 'Donor') +
  ggtitle('Donor ID') +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=28),legend.title = element_text(size=28),
        plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))+
  scale_color_manual(values=palette.colors(palette = "R4"))
dev.off()

pdf('FSI_CD4_QCData_res.0.6.pdf', width = 18, height=10)
FeaturePlot(object = SetIdent(FSI_CD4, value = "RNA_snn_res.0.6"), reduction = "umap", pt.size=.4,
            features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
            label = T, label.size = 10) + DimPlot(FSI_CD4, reduction='umap', pt.size=.4, group.by ='Phase')+
  plot_layout(ncol=3) & 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=22))
dev.off()

pdf('FSI_CD4_QCData_Vln_res.0.6.pdf', width = 14, height=5)
VlnPlot(object = SetIdent(FSI_CD4, value = "RNA_snn_res.0.6"),  
        features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
        cols=c('purple','lightskyblue',"#CC79A7", "#E69F00","#009E73" ,"#F0E442" ,"black")) +
  plot_layout(nrow=1, ncol=4) & 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

FSI_CD4$CD45RA_ADT <- FSI_CD4@assays$ADT$data[3,]

pdf('FSI_CD4_Clustree_colouredbyCD45RA_ADT.pdf')
clustree(FSI_CD4, node_colour ='CD45RA_ADT', node_colour_aggr = 'median')
dev.off()

FSI_CD4@meta.data$FSI_CD4_res.0.6_clusters <- as.factor(FSI_CD4@meta.data$RNA_snn_res.0.6)
levels(FSI_CD4@meta.data$FSI_CD4_res.0.6_clusters) <- c('0: Activated', '1: Naive-like/T-CM - activation/tissue', 
                                                        '2: Th1 - 1', '3: Th1 - 2', 
                                                        '4: Th2/early development', '5: Th17', '6: Naive-like/T-CM - quiescence')
FSI_CD4@meta.data$FSI_CD4_res.0.6_clusters <- factor(FSI_CD4@meta.data$FSI_CD4_res.0.6_clusters,
                                                     levels=c('1: Naive-like/T-CM - activation/tissue', '6: Naive-like/T-CM - quiescence',
                                                              '4: Th2/early development', '0: Activated', 
                                                              '2: Th1 - 1', '3: Th1 - 2', 
                                                              '5: Th17'))
pdf('FSI_CD4_Proportionplot_total_res.0.6.pdf')
ggplot(FSI_CD4@meta.data, aes(x=orig.ident, fill=FSI_CD4_res.0.6_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  geom_text(stat = 'count', color='grey30',
            position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
  labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal intestine - CD4 subset')+
  theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('lightskyblue',"black", 'purple',"#009E73" ,"#CC79A7","#E69F00","#F0E442" ))
dev.off()

FSI_CD4@meta.data$FSI_CD4_res.0.8_clusters <- as.factor(FSI_CD4@meta.data$RNA_snn_res.0.8)
levels(FSI_CD4@meta.data$FSI_CD4_res.0.8_clusters) <-c('0: Th1 - 1', '1: Th1 - 2',
                                                       '2: Naive-like/T-CM - activation/tissue', '3: Th2/early development',
                                                       '4: Th17','5: Activated - 1', '6: Activated - 2',
                                                       '7: Naive-like/T-CM - quiescence', '8: Antiviral')
FSI_CD4@meta.data$FSI_CD4_res.0.8_clusters <- factor(FSI_CD4@meta.data$FSI_CD4_res.0.8_clusters,
                                                     levels=c('2: Naive-like/T-CM - activation/tissue','7: Naive-like/T-CM - quiescence',
                                                              '3: Th2/early development',
                                                              '5: Activated - 1','6: Activated - 2',
                                                              '0: Th1 - 1', '1: Th1 - 2', '8: Antiviral',
                                                              '4: Th17'))
pdf('FSI_CD4_Proportionplot_total_res.0.8.pdf')
ggplot(FSI_CD4@meta.data, aes(x=orig.ident, fill=FSI_CD4_res.0.8_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  geom_text(stat = 'count', color='grey40',
            position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
  labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal intestine - CD4 subset')+
  theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('lightskyblue',"black","#009E73" ,"purple", 'purple4',"#CC79A7","#E69F00", "#D55E00","#F0E442"))
dev.off()


pdf('FSI_CD4_Proportionplot_ClustersperDonor_res.0.6.pdf')
ggplot(FSI_CD4@meta.data, aes(x=Donor, fill=FSI_CD4_res.0.6_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='', subtitle='')+
  theme(plot.title = element_text(hjust=0.5, vjust=12,size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('lightskyblue',"black", 'purple',"#009E73" ,"#CC79A7","#E69F00","#F0E442" ))+
  coord_flip()
dev.off()

pdf('FSI_CD4_Proportionplot_ClustersperDonor_res.0.8.pdf')
ggplot(FSI_CD4@meta.data, aes(x=Donor, fill=FSI_CD4_res.0.8_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='', subtitle='')+
  theme(plot.title = element_text(hjust=0.5, vjust=12,size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c('lightskyblue',"black","#009E73" ,"purple", 'purple4',"#CC79A7","#E69F00", "#D55E00","#F0E442"))+
  coord_flip()
dev.off()


pdf('FSI_CD4_Proportionplot_DonorsperCluster_res.0.6.pdf')
ggplot(FSI_CD4@meta.data, aes(fill=Donor, x=RNA_snn_res.0.6)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  labs(fill = "Donor ID", title='', subtitle='')+
  theme(plot.title = element_text(hjust=0.5, vjust=12,size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=palette.colors(palette = "R4"))
dev.off()

pdf('FSI_CD4_Proportionplot_DonorsperCluster_res.0.8.pdf')
ggplot(FSI_CD4@meta.data, aes(fill=Donor, x=RNA_snn_res.0.8)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  labs(fill = "Donor ID", title='', subtitle='')+
  theme(plot.title = element_text(hjust=0.5, vjust=12,size=15), plot.subtitle = element_text(hjust=0.5, vjust=14),
        axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=palette.colors(palette = "R4"))
dev.off()


#### 10. CD4 - DEG ####
DefaultAssay(object = FSI_CD4) <- "RNA"
FSI_CD4 <- SetIdent(FSI_CD4, value = "RNA_snn_res.0.3")

##Find Markers that are specific for each cluster
Batchedmarkers.mast_i_RNA_data_0.3=FindAllMarkers(FSI_CD4, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_0.3$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.3$pct.1/Batchedmarkers.mast_i_RNA_data_0.3$pct.2

## Create list
listDEgenes_i_RNA_MAST_0.3<- split(Batchedmarkers.mast_i_RNA_data_0.3, f=Batchedmarkers.mast_i_RNA_data_0.3$cluster)

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
listDEgenes_i_RNA_MAST_0.3 <-lapply(listDEgenes_i_RNA_MAST_0.3, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
listDEgenes_i_RNA_MAST_0.3<-lapply(listDEgenes_i_RNA_MAST_0.3,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

