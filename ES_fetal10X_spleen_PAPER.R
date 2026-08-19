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

### 1. Setting working directory and loading data----------------------------------------------------------------
setwd("~/PhD/Fetal 10X/sample_filtered_feature_bc_matrix_FSP/sample_filtered_feature_bc_matrix")

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
FSP_data <- Read10X(data.dir = getwd()) 

### 2. Seurat object  -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data
FSP <- CreateSeuratObject(counts = FSP_data$`Gene Expression`, project = "FSP")
FSP

FSP_ADT <- CreateAssayObject(counts=FSP_data$`Antibody Capture`)
FSP[['ADT']] <- FSP_ADT

rownames(FSP[['ADT']])
rownames(FSP[['RNA']])

rm(FSP_data)
rm(FSP_ADT)

### 3. Add info to metadata ----------------------------------------------------------------------------------
seurat_G <- readRDS("~/all_ab_old_harmony_standarized.rds")

colnames(seurat_G@meta.data)
rownames(seurat_G@meta.data)

## extract only the FSP data and remove FSP from rownames
FSP_seurat_G <- subset(seurat_G, subset=orig.ident=='FSP_ab') 
rownames(FSP_seurat_G@meta.data) <- str_remove(rownames(FSP_seurat_G@meta.data), 'FSP_')
rownames(FSP_seurat_G@meta.data)

# donor info = 'donor_id_test' column
FSP@meta.data <- cbind(FSP@meta.data, 
                       FSP_seurat_G@meta.data[match(rownames(FSP@meta.data),
                                                    rownames(FSP_seurat_G@meta.data)),
                                              c(12)])
colnames(FSP@meta.data)[6] <- 'Donor'


### 4. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
setwd("~/")
DefaultAssay(FSP) <- 'RNA'
  
## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  FSP <- PercentageFeatureSet(FSP, pattern = "^MT-", col.name = 'percent.mt')
  FSP <- PercentageFeatureSet(FSP, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(FSP@assays$RNA@counts)[grep('MT', rownames(FSP@assays$RNA@counts))]) #somewhere in name
  length(rownames(FSP@assays$RNA@counts)[grep('^MT\\.', rownames(FSP@assays$RNA@counts))]) #start of name followed by .
  length(rownames(FSP@assays$RNA@counts)[grep('^MT-', rownames(FSP@assays$RNA@counts))]) #start of name followed by -
  length(rownames(FSP@assays$RNA@counts)[grep('^RP', rownames(FSP@assays$RNA@counts))]) #start of name followed by -
  
  
  # Show QC metrics for the first 5 cells
  head(FSP@meta.data, 5)
  
  
  # Visualize QC metrics as a violin plot
  pdf(file="FSP_beforeQC_Vlns.pdf", width = 20, height = 8)
  VlnPlot(FSP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="FSP_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(FSP, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)
  plot1.z <- FeatureScatter(FSP, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.00001)+
    scale_x_continuous(limits=c(0,5000))
  plot2 <- FeatureScatter(FSP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.00001)
  plot3 <- FeatureScatter(FSP, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size = 0.00001)
  plot1 + plot1.z+  plot2 + plot3 + plot_layout(ncol=4)
  dev.off()
  
  #With filtering cut-offs
  pdf(file="FSP_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
  plot1 <- ggplot(FSP@meta.data, aes(FSP$nCount_RNA, FSP$nFeature_RNA))+geom_point(size=.5)+theme_bw()+
    geom_vline(xintercept = 32000)+ annotate("text",x=33000,y=1,label=c("32000"),hjust=0, size=2.8)+
    geom_vline(xintercept = 1200)+ annotate("text",x=1100,y=1,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 750)+ annotate("text",x=600,y=1,label=c("750"),hjust=0, size=2.8)
  plot2 <- ggplot(FSP@meta.data, aes(FSP$nCount_RNA, FSP$percent.mt))+geom_point(size=.5)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=15,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=3000,y=8,label=c("7"),vjust=0, size=2.8)+
    xlim(0,5000)
  plot3 <- ggplot(FSP@meta.data, aes(FSP$nCount_RNA, FSP$percent.mt))+geom_point(size=.5)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=20,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=10000,y=9,label=c("7"),vjust=0, size=2.8)
  plot4 <- ggplot(FSP@meta.data, aes(FSP$nFeature_RNA, FSP$percent.mt))+geom_point(size=.5)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 750)+ annotate("text",x=900,y=20,label=c("750"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=10000,y=9,label=c("7"),vjust=0, size=2.8)+
    xlim(0,1000)
  plot5 <- ggplot(FSP@meta.data, aes(FSP$nCount_RNA, FSP$percent.ribo))+geom_point(size=.5)+theme_bw()+
    geom_vline(xintercept = 1200)+ annotate("text",x=1400,y=50,label=c("1200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 7)+ annotate("text",x=50,y=10,label=c("7"),hjust=0, size=2.8)+xlim(0,5000)
  plot6 <- ggplot(FSP@meta.data, aes(FSP$percent.mt, FSP$percent.ribo))+geom_point(size=.5)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5+plot6
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(FSP, expression = nCount_RNA < 32000 & nFeature_RNA > 750 & 
                           nCount_RNA > 1200 &percent.mt < 7 &percent.ribo >7)
  length(selected) #How many cells are left
  
  # Filter cells
  FSP_filtered <- subset(FSP, subset = nCount_RNA < 32000 & nFeature_RNA > 750 & 
                           nCount_RNA > 1200 &percent.mt < 7 &percent.ribo >7)
  
  # Visualize QC metrics post-filtering
  pdf(file="FSP_afterQC_Vlns.pdf", width = 20, height = 8)
  VlnPlot(FSP_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 5, pt.size = 0.000000001)
  dev.off()
  
  
## remove doublets
  FSP_filtered$doublets <- ifelse(FSP_filtered$Donor=='doublet', T, F)
  FSP_filtered@meta.data[is.na(FSP_filtered$doublets),'doublets'] <- F
  FSP_filtered <- subset(FSP_filtered, subset= doublets==F)
  
## unassigned donor = NA
  FSP_filtered@meta.data[grep('unassigned', FSP_filtered$Donor), 'Donor'] <- NA
 
### 5. Normalization and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  FSP_filtered <- NormalizeData(FSP_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  FSP_filtered <- CellCycleScoring(object = FSP_filtered, s.features = s.genes, g2m.features = g2m.genes, 
                                   set.ident = FALSE)
  head(x = FSP_filtered@meta.data)
  
  FSP_filtered <- ScaleData(FSP_filtered)
  FSP_filtered <- FindVariableFeatures(FSP_filtered)
  
  ## remove TCR genes from variable features 
  VariableFeatures(FSP_filtered) <- VariableFeatures(FSP_filtered)[!grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FSP_filtered))]
  
  FSP_filtered <- RunPCA(FSP_filtered)

### 6. ADT data ------------------------------------------------------------------------------
  DefaultAssay(FSP_filtered) <- 'ADT'
  FSP_filtered <- NormalizeData(FSP_filtered, normalization.method = 'CLR', margin = 2)
  
  rownames(FSP_filtered@assays$ADT$data)
  FSP_filtered$CD4_ADT <- FSP_filtered@assays$ADT$data[6,]
  FSP_filtered$CD8_ADT <- FSP_filtered@assays$ADT$data[7,]
  FSP_filtered$CD45RA_ADT <- FSP_filtered@assays$ADT$data[3,]
  
  DefaultAssay(FSP_filtered) <- 'RNA'

#### Save and load -----------------------------------------------------------------------------------------------------------
  
  setwd("~/PhD/Fetal 10X/UsedObjects")
  save(FSP_filtered, file = 'FSP_filtered_simple.Robj')
  save(FSP_CD4, file = 'FSP_CD4subset_wTreg.Robj')
  
  ## load
  setwd("~/PhD/Fetal 10X/UsedObjects")
  load('FSP_filtered_simple.Robj')
  load('FSP_CD4subset_wTreg.Robj')
  
### 7. Clustering ----------------------------------------------------------------------------------------------------------------------------
  DefaultAssay(FSP_filtered) <- 'RNA'
  FSP_filtered <- RunUMAP(FSP_filtered, dims = 1:30)
  
  FSP_filtered <- FindNeighbors(object = FSP_filtered, reduction = "pca", dims = 1:30)
  FSP_filtered <- FindClusters(object = FSP_filtered, resolution = seq(0.1:1, by=0.1))
  
  DimPlot(object = FSP_filtered, reduction = "umap", pt.size=.2, group.by = 'RNA_snn_res.0.4')
  
  
  clustree(FSP_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')
  
  FeaturePlot(FSP_filtered, 'CD4.1',reduction='umap', pt.size=.2)
  FeaturePlot(FSP_filtered, 'CD8a',reduction='umap', pt.size=.2)
  FeaturePlot(FSP_filtered, 'CD45RA',reduction='umap', pt.size=.2)
  FeaturePlot(FSP_filtered, 'CD27.1',reduction='umap', pt.size=.2)
  DimPlot(FSP_filtered, reduction='umap', group.by = 'Phase')
  
  ### 8. Figures ---------------------------------------------------------------------------------------------------
  setwd("")
  donor_colors <- c(palette.colors(palette = "R4")[1:4],palette.colors(palette = "R4")[6:7])
  
  FSP_filtered@meta.data$FSP_res.0.4_clusters <- as.factor(FSP_filtered@meta.data$RNA_snn_res.0.4)
  levels(FSP_filtered@meta.data$FSP_res.0.4_clusters) <- c('0: CD8 Naive-like', '1: CD4 Naive/TCM 1', 
                                                           '2: CD4 Naive/TCM 2 - stressed','3: CD8 effector/memory', 
                                                           '4: Treg', '5: Proliferating - S/G2M Phase', '6: CD4/CD8 early memory/Th17',
                                                           '7: CD4/CD8 activated/Th1', 
                                                           '8: B cell contamination')
  FSP_filtered@meta.data$FSP_res.0.4_clusters <- factor(FSP_filtered@meta.data$FSP_res.0.4_clusters,
                                                        levels=c( '8: B cell contamination', '5: Proliferating - S/G2M Phase', 
                                                                  '0: CD8 Naive-like', '3: CD8 effector/memory',
                                                                  '7: CD4/CD8 activated/Th1', '6: CD4/CD8 early memory/Th17',
                                                                  '1: CD4 Naive/TCM 1',  '2: CD4 Naive/TCM 2 - stressed',
                                                                  '4: Treg'))
  
  FigS1C <- ((FeaturePlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.5,
                          label.size = 12,
                          features = c('CD4.1'), cols = c('gold', 'purple'),
                          label = T)+ggtitle('CD4 surface expression')) + 
               (FeaturePlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.5,
                            label.size = 12,
                            features = c('CD8a'), cols = c('gold', 'purple'),
                            label = T)+ggtitle('CD8A surface expression'))+
               DimPlot(FSP_filtered, reduction='umap', pt.size=.5, group.by ='Phase', cols=c('grey20', 'violet', 'seagreen1')) + 
               DimPlot(FSP_filtered, reduction='umap', pt.size=.5, group.by ='Donor', cols=donor_colors)  &
               theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                     legend.text = element_text(size=26),legend.title = element_text(size=28),
                     legend.key.size = unit(20,'points'),
                     plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=2),plot.subtitle = element_text(size=28, hjust=0.5),
                     plot.margin = margin(20,20,20,20),
                     text=element_text(size=8))) +
    (VlnPlot(object = SetIdent(FSP_filtered, value = "FSP_res.0.4_clusters"),  
             features = c('CD45RA'), pt.size=0,
             cols=c('grey','red4','mistyrose', 'plum3',
                    'cadetblue2','deepskyblue3','blue', '#669999','turquoise')) +
       scale_x_discrete(labels=c('8','5','0','3','7','6','1','2','4'))+
       ggtitle('CD45RA surface expression')+
       theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=22), 
             axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
             legend.text = element_text(size=14),
             plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),
             plot.margin = margin(10,10,10,10),
             text=element_text(size=14))+
       xlab('Cluster ID')+
       guides(fill='none'))+
    (ggplot(FSP_filtered@meta.data, aes(x=orig.ident, fill=FSP_res.0.4_clusters)) + theme_classic() +
       geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
       geom_text(stat = 'count', size=6,
                 position = position_fill(vjust=.5), aes(color=FSP_res.0.4_clusters,
                                                         label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
       labs(fill = "ClusterID", title='Proportion of T cell clusters')+
       theme(axis.text.x = element_blank(), axis.title.x = element_text(size=22), 
             axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
             legend.text = element_text(size=24),legend.title = element_text(size=26),
             plot.title = element_text(size=30, face='bold', hjust=0.5, vjust=3),plot.subtitle = element_text(size=28, hjust=0.5),
             plot.margin = margin(40,40,40,40),
             text=element_text(size=8))+
       scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
       scale_fill_manual(values=c('grey','red4','mistyrose', 'plum3',
                                  'cadetblue2','deepskyblue3','blue', '#669999','turquoise'))+
       scale_color_manual(values=c('grey10','grey80', rep('grey20',4), 'grey60','grey10', 'grey10'))+
       guides(color='none'))+ 
    plot_layout(ncol=2)+
    plot_annotation(title = "Fetal spleen", 
                    theme=theme(plot.title=element_text(size=50, hjust=0.35, vjust=4), 
                                plot.margin=margin(60,30,10,30)))
  FigS1C
  ggsave(file='FSP_overview_res.0.4_SupplFig1C.pdf', FigS1C, width=20, height=18)
  
  
  pdf('FSP_SurfaceData_res.0.4.pdf', width = 18, height=18)
  FeaturePlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.2,
              features =FSP_filtered@assays[["ADT"]]@counts@Dimnames[[1]], 
              label = T, label.size = 14) + plot_layout(ncol=3) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
          text=element_text(size=10))
  dev.off()
  
  pdf('FSP_SurfaceData_Vln_res.0.4.pdf', width = 18, height=18)
  VlnPlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"),  
          features = FSP_filtered@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
          cols=c('mistyrose',  'black','grey','plum3','turquoise','red4',"gold",'coral','#669999'),
          group.by = 'RNA_snn_res.0.4') + plot_layout(ncol=3)& 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  pdf('FSP_Clustering_res.0.4.pdf', height=12, width=12)
  DimPlot(object = FSP_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.4', 
          label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
          cols=c('mistyrose', 'black','grey','plum3','turquoise','red4',"gold",'coral','#669999')) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen', color='Cluster ID')+
    scale_color_manual(values=c('mistyrose', 'black','grey','plum3','turquoise','red4',"gold",'coral','#669999'))
  dev.off()
  
  pdf('FSP_Clustering_res.0.4_clusternames.pdf', height=12, width=18)
  DimPlot(object = FSP_filtered, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.4', 
          label = T, label.size = 22,label.color = 'grey50', repel=T,label.box = T,
          cols=c('mistyrose',  'black','grey','plum3','turquoise','red4',"gold",'coral','#669999')) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen', color='Cluster ID')+
    scale_color_manual(values=c('mistyrose',  'black','grey','plum3','turquoise','red4',"gold",'coral','#669999'),
                       labels= c('0: CD8 Naive-like', '1: CD4 Naive/TCM 1', 
                                 '2: CD4 Naive/TCM 2 - stressed','3: CD8 effector/memory', 
                                 '4: Treg', '5: Proliferating - S/G2M Phase', '6: CD4/CD8 early development/Th17',
                                 '7: CD4/CD8 activated/Th1', 
                                 '8: B cell contamination'))
  dev.off()
  
  pdf('FSP_Donor.pdf', height=12, width=12)
  DimPlot(object = FSP_filtered, reduction = "umap", pt.size=1, group.by = 'Donor') +
    ggtitle('Donor ID') +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    scale_color_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  pdf('FSP_QCData_res.0.4.pdf', width = 18, height=10)
  FeaturePlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.4,
              features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
              label = T, label.size = 10) + DimPlot(FSP_filtered, reduction='umap', pt.size=.4, group.by ='Phase')+
    plot_layout(ncol=3) & 
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()
  
  pdf('FSP_QCData_Vln_res.0.4.pdf', width = 14, height=5)
  VlnPlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"),  
          features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
          cols=c('mistyrose',  'black','grey','plum3','turquoise','red4',"gold",'coral','#669999')) +
    plot_layout(nrow=1, ncol=4) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  pdf('FSP_overview_res.0.4.pdf', width = 18, height=16)
  FeaturePlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.5,
              label.size = 12,
              features = c('CD4.1', 'CD8a', 'CD45RA'),
              label = T) + DimPlot(FSP_filtered, reduction='umap', pt.size=.5, group.by ='Phase') + plot_layout(ncol=2) &
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()
  
  pdf('FSP_CD45RA_res.0.4.pdf', width = 6, height=6)
  FeaturePlot(object = SetIdent(FSP_filtered, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.5,
              label.size = 8,
              features = c('CD45RA'),
              label = T,label.color = c('black','blue','black','red4','blue','turquoise3','blue','blue','grey'))+
    scale_color_gradient(limits=c(0,4),low="lightgrey", high= "blue")+
    ggtitle('CD45RA surface expression')+
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()
  
  
  pdf('FSP_Clustree_colouredbyCD4ADT.pdf')
  clustree(FSP_filtered, node_colour ='CD4_ADT', node_colour_aggr = 'median') + scale_color_gradient(low = 'purple', high='gold')
  dev.off()
  
  
  FSP_filtered@meta.data$FSP_res.0.4_clusters <- as.factor(FSP_filtered@meta.data$RNA_snn_res.0.4)
  levels(FSP_filtered@meta.data$FSP_res.0.4_clusters) <- c('0: CD8 Naive-like', '1: CD4 Naive/TCM 1', 
                                                           '2: CD4 Naive/TCM 2 - stressed','3: CD8 effector/memory', 
                                                           '4: Treg', '5: Proliferating - S/G2M Phase', '6: CD4/CD8 early development/Th17',
                                                           '7: CD4/CD8 activated/Th1', 
                                                           '8: B cell contamination')
  FSP_filtered@meta.data$FSP_res.0.4_clusters <- factor(FSP_filtered@meta.data$FSP_res.0.4_clusters,
                                                        levels=c( '8: B cell contamination', '5: Proliferating - S/G2M Phase', 
                                                                 '0: CD8 Naive-like', '3: CD8 effector/memory',
                                                                 '7: CD4/CD8 activated/Th1', '6: CD4/CD8 early development/Th17',
                                                                 '1: CD4 Naive/TCM 1',  '2: CD4 Naive/TCM 2 - stressed',
                                                                 '4: Treg'))
  pdf('FSP_Proportionplot_total_res.0.4.pdf')
  ggplot(FSP_filtered@meta.data, aes(x=orig.ident, fill=FSP_res.0.4_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    geom_text(stat = 'count', col='grey40',size=4,
              position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
    labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal Spleen')+
    theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('#669999','red4','mistyrose', 'plum3',
                               'coral',"gold",  'black','grey','turquoise'))
  dev.off()
  
  pdf('FSP_Proportionplot_total_res.0.4_simplecolors.pdf', height=9, width=10)
  ggplot(FSP_filtered@meta.data, aes(x=orig.ident, fill=FSP_res.0.4_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    geom_text(stat = 'count', size=5,
              position = position_fill(vjust=.5), aes(color=FSP_res.0.4_clusters,
                                                      label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
    labs(fill = "ClusterID", title='Fetal spleen', subtitle='Proportion of T cell clusters')+
    theme(axis.text.y = element_text(size=16), axis.title.y = element_text(size=18,vjust=3), 
          axis.text.x = element_text(size=20),
          legend.text = element_text(size=16), legend.title = element_text(size=18),
          plot.title = element_text(size=28, hjust=0.5, , vjust=4), plot.subtitle = element_text(size=24, hjust=0.5, , vjust=4),
          plot.margin = margin(30,30,30,30),
          text=element_text(size=16))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('grey','red4','mistyrose', 'plum3',
                               'cadetblue2','deepskyblue3','blue', '#669999','turquoise'))+
    scale_color_manual(values=c('grey10','grey80', rep('grey20',4), 'grey60','grey10', 'grey10'))+
    guides(color='none')
  dev.off()
  
  pdf('FSP_Proportionplot_ClustersperDonor_res.0.4.pdf', width=8,height=5)
  ggplot(FSP_filtered@meta.data, aes(x=Donor, fill=FSP_res.0.4_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    labs(fill = "Cluster", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('#669999','red4','mistyrose', 'plum3',
                               'coral',"gold",  'black','grey','turquoise'))+
    coord_flip()
  dev.off()
  
  pdf('FSP_Proportionplot_DonorsperCluster_res.0.4.pdf')
  ggplot(FSP_filtered@meta.data, aes(fill=Donor, x=RNA_snn_res.0.4)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "Donor ID", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  pdf('FSP_CD45RA_Vln_res.0.4.pdf', width = 12, height=6)
  VlnPlot(object = SetIdent(FSP_filtered, value = "FSP_res.0.4_clusters"),  
          features = c('CD45RA'), pt.size=0,
          cols=c('grey','red4','mistyrose', 'plum3',
                 'cadetblue2','deepskyblue3','blue', '#669999','turquoise')) +
    ggtitle('CD45RA surface expression')+
    theme(axis.text = element_text(size=18), axis.title = element_text(size=20), legend.text = element_text(size=18),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10))+
    scale_x_discrete(labels=c('8','5','0','3','7','6','1','2','4'))
  dev.off()
  
  ### 9. DEG -----------------------------------------------------------------------------------------------------------
  #### ROC/MAST for cluster defining markers
  DefaultAssay(object = FSP_filtered) <- "RNA"
  FSP_filtered <- SetIdent(FSP_filtered, value = "RNA_snn_res.0.4")
  
  ##Find Markers that are specific for each cluster
  Batchedmarkers.mast_i_RNA_data_0.4=FindAllMarkers(FSP_filtered, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                    min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
  ##optional: add pct.fold = how large is the absolute difference in percentage?
  Batchedmarkers.mast_i_RNA_data_0.4$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.4$pct.1/Batchedmarkers.mast_i_RNA_data_0.4$pct.2
  
  ## Create list
  listDEgenes_i_RNA_MAST_0.4 <- split(Batchedmarkers.mast_i_RNA_data_0.4, f=Batchedmarkers.mast_i_RNA_data_0.4$cluster)
  
  ## Filter on adj.P-value
  ##change name according to test used (MAST, roc, negbinom, et.c)
  listDEgenes_i_RNA_MAST_0.4 <-lapply(listDEgenes_i_RNA_MAST_0.4, function(x){dplyr::filter(x, p_val_adj<0.05)})
  ## Sort on logFC
  listDEgenes_i_RNA_MAST_0.4 <-lapply(listDEgenes_i_RNA_MAST_0.4,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})
  
  
  ##save as Robj
  setwd("")
  save(listDEgenes_i_RNA_MAST_0.4, file='SupplData2_listDEG_fetal_FSP_scRNA_0.4.Robj')
  
  ## Write to Excel
  library('openxlsx')
  write.xlsx(listDEgenes_i_RNA_MAST_0.4, file='SupplData2_listDEG_fetal_FSP_scRNA_0.4.xlsx')
  detach("package:openxlsx", unload=TRUE)

  #### 10.2  subset naive and memory CD4 T cells - with Treg ----------------------------------------------------------------------------------
  FSP_filtered<- SetIdent(FSP_filtered, value = "RNA_snn_res.0.4")
  FSP_CD4 <- subset(FSP_filtered,  idents= c('1','2','4','6','7'), subset=CD8_ADT<1&CD4_ADT>1)
  rm(FSP_filtered)
  
  ## recluster
  FSP_CD4 <- FindVariableFeatures(FSP_CD4)
  VariableFeatures(FSP_CD4) <- VariableFeatures(FSP_CD4)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(FSP_CD4))]
  FSP_CD4 <- ScaleData(FSP_CD4)
  FSP_CD4 <- RunPCA(FSP_CD4)
  FSP_CD4 <- RunUMAP(FSP_CD4, dims = 1:30)
  FSP_CD4 <- FindNeighbors(object = FSP_CD4, reduction = "pca", dims = 1:30)
  FSP_CD4 <- FindClusters(object = FSP_CD4, resolution = seq(0.1:1, by=0.1))
  
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=.5, group.by = 'RNA_snn_res.0.3')+
    DimPlot(object = FSP_CD4, reduction = "umap", pt.size=.5, group.by = 'FSP_res.0.4_clusters')
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=.5, group.by = 'Donor')
  
  clustree(FSP_CD4)
  
  FeaturePlot(FSP_CD4, 'CD4.1',reduction='umap', pt.size=1)
  FeaturePlot(FSP_CD4, 'CD8a',reduction='umap', pt.size=1)
  FeaturePlot(FSP_CD4, 'CD45RA',reduction='umap', pt.size=1)
  FeaturePlot(FSP_CD4, 'CD27.1',reduction='umap', pt.size=.2)
  DimPlot(FSP_CD4, reduction='umap', group.by = 'Phase', pt.size=1)
  
  #### 10. CD4 - Figures ####
  setwd("")
  
  pdf('FSP_CD4_SurfaceData_res.0.2.pdf', width = 18, height=18)
  FeaturePlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.2,
              features =FSP_CD4@assays[["ADT"]]@counts@Dimnames[[1]], 
              label = T, label.size = 14) + plot_layout(ncol=3) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
          text=element_text(size=10))
  dev.off()
  
  pdf('FSP_CD4_SurfaceData_Vln_res.0.2.pdf', width = 18, height=18)
  VlnPlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.2"),  
          features = FSP_CD4@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
          cols=c( 'black','grey',"gold","coral"),
          group.by = 'RNA_snn_res.0.2') + plot_layout(ncol=3)& 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  
  pdf('FSP_CD4_SurfaceData_res.0.4.pdf', width = 18, height=18)
  FeaturePlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.4"), reduction = "umap", pt.size=.2,
              features =FSP_CD4@assays[["ADT"]]@counts@Dimnames[[1]], 
              label = T, label.size = 14) + plot_layout(ncol=3) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20), plot.margin = margin(10,10,10,10),
          text=element_text(size=10))
  dev.off()
  
  pdf('FSP_CD4_SurfaceData_Vln_res.0.4.pdf', width = 18, height=18)
  VlnPlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.4"),  
          features = FSP_CD4@assays[["ADT"]]@counts@Dimnames[[1]], pt.size=0,
          cols=c( 'blue','black','grey','antiquewhite4',"gold",'#669999',"coral"),
          group.by = 'RNA_snn_res.0.4') + plot_layout(ncol=3)& 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  pdf('FSP_CD4_Clustering_res.0.2.pdf', height=12, width=12)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.2',
          label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
          cols=c( 'black','grey',"gold","coral")) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen - CD4 subset', color='Cluster ID')
  dev.off()
  
  pdf('FSP_CD4_Clustering_res.0.2_clusternames.pdf', height=12, width=18)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.2', 
          label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
          cols=c( 'black','grey',"gold","coral")) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen - CD4 subset', color='Cluster ID')+
    scale_color_manual(values=c( 'black','grey',"gold","coral"),
                       labels=c('0: Naive/TCM 1', 
                                '1: Naive/TCM 2',
                                '2: early development/Th17', '3: activated/Th1'))
  dev.off()
  
  
  pdf('FSP_CD4_Clustering_res.0.4.pdf', height=12, width=12)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.4', 
          label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
          cols=c( 'blue','black','grey','antiquewhite4',"gold",'#669999',"coral" )) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen - CD4 subset', color='Cluster ID')
  dev.off()
  
  pdf('FSP_CD4_Clustering_res.0.4_clusternames.pdf', height=12, width=18)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.4', label = T, 
          label.size = 22,label.color = 'grey60', repel=T,label.box = T,
          cols=c( 'blue','black','grey','antiquewhite4',"gold",'#669999',"coral")) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen - CD4 subset', color='Cluster ID')+
    scale_color_manual(values=c('blue','black','grey','antiquewhite4',"gold",'#669999',"coral"),
                       labels=c('0: Naive/TCM 1', 
                                '1: Naive/TCM 2',
                                '2: Naive/TCM 3', '3: Naive/TCM - stressed',
                                '4: early development/Th17', '5: Naive/TCM - low quality',
                                '6: activated/Th1'))
  dev.off()
  
  pdf('FSP_CD4_Clustering_res.0.3.pdf', height=12, width=12)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'RNA_snn_res.0.3', 
          label = T, label.size = 22,label.color = 'grey60', repel=T,label.box = T,
          cols=c( 'blue','grey','black','antiquewhite4',"gold","coral")) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen - CD4 subset', color='Cluster ID')
  dev.off()
  
  
  pdf('FSP_CD4_Oldclusters_totalFSPres0.4.pdf', height=12, width=18)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'FSP_res.0.4_clusters', 
          cols=c("coral","gold", 'black','grey')) +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    labs(title='Fetal Spleen - CD4', color='Cluster ID - total FSP clusters')
  dev.off()
  
  pdf('FSP_CD4_Donor.pdf', height=12, width=12)
  DimPlot(object = FSP_CD4, reduction = "umap", pt.size=1, group.by = 'Donor') +
    ggtitle('Donor ID') +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=28),legend.title = element_text(size=28),
          plot.title = element_text(size=34, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))+
    scale_color_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  pdf('FSP_CD4_QCData_res.0.2.pdf', width = 18, height=10)
  FeaturePlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.2"), reduction = "umap", pt.size=.4,
              features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'),
              label = T, label.size = 10) + DimPlot(FSP_CD4, reduction='umap', pt.size=.4, group.by ='Phase')+
    plot_layout(ncol=3) & 
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          plot.title = element_text(size=28, hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
          plot.margin = margin(20,20,20,20),
          text=element_text(size=22))
  dev.off()
  
  pdf('FSP_CD4_QCData_Vln_res.0.2.pdf', width = 14, height=5)
  VlnPlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.2"),  
          features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
          cols=c('black','grey',"gold","coral")) +
    plot_layout(nrow=1, ncol=4) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  pdf('FSP_CD4_QCData_Vln_res.0.4.pdf', width = 14, height=5)
  VlnPlot(object = SetIdent(FSP_CD4, value = "RNA_snn_res.0.4"),  
          features = c('nFeature_RNA','nCount_RNA', 'percent.mt', 'percent.ribo'), pt.size=0,
          cols=c('blue','black','grey','antiquewhite4',"gold",'#669999',"coral")) +
    plot_layout(nrow=1, ncol=4) & 
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
          plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
          plot.margin = margin(10,10,10,10),
          text=element_text(size=14))
  dev.off()
  
  
  FSP_CD4$CD45RA_ADT <- FSP_CD4@assays$ADT$data[3,]
  
  pdf('FSP_CD4_Clustree_colouredbyCD45RA_ADT.pdf')
  clustree(FSP_CD4, node_colour ='CD45RA_ADT', node_colour_aggr = 'median')
  dev.off()
  
  
  FSP_CD4@meta.data$FSP_CD4_res.0.2_clusters <- as.factor(FSP_CD4@meta.data$RNA_snn_res.0.2)
  levels(FSP_CD4@meta.data$FSP_CD4_res.0.2_clusters) <- c('0: Naive/TCM 1', 
                                                          '1: Naive/TCM - stressed',
                                                          '2: early development/Th17', '3: activated/Th1')
  FSP_CD4@meta.data$FSP_CD4_res.0.2_clusters <- factor(FSP_CD4@meta.data$FSP_CD4_res.0.2_clusters,
                                                       levels=c('0: Naive/TCM 1', 
                                                                '1: Naive/TCM - stressed',
                                                                '2: early development/Th17', '3: activated/Th1'))
  pdf('FSP_CD4_Proportionplot_total_res.0.2.pdf')
  ggplot(FSP_CD4@meta.data, aes(x=orig.ident, fill=FSP_CD4_res.0.2_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    geom_text(stat = 'count', color='grey30',
              position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
    labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal Spleen - CD4 subset')+
    theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('black','grey',"gold","coral" ))
  dev.off()
  
  FSP_CD4@meta.data$FSP_CD4_res.0.4_clusters <- as.factor(FSP_CD4@meta.data$RNA_snn_res.0.4)
  levels(FSP_CD4@meta.data$FSP_CD4_res.0.4_clusters) <- c('0: Naive/TCM 1', 
                                                          '1: Naive/TCM 2',
                                                          '2: Naive/TCM 3', '3: Naive/TCM - stressed',
                                                          '4: early development/Th17', '5: Naive/TCM - low quality',
                                                          '6: activated/Th1')
  FSP_CD4@meta.data$FSP_CD4_res.0.4_clusters <- factor(FSP_CD4@meta.data$FSP_CD4_res.0.4_clusters,
                                                       levels=c('0: Naive/TCM 1', 
                                                                '1: Naive/TCM 2',
                                                                '2: Naive/TCM 3', 
                                                                '3: Naive/TCM - stressed','5: Naive/TCM - low quality',
                                                                '4: early development/Th17', 
                                                                '6: activated/Th1'))
  pdf('FSP_CD4_Proportionplot_total_res.0.4.pdf')
  ggplot(FSP_CD4@meta.data, aes(x=orig.ident, fill=FSP_CD4_res.0.4_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    geom_text(stat = 'count', color='grey30',
              position = position_fill(vjust=.5), aes(label = scales::percent(accuracy=0.1,prop.table(stat(count))))) +
    labs(fill = "Cluster", title='Proportion of T cell clusters', subtitle='Fetal Spleen - CD4 subset')+
    theme(plot.title = element_text(hjust=0.5, vjust=6,size=16), plot.subtitle = element_text(hjust=0.5, vjust=6, size=14),
          axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c( 'blue','black','grey','antiquewhite4','#669999',"gold","coral"))
  dev.off()
  
  
  
  pdf('FSP_CD4_Proportionplot_ClustersperDonor_res.0.2.pdf', width=8, height=5)
  ggplot(FSP_CD4@meta.data, aes(x=Donor, fill=FSP_CD4_res.0.2_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    labs(fill = "Cluster", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('black','grey',"gold","coral"  ))+
    coord_flip()
  dev.off()
  
  pdf('FSP_CD4_Proportionplot_ClustersperDonor_res.0.4.pdf', width=8, height=5)
  ggplot(FSP_CD4@meta.data, aes(x=Donor, fill=FSP_CD4_res.0.4_clusters)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    labs(fill = "Cluster", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=c('blue','#669999','black','grey','antiquewhite4',"gold","coral"))+
    coord_flip()
  dev.off()
  
  
  pdf('FSP_CD4_Proportionplot_DonorsperCluster_res.0.2.pdf')
  ggplot(FSP_CD4@meta.data, aes(fill=Donor, x=RNA_snn_res.0.2)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "Donor ID", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  pdf('FSP_CD4_Proportionplot_DonorsperCluster_res.0.4.pdf')
  ggplot(FSP_CD4@meta.data, aes(fill=Donor, x=RNA_snn_res.0.4)) + theme_classic() +
    geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
    labs(fill = "Donor ID", title='', subtitle='')+
    theme(axis.title.y = element_text(vjust=2.5),plot.margin = margin(45,6,6,6))+
    scale_fill_manual(values=palette.colors(palette = "R4"))
  dev.off()
  
  
  #### 10. CD4 - DEG ####
  DefaultAssay(object = FSP_CD4) <- "RNA"
  FSP_CD4 <- SetIdent(FSP_CD4, value = "RNA_snn_res.0.3")
  
  ##Find Markers that are specific for each cluster
  Batchedmarkers.mast_i_RNA_data_0.3=FindAllMarkers(FSP_CD4, test.use = "MAST", slot='data',logfc.threshold = 0.1,
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
 