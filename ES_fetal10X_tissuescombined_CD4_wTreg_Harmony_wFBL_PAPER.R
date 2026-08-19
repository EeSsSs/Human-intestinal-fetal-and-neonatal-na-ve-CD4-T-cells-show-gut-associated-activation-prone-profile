#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
library(ggbreak)
library(Seurat)
library(sctransform)
library(ggplot2)
library(Matrix)
library(patchwork)
library(reshape2)
library('dplyr')
library('ggplot2')
library('MAST')
library('clustree')
library('patchwork')
library('AUCell')
library(readxl)
library(writexl)
library(DESeq2)
library(tibble)
library(apeglm)
library(glmGamPoi)
library(reticulate)
library(pals)
library(rcartocolor)
library(Seurat.utils)
library(stringr)
library(scRepertoire)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ampir)
library(purrr)
library(tidyr)
library(cowplot)
library(pheatmap)
library(forcats)

# library(batchelor)
# library(SeuratWrappers)

### 1. Setting working directory and loading data / saving data ----------------------------------------------------------------
setwd("~/")
load('FSI_CD4subset_wTreg.Robj')
load('FSP_CD4subset_wTreg.Robj')
load('FLN_CD4subset_wTreg.Robj')
load('FBL_CD4subset_wTreg.Robj')


#### 3. integrate ---------------------------------------------------------------------------------

FLN_CD4 <- SetIdent(FLN_CD4, value = "orig.ident")
FSP_CD4 <- SetIdent(FSP_CD4, value = "orig.ident")
FSI_CD4 <- SetIdent(FSI_CD4, value = "orig.ident")
FBL_CD4 <- SetIdent(FBL_CD4, value = "orig.ident")

colnames(FSP_CD4@meta.data) 
colnames(FSP_CD4@meta.data) [27] <- 'TotalTcells_clusterID'
FSP_CD4$TotalTcells_clusterID <- paste0('FSP - ', FSP_CD4$TotalTcells_clusterID)
FSP_CD4@meta.data <- FSP_CD4@meta.data[,-c(28:38)]

colnames(FLN_CD4@meta.data)
colnames(FLN_CD4@meta.data) [27] <- 'TotalTcells_clusterID'
FLN_CD4$TotalTcells_clusterID <- paste0('FLN - ', FLN_CD4$TotalTcells_clusterID)
FLN_CD4@meta.data <- FLN_CD4@meta.data[,-c(28:38)]

colnames(FSI_CD4@meta.data)
colnames(FSI_CD4@meta.data) [25] <- 'TotalTcells_clusterID'
FSI_CD4$TotalTcells_clusterID <- paste0('FSI - ', FSI_CD4$TotalTcells_clusterID)
FSI_CD4@meta.data <- FSI_CD4@meta.data[,-c(26:36)]

colnames(FBL_CD4@meta.data)
colnames(FBL_CD4@meta.data) [26] <- 'TotalTcells_clusterID'
FBL_CD4$TotalTcells_clusterID <- paste0('FBL - ', FBL_CD4$TotalTcells_clusterID)

### merge to get the layers
tissues.merged <- merge(FSP_CD4, c(FLN_CD4,FSI_CD4, FBL_CD4), merge.data=T)

tissues.merged <- NormalizeData(tissues.merged)
tissues.merged <- FindVariableFeatures(tissues.merged)
VariableFeatures(tissues.merged) <- VariableFeatures(tissues.merged)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(tissues.merged))]
tissues.merged <- ScaleData(tissues.merged)
tissues.merged <- RunPCA(tissues.merged)


###

tissues.integrated <- IntegrateLayers(
  object = tissues.merged, method = HarmonyIntegration,
  orig.reduction = 'pca', new.reduction = "harmony",
  verbose = FALSE
)

CD4_tissues_wTreg <- tissues.integrated


#### 4. exclude d010 and Cluster ####
####
DefaultAssay(CD4_tissues_wTreg) <- 'RNA'

CD4_tissues_wTreg@meta.data[is.na(CD4_tissues_wTreg$Donor),'Donor'] <- 'NA'
CD4_tissues_wTreg_2 <- subset(CD4_tissues_wTreg, subset=Donor!='d010')
CD4_tissues_wTreg_2@meta.data[CD4_tissues_wTreg_2$Donor=='NA','Donor'] <- NA
## recluster
CD4_tissues_wTreg_2 <- FindVariableFeatures(CD4_tissues_wTreg_2)
VariableFeatures(CD4_tissues_wTreg_2) <- VariableFeatures(CD4_tissues_wTreg_2)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(CD4_tissues_wTreg_2))]
CD4_tissues_wTreg_2 <- ScaleData(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2 <- RunPCA(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2 <- RunUMAP(CD4_tissues_wTreg_2, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
CD4_tissues_wTreg_2 <- FindNeighbors(object = CD4_tissues_wTreg_2, reduction = "harmony", dims = 1:30)
CD4_tissues_wTreg_2 <- FindClusters(CD4_tissues_wTreg_2, resolution = seq(0.1:1, by=0.1),
                            cluster.name = c('harmony_clusters_res.0.1', 'harmony_clusters_res.0.2',
                                             'harmony_clusters_res.0.3','harmony_clusters_res.0.4','harmony_clusters_res.0.5',
                                             'harmony_clusters_res.0.6','harmony_clusters_res.0.7','harmony_clusters_res.0.8',
                                             'harmony_clusters_res.0.9','harmony_clusters_res.1.0'))

DimPlot(CD4_tissues_wTreg_2,reduction = "umap.harmony" ,
        group.by = 'harmony_clusters_res.0.7')

#### merge all classical naive clusters --------
DimPlot(CD4_tissues_wTreg_2,reduction = "umap.harmony" ,
        group.by = 'harmony_clusters_res.0.7', pt.size=1)+
  scale_color_manual(values=c(palette.colors(palette = "R4"),
                              'darkgreen','maroon2', 'orange2', 'violet', 'gold','skyblue'))

mapping <- c('0' = '1', '1' = '0', '2' = '0', '3' = '0', '4' = '4',
             '5' = '2', '6' = '6', '7' = '3', '8' = '5', '9' = '7', '10' = '8')
old_clusters <- as.character(CD4_tissues_wTreg_2$harmony_clusters_res.0.7)
new_clusters <- mapping[old_clusters]
names(new_clusters) <- colnames(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2$tissues_res.0.7_manual <- factor(new_clusters)

DimPlot(CD4_tissues_wTreg_2,reduction = "umap.harmony" ,
        group.by = 'tissues_res.0.7_manual', pt.size=1)+
  scale_color_manual(values=c(palette.colors(palette = "R4"),
                              'darkgreen','maroon2', 'orange2', 'violet', 'gold','skyblue'))


mapping_tissue <- c('FBL'='PB','FLN' = 'MLN', 'FSI' = 'SI', 'FSP' = 'SP')
old_tissues <- as.character(CD4_tissues_wTreg_2$orig.ident)
new_tissues <- mapping_tissue[old_tissues]
names(new_tissues) <- colnames(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2$orig.ident <- factor(new_tissues)

#### 5. Figures Paper -------------------------------------------------------------
setwd("")
cluster_colors <- c('black','coral4', 'turquoise','deepskyblue2','orange','maroon','lightblue3','gold', 'seagreen')
tissue_colors <- c('turquoise', 'gold', 'maroon','black')
donor_colors <- c(palette.colors(palette = "R4")[1:4],palette.colors(palette = "R4")[6:7])

DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)

CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_clusters <- as.factor(CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_manual)
levels(CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_clusters) <- c( '0: Classical Naive', '1: Tissue Naive-like',
                                                                     '2: Treg 1','3: Treg 2',
                                                                     '4: Th1','5: Circulating Memory',
                                                                     '6: Early activated/stressed',
                                                                     '7: Th17', '8: Activated/regulatory')
CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_clusters <- factor(CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_clusters,
                                                                 levels=c( '0: Classical Naive', '1: Tissue Naive-like',
                                                                           '6: Early activated/stressed',
                                                                           '5: Circulating Memory',
                                                                           '4: Th1', '7: Th17', 
                                                                           '8: Activated/regulatory', '2: Treg 1', '3: Treg 2'))


#### Fig S3 ####
CD4_tissues_wTreg_2@meta.data$tissues_res.0.8_clusters <- as.factor(CD4_tissues_wTreg_2@meta.data$harmony_clusters_res.0.8)
levels(CD4_tissues_wTreg_2@meta.data$tissues_res.0.8_clusters) <- c( '0: Classical Naive - 1', '1:  Classical Naive - 2', 
                                                                     '2: Tissue Naive-like','3: Classical Naive - 3',
                                                                     '4: Early activated/stressed','5: Treg 1',
                                                                     '6: Circulating Memory',
                                                                     '7: Th1','8: Treg 2',
                                                                     '9: Th17',  '10: Classical Naive - 4', '11: Activated/regulatory')
CD4_tissues_wTreg_2@meta.data$tissues_res.0.8_clusters <- factor(CD4_tissues_wTreg_2@meta.data$tissues_res.0.8_clusters,
                                                                 levels=c( '0: Classical Naive - 1','1:  Classical Naive - 2', 
                                                                           '3: Classical Naive - 3', '10: Classical Naive - 4', 
                                                                           '2: Tissue Naive-like',
                                                                           '4: Early activated/stressed',
                                                                           '6: Circulating Memory',
                                                                           '7: Th1', '9: Th17', 
                                                                           '11: Activated/regulatory', '5: Treg 1', '8: Treg 2'))

figS3b <- DimPlot(object = CD4_tissues_wTreg_2, reduction = "umap.harmony", pt.size=1, group.by = 'orig.ident', 
                 cols=tissue_colors, shuffle=T) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,40),
        text=element_text(size=8))+
  labs(title=expression(bold('Tissue origin of fetal CD4'^+''*' T cells')), color='')+
  guides(color=guide_legend(override.aes=list(size=4)))+
  scale_color_manual(values = tissue_colors, labels=c('MLN', 'Small intestine', 'Spleen', 'Blood'))

figS3c <- DimPlot(object = CD4_tissues_wTreg_2, reduction = "umap.harmony", pt.size=1, group.by = 'harmony_clusters_res.0.8', 
                 label = T, label.size = 12,label.color = c('grey60','grey60','grey60','grey30','grey30',
                                                            'grey30','grey60','grey30','grey30','grey30','grey30','grey30'),
                 repel=T,label.box = T,cols=cluster_colors_3) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))+
  labs(title=expression(bold('Cluster ID of fetal CD4'^+''*' T cells')), color='')+
  guides(color=guide_legend(override.aes=list(size=7)))+
  scale_color_manual(values=cluster_colors_3, 
                     labels= c( '0: Classical Naive - 1', '1:  Classical Naive - 2', 
                                '2: Tissue Naive-like','3 - Classical Naive - 3',
                                '4: Early activated/stressed','5: Treg 1',
                                '6: Circulating Memory',
                                '7: Th1','8: Treg 2',
                                '9: Th17',  '10:Classical Naive - 4', '11: Activated/regulatory'))

percentages_0.5 <- CD4_tissues_wTreg_2@meta.data %>% group_by(orig.ident, tissues_res.0.8_clusters) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup()

percentages_0.5$text_color <- c('grey60','grey60','grey60','grey30','grey30','grey30','grey30',
                                'grey30','grey30','grey30','grey30','grey30',
                                'grey60','grey60','grey60','grey60','grey30','grey30','grey30',
                                'grey30','grey30','grey30',
                                'grey60','grey60','grey30','grey30','grey30','grey30','grey30','grey30','grey30',
                                'grey30','grey30',
                                'grey60','grey60','grey60','grey60','grey30','grey30','grey30',
                                'grey30','grey30','grey30','grey30')

CD4_tissues_wTreg_2$orig.ident <- factor(CD4_tissues_wTreg_2$orig.ident, levels=c('MLN', 'SI', 'SP', 'PB'))
figS3d <- ggplot(CD4_tissues_wTreg_2@meta.data, aes(x=orig.ident, fill=tissues_res.0.8_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  geom_text_repel( data = percentages_0.5, force=0.01,
                   min.segment.length = unit(0.05, 'lines'),
                   aes(x = orig.ident, y = freq,  
                       label = scales::percent(freq, accuracy = 0.1),color=text_color,
                       group = tissues_res.0.8_clusters), 
                   position = position_fill(vjust = 0.5), 
                   size = 5) +
  labs(fill = "Cluster ID", title=expression(bold('Tissue distribution of fetal CD4'^+''*' T cell clusters')), subtitle='')+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,50),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=cluster_colors_3[c(1,2,4,11,3,5,7,8,10,12,6,9)], labels=c('0','1','3','10','2','4','6','7','9','11','5','8'))+
  scale_color_manual(values=c('black', 'grey60'))+
  guides(color='none')

figS3 <- figS3b + figS3c + figS3d + plot_layout(ncol = 2, nrow=2)
ggsave("CD4_tissues_wTreg_noD010_clusteringoverview_figS3_2_res.0.8.pdf", figS3, width = 22, height = 15)


#### 6. DEG ####
#### ROC/MAST for cluster defining markers 
DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2 <- SetIdent(CD4_tissues_wTreg_2, value = "harmony_clusters_res.0.8")
CD4_tissues_wTreg_2 <- SetIdent(CD4_tissues_wTreg_2, value = "tissues_res.0.7_manual")

##Find Markers that are specific for each cluster
Batchedmarkers.mast_i_RNA_data_0.8=FindAllMarkers(CD4_tissues_wTreg_2, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_0.8$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.8$pct.1/Batchedmarkers.mast_i_RNA_data_0.8$pct.2

## Create list
listDEgenes_i_RNA_MAST_0.8<- split(Batchedmarkers.mast_i_RNA_data_0.8, f=Batchedmarkers.mast_i_RNA_data_0.8$cluster)

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
listDEgenes_i_RNA_MAST_0.8 <-lapply(listDEgenes_i_RNA_MAST_0.8, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
listDEgenes_i_RNA_MAST_0.8<-lapply(listDEgenes_i_RNA_MAST_0.8,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

