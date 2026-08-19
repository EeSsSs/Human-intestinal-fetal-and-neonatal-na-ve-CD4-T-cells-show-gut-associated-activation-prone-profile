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
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(writexl)
library(openxlsx)
library(xlsx)
library(readxl)
library(ggstar)

#### 0. load files ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/")
load('EA023SCT_pctmtribo_nocellcylce.Robj')
load('EA025SCT_pctmtribo_nocellcylce.Robj')
load('EA028SCT_pctmtribo_nocellcylce.Robj')

load('Naive.gutblood.merged_SCTmtrb.Robj')

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/")
load('EA024SCT_pctmtribo_nocellcylce.Robj')
load('EA026SCT_pctmtribo_nocellcylce.Robj')
load('EA029SCT_pctmtribo_nocellcylce.Robj')
load('EA030SCT_pctmtribo_nocellcylce.Robj')


## only for part 1
list_Naive <- list(EA023_mtrb, EA024_mtrb, EA025_mtrb, EA026_mtrb, EA028_mtrb, EA029_mtrb, EA030_mtrb)


#### 2.2 combine with blood ####
## merge without batch-correction
Naive.gutblood.merged <- merge(EA029_mtrb, c(EA024_mtrb, EA026_mtrb, EA023_mtrb, EA025_mtrb, EA028_mtrb, EA030_mtrb), merge.data=T)


setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
save(Naive.gutblood.merged, file='Naive.gutblood.merged_SCTmtrb_children.Robj')
load(file = "T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Naive.gutblood.merged_SCTmtrb_children.Robj")

## after merge, you have to set features again: variablefeatures are just all SCT features
Naive.gutblood.merged <- RunPCA(Naive.gutblood.merged, features = rownames(Naive.gutblood.merged@assays[["SCT"]]@scale.data))
Naive.gutblood.merged <- RunUMAP(Naive.gutblood.merged, dims = 1:30)

##clustering
Naive.gutblood.merged <- FindNeighbors(object = Naive.gutblood.merged, reduction = "pca", dims = 1:30)
Naive.gutblood.merged <- FindClusters(object = Naive.gutblood.merged, resolution = seq(0,1.3,0.1))

##add group
Naive.gutblood.merged$group.ident <- Naive.gutblood.merged$orig.ident
Naive.gutblood.merged@meta.data[Naive.gutblood.merged$orig.ident=='EA024'|Naive.gutblood.merged$orig.ident=='EA026'|
                                  Naive.gutblood.merged$orig.ident=='EA029'|Naive.gutblood.merged$orig.ident=='EA030',
                       'group.ident'] <- 'Blood'
Naive.gutblood.merged@meta.data[Naive.gutblood.merged$orig.ident=='EA023'|Naive.gutblood.merged$orig.ident=='EA025'|
                                  Naive.gutblood.merged$orig.ident=='EA028',
                                'group.ident'] <- 'Gut'
Naive.gutblood.merged$donor.ident <- Naive.gutblood.merged$orig.ident
Naive.gutblood.merged@meta.data[Naive.gutblood.merged$orig.ident=='EA023'|Naive.gutblood.merged$orig.ident=='EA024',
                                'donor.ident'] <- 'HINT139'
Naive.gutblood.merged@meta.data[Naive.gutblood.merged$orig.ident=='EA025'|Naive.gutblood.merged$orig.ident=='EA026',
                                'donor.ident'] <- 'HINT140'
Naive.gutblood.merged@meta.data[Naive.gutblood.merged$orig.ident=='EA028'|Naive.gutblood.merged$orig.ident=='EA029',
                                'donor.ident'] <- 'HINT148'
Naive.gutblood.merged@meta.data[Naive.gutblood.merged$orig.ident=='EA030',
                                'donor.ident'] <- 'HINT149'
#### visualization - adult ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
DefaultAssay(Naive.gutblood.merged) <- 'RNA'

Fig_clustree.merged <- clustree(Naive.gutblood.merged)
ggsave(plot=Fig_clustree.merged, file='GutBlood_Merged_Clustree_SCTonmtrb_children.pdf', width=8, height=12)

Fig_Dimplot.merged_OrigID <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", group.by = 'orig.ident', pt.size=0.9)+
  scale_color_manual(values=c('pink','green', 'violet',  'seagreen3','purple4','seagreen1','lightblue'))+
                                      ggtitle('Sample ID', sub= 'Children Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_groupID <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", group.by = 'group.ident', pt.size=0.9)+
  scale_color_manual(values=c('blue3',  'orange'))+
                                      ggtitle('Gut vs Blood', sub= 'Children Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_DonorID <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", group.by = 'donor.ident', pt.size=0.9)+
  scale_color_manual(values=c( 'seagreen3','purple','skyblue','gold'))+
                                      ggtitle('Donor ID', sub= 'Children Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Naive.gutblood.merged@active.ident <- Naive.gutblood.merged$SCT_snn_res.0.6
Fig_Dimplot.merged_res0.6 <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", pt.size=0.9,
                                     cols=c('black','skyblue2'))+
       ggtitle('Clustering resolution 0.6', sub= 'Children Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Naive.gutblood.merged@active.ident <- Naive.gutblood.merged$SCT_snn_res.0.8
Fig_Dimplot.merged_res0.8 <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", pt.size=0.9,
                                     cols=c('black','skyblue2','orange'))+
  ggtitle('Clustering resolution 0.8', sub= 'Children Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Naive.gutblood.merged@active.ident <- Naive.gutblood.merged$SCT_snn_res.0.9
Fig_Dimplot.merged_res0.9 <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", pt.size=0.9)+
                                ggtitle('Clustering resolution 0.9', sub= 'SCT on pct.mtrb')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5))
                              
Fig_Dimplot.merged_Phase <- DimPlot(object = Naive.gutblood.merged, reduction = "umap", group.by = 'Phase', pt.size=0.9)+
  scale_color_colorblind()
Fig_Featureplot.merged_nFeature <- FeaturePlot(Naive.gutblood.merged, features='nFeature_RNA', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_mtrb <- FeaturePlot(Naive.gutblood.merged, features='percent.mt', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_rb <- FeaturePlot(Naive.gutblood.merged, features='percent.ribo', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Naive.gutblood.merged$flowCD31_bin <- Naive.gutblood.merged$flowCD31 > 300
Fig_Featureplot.merged_CD31 <- DimPlot(Naive.gutblood.merged, group.by ='flowCD31_bin', pt.size=0.9, order=T) + 
  scale_color_manual(values=c('orange2', 'skyblue')) + ggtitle('surface CD31+')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_clustering_ext <- Fig_Dimplot.merged_OrigID + Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.8+Fig_Dimplot.merged_res0.9+Fig_Featureplot.merged_CD31+
  Fig_Featureplot.merged_mtrb+ Fig_Featureplot.merged_rb + Fig_Featureplot.merged_nFeature+Fig_Dimplot.merged_Phase+
  plot_layout(ncol=5)

Fig_clustering <- Fig_Dimplot.merged_DonorID + Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.8+
  plot_layout(ncol=4)

ggsave('GutBlood_Merged_ClusteringOverview_SCTonmtrb_children.pdf',Fig_clustering_ext, width=24, height=10)
ggsave('GutBlood_Merged_ClusteringOverview_simple_SCTonmtrb_children.pdf',Fig_clustering, width=24, height=6)

Fig_Vlnplot.merged_CD31 <- VlnPlot(Naive.gutblood.merged, group.by = 'SCT_snn_res.0.6', feature='flowCD31')+
  VlnPlot(Naive.gutblood.merged, group.by = 'SCT_snn_res.0.8', feature='flowCD31')+
  VlnPlot(Naive.gutblood.merged, group.by = 'SCT_snn_res.0.9', feature='flowCD31')
ggsave('GutBlood_Merged_CD31Vlns_SCTonmtrb_children.pdf',Fig_Vlnplot.merged_CD31, width=14, height=20)

Fig_barplot_res0.6 <- ggplot(Naive.gutblood.merged@meta.data, aes(x=group.ident, fill=SCT_snn_res.0.6)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fiEA = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.5 - SCT on pct.mtrb')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))

Fig_barplot_res0.8 <- ggplot(Naive.gutblood.merged@meta.data, aes(x=group.ident, fill=SCT_snn_res.0.8)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fiEA = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.8 - SCT on pct.mtrb  ')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))

Fig_barplot_res0.9 <- ggplot(Naive.gutblood.merged@meta.data, aes(x=group.ident, fill=SCT_snn_res.0.9)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fiEA = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.9 - SCT on pct.mtrb  - noEA010 ')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))

Fig_barplots <- Fig_barplot_res0.6 + Fig_barplot_res0.8 + Fig_barplot_res0.9 + plot_layout(ncol=3)
ggsave(plot=Fig_barplots, filename='GutBlood_Merged_Barplot_groupsperclusters_SCTonmtrb_children.pdf',width=16, height=5)


#### visualize potential gut-enriched naive T cell genes
DefaultAssay(object = Naive.gutblood.merged) <- "RNA"
Naive.gutblood.merged@active.ident <- Naive.gutblood.merged$SCT_snn_res.0.8
candidate_genes <- c('CXCR4','TNFAIP3','CTLA4','SOCS3', 
                     'KLF6', 'IL2RA', 'CCL5', 'ICOS',
                     'BTG1', 'LEF1','KLF2','SOX4',
                     'CREM','SRGN', 'SLC2A3',  'PTGER4',
                     'GPR183', 'CD69', 'LDLR', 'HIF1A')
FeaturePlot(Naive.gutblood.merged, features=candidate_genes)
VlnPlot(Naive.gutblood.merged, features=candidate_genes)

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
Vln_GutCandidates <- VlnPlot(Naive.gutblood.merged, features=candidate_genes, pt.size=0.00001) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
Feature_GutCandidates <- FeaturePlot(Naive.gutblood.merged, features=candidate_genes, pt.size=1) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
ggsave(plot=Vln_GutCandidates, filename='CandidateGenes_Vlnplot_Gutenriched_NaiveTcells_children.pdf', height = 16, width=18)
ggsave(plot=Feature_GutCandidates, filename='CandidateGenes_FeaturePlot_Gutenriched_NaiveTcells_children.pdf', height = 16, width=18)


#### heatmap
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
load('ListDE_Naive.gutblood.merged.SCTonmtrb_res0.8_allmarkers_MASTRNAdata_children.Robj')
top50 <- lapply(ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata, function(df) {df %>% subset(avg_log2FC>0) %>% top_n(50, avg_log2FC)})
top50 <- bind_rows(top50)

Naive.gutblood.merged <- ScaleData(Naive.gutblood.merged)
c <- DoHeatmap(
  Naive.gutblood.merged,
  features = top50$gene,
  cells = NULL,
  group.by = "SCT_snn_res.0.8",
  group.bar = TRUE,
  disp.min = -2.5,
  disp.max = NULL,
  label = TRUE,
  size = 6,
  hjust = 0,
  angle = 0,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.01,
  combine = TRUE
) + 
  theme(text=element_text(size=18), legend.text = element_text(size=18), 
        legend.title = element_text(size=18), 
        plot.margin = margin(0,0,10,0), legend.key.size = unit(.5, 'cm'))
c
#### Differential gene expression - adult ####
DefaultAssay(object = Naive.gutblood.merged) <- "RNA"
Naive.gutblood.merged <- JoinLayers(Naive.gutblood.merged)

Naive.gutblood.merged@active.ident <- Naive.gutblood.merged$SCT_snn_res.0.8

##Find Markers that are specific for each cluster
Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata=FindAllMarkers(Naive.gutblood.merged, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)

## Create list
ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata<- split(Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata, 
                                                          f=Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata$cluster)
## Filter on adj.P-value
ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata <-lapply(ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata<-lapply(ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

##save as Robj
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
save(ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata, file='ListDE_Naive.gutblood.merged.SCTonmtrb_res0.8_allmarkers_MASTRNAdata_children.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(ListDE_Naive.gutblood.merged_res0.8_allmarkers_MASTRNAdata, file='ListDE_Naive.gutblood.merged.SCTonmtrb_res0.8_allmarkers_MASTRNAdata_children.xlsx')


#### 2.4 PAPER combine children with adult ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/")
load('EA012SCT_pctmtribo_nocellcylce.Robj')
load('EA014SCT_pctmtribo_nocellcylce.Robj')
load('EA016SCT_pctmtribo_nocellcylce.Robj')

load('EA021SCT_pctmtribo_nocellcylce.Robj')

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/")
load('EA013SCT_pctmtribo_nocellcylce.Robj')
load('EA015SCT_pctmtribo_nocellcylce.Robj')
load('EA017SCT_pctmtribo_nocellcylce.Robj')
load('EA022SCT_pctmtribo_nocellcylce.Robj')

### 
Naive.gutblood.merged.all <- merge(EA012_mtrb, c(EA013_mtrb, EA014_mtrb, EA015_mtrb, EA016_mtrb, EA017_mtrb, EA021_mtrb, EA022_mtrb,
                                            EA023_mtrb, EA024_mtrb, EA025_mtrb, EA026_mtrb, EA028_mtrb, EA029_mtrb, EA030_mtrb))

#### load merged children/adult object ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
save(Naive.gutblood.merged.all, file='Naive.gutblood.merged.all_SCTmtrb_childrenadults.Robj')
load(file = "T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Naive.gutblood.merged.all_SCTmtrb_childrenadults.Robj")

#### prepare children/adult merged object ####
## after merge, you have to set features again: variablefeatures are just all SCT features
Naive.gutblood.merged.all <- RunPCA(Naive.gutblood.merged.all, features = rownames(Naive.gutblood.merged.all@assays[["SCT"]]@scale.data))
Naive.gutblood.merged.all <- RunUMAP(Naive.gutblood.merged.all, dims = 1:30)

##clustering
Naive.gutblood.merged.all <- FindNeighbors(object = Naive.gutblood.merged.all, reduction = "pca", dims = 1:30)
Naive.gutblood.merged.all <- FindClusters(object = Naive.gutblood.merged.all, resolution = seq(0,1.3,0.1))

##add group
Naive.gutblood.merged.all$group.ident <- Naive.gutblood.merged.all$orig.ident
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$orig.ident=='EA013'|Naive.gutblood.merged.all$ orig.ident=='EA015'|
                                      Naive.gutblood.merged.all$ orig.ident=='EA017'|Naive.gutblood.merged.all$ orig.ident=='EA022'|
                                      Naive.gutblood.merged.all$orig.ident=='EA024'|Naive.gutblood.merged.all$orig.ident=='EA026'|
                                  Naive.gutblood.merged.all$orig.ident=='EA029'|Naive.gutblood.merged.all$orig.ident=='EA030',
                                'group.ident'] <- 'Blood'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$ orig.ident=='EA012'|Naive.gutblood.merged.all$ orig.ident=='EA014'|
                                      Naive.gutblood.merged.all$ orig.ident=='EA016'|Naive.gutblood.merged.all$ orig.ident=='EA021'|
                                      Naive.gutblood.merged.all$orig.ident=='EA023'|Naive.gutblood.merged.all$orig.ident=='EA025'|
                                  Naive.gutblood.merged.all$orig.ident=='EA028',
                                'group.ident'] <- 'Colon'
Naive.gutblood.merged.all$donor.ident <- Naive.gutblood.merged.all$orig.ident
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$orig.ident=='EA023'|Naive.gutblood.merged.all$orig.ident=='EA024',
                                'donor.ident'] <- 'HINT139'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$orig.ident=='EA025'|Naive.gutblood.merged.all$orig.ident=='EA026',
                                'donor.ident'] <- 'HINT140'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$orig.ident=='EA028'|Naive.gutblood.merged.all$orig.ident=='EA029',
                                'donor.ident'] <- 'HINT148'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$orig.ident=='EA030',
                                'donor.ident'] <- 'HINT149'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$ orig.ident=='EA012'|Naive.gutblood.merged.all$ orig.ident=='EA013',
                                'donor.ident'] <- 'HINT129'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$ orig.ident=='EA014'|Naive.gutblood.merged.all$ orig.ident=='EA015',
                                'donor.ident'] <- 'HINT130'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$ orig.ident=='EA016'|Naive.gutblood.merged.all$ orig.ident=='EA017',
                                'donor.ident'] <- 'HINT131'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$ orig.ident=='EA021'|Naive.gutblood.merged.all$ orig.ident=='EA022',
                                'donor.ident'] <- 'HINT136'
Naive.gutblood.merged.all$age.ident <- Naive.gutblood.merged.all$orig.ident
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$orig.ident=='EA012'|Naive.gutblood.merged.all$ orig.ident=='EA013'|
                                      Naive.gutblood.merged.all$ orig.ident=='EA014'|Naive.gutblood.merged.all$ orig.ident=='EA015'|
                                      Naive.gutblood.merged.all$orig.ident=='EA016'|Naive.gutblood.merged.all$orig.ident=='EA017'|
                                      Naive.gutblood.merged.all$orig.ident=='EA021'|Naive.gutblood.merged.all$orig.ident=='EA022',
                                    'age.ident'] <- 'Children'
Naive.gutblood.merged.all@meta.data[Naive.gutblood.merged.all$ orig.ident=='EA023'|Naive.gutblood.merged.all$ orig.ident=='EA024'|
                                      Naive.gutblood.merged.all$ orig.ident=='EA025'|Naive.gutblood.merged.all$ orig.ident=='EA026'|
                                      Naive.gutblood.merged.all$orig.ident=='EA028'|Naive.gutblood.merged.all$orig.ident=='EA029'|
                                      Naive.gutblood.merged.all$orig.ident=='EA030',
                                    'age.ident'] <- 'Adults'
#### Surface data Sorting (see also T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/FACS verification/Graphs/scRNAseq_GutProfile_FACSverification_P1.R) ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut sorting")
SortFACSdata_HINTdonors_final <- read_excel("SortFACSdata_HINTdonors_final.xlsx")
colnames(SortFACSdata_HINTdonors_final)[1] <- 'X'
SortFACSdata_HINTdonors_final$Sequenced <- as.logical(SortFACSdata_HINTdonors_final$Sequenced )
SortFACSdata_HINTdonors_final[SortFACSdata_HINTdonors_final$Donor=='HINT139'&SortFACSdata_HINTdonors_final$Tissue=='Gut',
                              7:12] <- NA
SortFACSdata_HINTdonors_final[SortFACSdata_HINTdonors_final$Donor=='HINT139'&SortFACSdata_HINTdonors_final$Tissue=='Gut',
                              13] <- 11.2
SortFACSdata_HINTdonors_final[SortFACSdata_HINTdonors_final$Donor=='HINT139'&SortFACSdata_HINTdonors_final$Tissue=='Gut',
                              14] <- 10.7
SortFACSdata_HINTdonors_final[SortFACSdata_HINTdonors_final$Donor=='HINT139'&SortFACSdata_HINTdonors_final$Tissue=='Gut',
                              17] <- 96.2
SortFACSdata_HINTdonors_final[SortFACSdata_HINTdonors_final$Donor=='HINT139'&SortFACSdata_HINTdonors_final$Tissue=='Gut',
                              16] <- 10.7*0.962
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/FACS verification")
load('FACSverification_simpleformat.RData')
HINT$Sequenced <- rep(F, length(HINT$X))

Naive.gutblood.FACS <- rbind(SortFACSdata_HINTdonors_final[,c(1:6,16)], HINT[,c(1,10,120:124)])
Naive.gutblood.FACS$Age <- round(as.numeric(Naive.gutblood.FACS$Age))

#### Paper Figures ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/Paper - T cell development & the gut")

## calculate %Naive of children/adult in the tissues
median(unlist(Naive.gutblood.FACS[Naive.gutblood.FACS$Donor!='HINT135'&Naive.gutblood.FACS$Age>18&
                                    Naive.gutblood.FACS$Tissue=='Blood','TrueNaive_ofCD4']))
median(unlist(Naive.gutblood.FACS[Naive.gutblood.FACS$Donor!='HINT135'&Naive.gutblood.FACS$Age>18&
                                    Naive.gutblood.FACS$Tissue=='Gut','TrueNaive_ofCD4']))
median(unlist(Naive.gutblood.FACS[Naive.gutblood.FACS$Donor!='HINT135'&Naive.gutblood.FACS$Age<18&
                                    Naive.gutblood.FACS$Tissue=='Blood','TrueNaive_ofCD4']))
median(unlist(Naive.gutblood.FACS[Naive.gutblood.FACS$Donor!='HINT135'&Naive.gutblood.FACS$Age<18&
                                    Naive.gutblood.FACS$Tissue=='Gut','TrueNaive_ofCD4']))


fig4a <- ggplot(Naive.gutblood.FACS[Naive.gutblood.FACS$Donor!='HINT135',], aes(Age,TrueNaive_ofCD4))+
  geom_star(aes(fill=Tissue,alpha=Sequenced, starshape=Donor), size=7, starstroke=1)+
  scale_starshape_manual(values=c(1, 13, 15, 11, 12, 14, 29, 2, 27, 3,4,5,9,22,25))+
  scale_fill_manual(values=c('blue','orange'),labels=c('Blood','Colon'))+
  scale_alpha_manual(values=c(0.5,1), labels=c('No','Yes'))+
  scale_y_sqrt(limits = c(NA,100),expand = c(0,0),breaks=c(0,0.5,1,5,10,50,100),
               labels=c('0%','0.5%','1%','5%','10%','50%','100%'))+
  scale_x_continuous(limits=c(0,50),expand=c(0,0))+
  ggtitle(expression(bold('Naive CD4'^+''*' T cells')),subtitle='Tissue distribution over age')+
  ylab(expression('% of CD4'^+''*' T Cells'))+
  xlab('Age (years)')+
  theme_bw()+
  theme(axis.line = element_line(linewidth = .75), panel.border = element_blank(),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=28),axis.title.y = element_text(size=32,vjust=1.5),
        axis.title.x = element_text(size=32),
        legend.text = element_text(size=32),legend.title = element_text(size=34),
        plot.title = element_text(size=38, face='bold', hjust=0.5),
        plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(20,80,20,20),panel.spacing = unit(2, "lines"),
        text=element_text(size=8), strip.text = element_text(size=32))+
  guides(starshape='none', fill=guide_legend(override.aes = list(size = 12)), alpha=guide_legend(override.aes = list(size = 12)))
fig4a
## calculate % cluster 2 in gut vs blood
sum(Naive.gutblood.merged.all$SCT_snn_res.0.6=='2'&Naive.gutblood.merged.all$group.ident=='Colon')/sum(Naive.gutblood.merged.all$group.ident=='Colon')*100
sum(Naive.gutblood.merged.all$SCT_snn_res.0.6=='2'&Naive.gutblood.merged.all$group.ident=='Blood')/sum(Naive.gutblood.merged.all$group.ident=='Blood')*100

fig4b <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", group.by = 'group.ident', 
                 pt.size=2)+
  scale_color_manual(values=c('blue3',  'orange'))+
  ggtitle(expression(bold('Tissue origin')), sub= expression('Children & Adult Naive CD4'^+''*' T cells'))+
  theme(axis.line = element_line(linewidth = .75), panel.border = element_blank(),
        axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.text = element_text(size=32),legend.title = element_text(size=34),
        plot.title = element_text(size=38, face='bold', hjust=0.5),
        plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(40,40,40,80)) + guides(color=guide_legend(override.aes=list(size=12)))

fig4c <- DimPlot(object = Naive.gutblood.merged.all, group.by = 'SCT_snn_res.0.6', reduction = "umap",
                 pt.size=2,
        cols=c('gold','darkgreen','magenta'))+
  ggtitle(expression(bold('Cluster ID')), sub= expression('Children & Adult Naive CD4'^+''*' T cells'))+
  theme(axis.line = element_line(linewidth = .75), panel.border = element_blank(),
        axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.text = element_text(size=32),legend.title = element_text(size=34),
        plot.title = element_text(size=38, face='bold', hjust=0.5),
        plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(40,40,40,80)) + guides(color=guide_legend(override.aes=list(size=12)))

## calculate % Cluster 2 in children vs adult gut
sum(Naive.gutblood.merged.all$SCT_snn_res.0.6=='2'&Naive.gutblood.merged.all$age.ident=='Children'&
      Naive.gutblood.merged.all$group.ident=='Colon')/sum(Naive.gutblood.merged.all$age.ident=='Children'&
                                                            Naive.gutblood.merged.all$group.ident=='Colon')*100
sum(Naive.gutblood.merged.all$SCT_snn_res.0.6=='2'&Naive.gutblood.merged.all$age.ident=='Adults'&
      Naive.gutblood.merged.all$group.ident=='Colon')/sum(Naive.gutblood.merged.all$age.ident=='Adults'&
                                                            Naive.gutblood.merged.all$group.ident=='Colon')*100

fig4d <- ggplot(Naive.gutblood.merged.all@meta.data %>% filter(group.ident=='Colon'), aes(x=SCT_snn_res.0.6, fill=age.ident)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab(expression("Fraction of cluster")) + 
  scale_fill_manual(values=c('black','pink'))+
  labs(fill = "Age group", title='Children vs Adult - cluster distribution', subtitle=expression('Colon Naive CD4'^+''*' T cells'))+
  theme(axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=28),axis.title.y = element_text(size=32,vjust=1.5),
        axis.title.x = element_text(size=32),
        legend.text = element_text(size=32),legend.title = element_blank(),
        plot.title = element_text(size=38, face='bold', hjust=0.5),plot.subtitle = element_text(size=36, hjust=0.5, vjust=2),
        plot.margin = margin(20,20,20,40),
        text=element_text(size=8))+ guides(fill=guide_legend(override.aes=list(size=12)))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())

fig4 <- fig4a + fig4b + fig4c + fig4d + plot_layout(ncol=4)
ggsave("GutBlood_Merged_res.0.6_childrenadults_clusteringoverview_fig4.pdf", fig4, width = 40, height = 10)

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
C2genes <- read_excel("SupplData10_Naive.gutblood.merged.all.SCTonmtrb_res0.6_manual_markers_MASTRNAdata_childrenadults_C2vcC10.xlsx")
Fig_Volcano_C2vsC01 <- ggplot(C2genes, 
                             aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(colour = abs(avg_log2FC)), size=6) +
  ggtitle(expression(bold('Differential gene expression in Colon-enriched cluster 2')), 
          subtitle = 'Clusters 0&1                                                                                                         Cluster 2') +
  geom_text_repel(aes(label=gene,x = avg_log2FC, y = -log10(p_val_adj)), 
                  size=10, direction='both', nudge_y = 0.25,
                  max.overlaps = 11)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_gradient(low = "gold", high = "blue") +
  scale_y_continuous(limits=c(0,50),expand=c(0,0))+
  scale_x_continuous(limits=c(-2,5.5))+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        axis.text = element_text(size=36),
        axis.title.x = element_text(size=40, vjust=-2),axis.title.y = element_text(size=40, vjust=3),
        legend.text = element_text(size=28),legend.title = element_text(size=30),
        plot.title = element_text(size=50, face='bold', hjust=0.5, vjust=2),
        plot.subtitle = element_text(size=52, hjust=0.5, vjust=1),
        plot.margin = margin(40,40,40,40))
Fig_Volcano_C2vsC01
ggsave(Fig_Volcano_C2vsC01,filename=('C2vsC01_res.0.6_manual_NaiveGutBlood_merged_childrenadults_Volcano_fig4e.pdf'), height=15, width=30)

DESeq2results_Naivegutblood_PCApseudobulk_paired_childrenadults_excl_clust2 <- read_excel("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels/SupplData11_DESeq2results_Naivegutblood_pseudobulk_childrenadults_excl.clust2_gutvsblood_paired.xlsx")
Fig_Volcano_ColonvsBlood <- ggplot(DESeq2results_Naivegutblood_PCApseudobulk_paired_childrenadults_excl_clust2, 
                                   aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(colour = abs(log2FoldChange)), size=6) +
  ggtitle(expression(bold('Pseudobulk differential gene expression between tissues in clusters 0&1')), 
          subtitle = 'Blood                                                                                                                Colon') +
  geom_text_repel(aes(label=Gene,x = log2FoldChange, y = -log10(padj)), 
                  size=10, direction='both', nudge_y = 0.25,
                  max.overlaps = 15)+
  geom_text_repel(data=DESeq2results_Naivegutblood_PCApseudobulk_paired_childrenadults_excl_clust2[DESeq2results_Naivegutblood_PCApseudobulk_paired_childrenadults_excl_clust2$Gene=='ICOS'|
                                                                                                     DESeq2results_Naivegutblood_PCApseudobulk_paired_childrenadults_excl_clust2$Gene=='STAT3',],
                  aes(label=Gene,x = log2FoldChange, y = -log10(padj)), 
                  size=7, direction='both', nudge_y = -0.35,
                  max.overlaps = 15)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_gradient(low = "gold", high = "blue") +
  scale_y_continuous(limits=c(0,30),expand=c(0,0))+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        axis.text = element_text(size=36),
        axis.title.x = element_text(size=40, vjust=-2),axis.title.y = element_text(size=40, vjust=3),
        legend.text = element_text(size=28),legend.title = element_text(size=30),
        plot.title = element_text(size=50, face='bold', hjust=0.5, vjust=2),
        plot.subtitle = element_text(size=52, hjust=0.5, vjust=1),
        plot.margin = margin(40,40,40,40))
Fig_Volcano_ColonvsBlood

ggsave(Fig_Volcano_ColonvsBlood,filename=('GutvsBlood_res.0.6_clusters01_NaiveGutBlood_merged_childrenadults_Volcano_fig4f.pdf'), height=15, width=30)

DefaultAssay(object = Naive.gutblood.merged.all) <- "RNA"
Naive.gutblood.merged.all <- JoinLayers(Naive.gutblood.merged.all)
Naive.gutblood.merged.all@active.ident <- Naive.gutblood.merged.all$SCT_snn_res.0.6
candidate_genes <- c('ICOS','TNFAIP3','CTLA4','SOCS3', 
                     'KLF6', 'CD69', 'CCL5', 'NR4A2',
                     'CXCR4','SLC2A3',  'PTGER4','GPR183',
                     'CD38')

fig4g <- VlnPlot(Naive.gutblood.merged.all, features=candidate_genes,
        split.by = 'group.ident', ncol=13, pt.size=0.000001, assay='RNA', slot='data', alpha=0.25,same.y.lims = T,
        cols=c('blue','orange'))&
  scale_x_discrete(labels=c('0','1','2')) &
  scale_y_continuous(expand=c(0,0))&
  xlab('Cluster ID')&
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18), legend.text = element_text(size=16),
        plot.title = element_text(size=22, hjust=0.5), plot.subtitle = element_text(size=20, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=16))

ggsave(fig4g,filename=('GutvsBlood_res.0.6_clustersVlnPlot_NaiveGutBlood_merged_childrenadults_fig4g.pdf'), height=4, width=30)

#### Figures - children/adult  - rest ####

##visualization
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
DefaultAssay(Naive.gutblood.merged.all) <- 'RNA'

Fig_clustree.merged <- clustree(Naive.gutblood.merged.all, prefix='SCT_snn_res.')
ggsave(plot=Fig_clustree.merged, file='GutBlood_Merged_Clustree_SCTonmtrb_childrenadults.pdf', width=8, height=12)

Fig_Dimplot.merged_OrigID <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", group.by = 'orig.ident', pt.size=0.9)+
  scale_color_manual(values=c('pink','green', 'violet',  'seagreen3','purple4','seagreen1','deeppink','darkgreen',
                                    'pink3','green3', 'magenta',  'seagreen4','purple2','seagreen2','palegreen1'))+
  ggtitle('Sample ID', sub= 'Children&Adults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_groupID <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", group.by = 'group.ident', pt.size=0.9)+
  scale_color_manual(values=c('blue3',  'orange'))+
  ggtitle('Gut vs Blood', sub= 'Children&Adult Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_ageID <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", group.by = 'age.ident', pt.size=0.9)+
  scale_color_manual(values=c('pink',  'black'))+
  ggtitle('Adult vs Children', sub= 'Children&Adult Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_DonorID <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", group.by = 'donor.ident', pt.size=0.9)+
  scale_color_manual(values=c( 'seagreen3','purple3','skyblue','gold','black','orange','cyan','magenta'))+
  ggtitle('Donor ID', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_Dimplot.merged_res0.6 <- DimPlot(object = Naive.gutblood.merged.all, group.by = 'SCT_snn_res.0.6', reduction = "umap", pt.size=2,
                                     cols=c('black','violet','orange'))+
  ggtitle('Clustering resolution 0.6', sub= 'Children&Adults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Naive.gutblood.merged.all@active.ident <- Naive.gutblood.merged.all$SCT_snn_res.0.9
Fig_Dimplot.merged_res0.9 <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", pt.size=0.9,
                                     cols=c('black','violet','cyan','orange'))+
  ggtitle('Clustering resolution 0.9', sub= 'SCT on pct.mtrb')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5))

Fig_Dimplot.merged_Phase <- DimPlot(object = Naive.gutblood.merged.all, reduction = "umap", group.by = 'Phase', pt.size=0.9)+
  scale_color_colorblind()
Fig_Featureplot.merged_nFeature <- FeaturePlot(Naive.gutblood.merged.all, features='nFeature_RNA', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_mtrb <- FeaturePlot(Naive.gutblood.merged.all, features='percent.mt', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_rb <- FeaturePlot(Naive.gutblood.merged.all, features='percent.ribo', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Naive.gutblood.merged.all$flowCD31_bin <- Naive.gutblood.merged.all$flowCD31 > 300
Fig_Featureplot.merged_CD31 <- DimPlot(Naive.gutblood.merged.all, group.by ='flowCD31_bin', pt.size=0.9, order=T) + 
  scale_color_manual(values=c('orange2', 'skyblue')) + ggtitle('surface CD31+')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_clustering_ext <- Fig_Dimplot.merged_OrigID + Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.9+Fig_Featureplot.merged_CD31+
  Fig_Featureplot.merged_mtrb+ Fig_Featureplot.merged_rb + Fig_Featureplot.merged_nFeature+Fig_Dimplot.merged_Phase+
  plot_layout(ncol=5)

Fig_clustering <- Fig_Dimplot.merged_DonorID + Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_ageID+
  plot_layout(ncol=4)

ggsave('GutBlood_Merged_ClusteringOverview_SCTonmtrb_childrenadults.pdf',Fig_clustering_ext, width=24, height=10)
ggsave('GutBlood_Merged_ClusteringOverview_simple_SCTonmtrb_childrenadults.pdf',Fig_clustering, width=24, height=6)

Fig_Vlnplot.merged_CD31 <- VlnPlot(Naive.gutblood.merged.all, group.by = 'SCT_snn_res.0.6', feature='flowCD31')+
  VlnPlot(Naive.gutblood.merged.all, group.by = 'SCT_snn_res.0.9', feature='flowCD31')
ggsave('GutBlood_Merged_CD31Vlns_SCTonmtrb_childrenadults.pdf',Fig_Vlnplot.merged_CD31, width=14, height=20)

Fig_barplot_res0.6 <- ggplot(Naive.gutblood.merged.all@meta.data, aes(x=group.ident, fill=SCT_snn_res.0.6)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  scale_fill_manual(values=c('black','skyblue2','orange'))+
  labs(fill = "Cluster", title='Proportionplot clusters/organ', subtitle='res 0.6 - SCT on pct.mtrb - Children&Adults')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))
  
Fig_barplot_res0.6_organ <- ggplot(Naive.gutblood.merged.all@meta.data, aes(x=SCT_snn_res.0.6, fill=group.ident)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    scale_fill_manual(values=c('blue','orange'))+
    labs(fill = "Organ", title='Proportionplot organ/clusters', subtitle='res 0.6 - SCT on pct.mtrb - Children&Adults')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
          axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
    scale_y_continuous(expand = c(0,0))
  
Fig_barplot_res0.6_kids <- ggplot(Naive.gutblood.merged.all@meta.data, aes(x=SCT_snn_res.0.6, fill=age.ident)) + theme_classic() +
    geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
    scale_fill_manual(values=c('pink','black'))+
  labs(fill = "Age group", title='Proportionplot age group/clusters', subtitle='res 0.6 - SCT on pct.mtrb - Children&Adults')+
    theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
          axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
    scale_y_continuous(expand = c(0,0))

Fig_barplots <- Fig_barplot_res0.6 + Fig_barplot_res0.6_organ +  Fig_barplot_res0.6_kids +plot_layout(ncol=3)
ggsave(plot=Fig_barplots, filename='GutBlood_Merged_Barplot_groupsperclusters_SCTonmtrb_childrenadults.pdf',width=16, height=5)


#### visualize potential gut-enriched naive T cell genes
DefaultAssay(object = Naive.gutblood.merged.all) <- "RNA"
Naive.gutblood.merged.all <- JoinLayers(Naive.gutblood.merged.all)
Naive.gutblood.merged.all@active.ident <- Naive.gutblood.merged.all$SCT_snn_res.0.6
candidate_genes <- c('CXCR4','TNFAIP3','CTLA4','SOCS3', 
                     'KLF6', 'IL2RA', 'CCL5', 'ICOS',
                     'BTG1', 'LEF1','KLF2','SOX4',
                     'CREM','SRGN', 'SLC2A3',  'PTGER4',
                     'GPR183', 'CD69', 'LDLR', 'HIF1A')
FeaturePlot(Naive.gutblood.merged.all, features=candidate_genes)
##joinlayers
VlnPlot(Naive.gutblood.merged.all, features=candidate_genes)

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
Vln_GutCandidates <- VlnPlot(Naive.gutblood.merged.all, features=candidate_genes, pt.size=0.00001) & 
  scale_fill_manual(values=c('black','skyblue2','orange'))&
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
Feature_GutCandidates <- FeaturePlot(Naive.gutblood.merged.all, features=candidate_genes, pt.size=1) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))

Vln_GutCandidates_organ <- VlnPlot(Naive.gutblood.merged.all, features=candidate_genes,pt.size=0.00001, group.by = 'group.ident') & 
  scale_fill_manual(values=c('blue','orange'))&
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))

Vln_GutCandidates <- 

Vln_GutCandidates <- VlnPlot(Naive.gutblood.merged.all, features=candidate_genes,
                             split.by = 'age.ident', pt.size=0.00001) & 
  scale_fill_manual(values=c('pink','black'))&
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
Feature_GutCandidates <- FeaturePlot(Naive.gutblood.merged.all, features=candidate_genes, pt.size=1) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))

Vln_GutCandidates_organ <- VlnPlot(Naive.gutblood.merged.all_not2, features=candidate_genes, 
                                   split.by = 'age.ident',pt.size=0.00001, group.by = 'group.ident') & 
  scale_fill_manual(values=c('pink','black'))&
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))

ggsave(plot=Vln_GutCandidates, filename='CandidateGenes_Vlnplot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)
ggsave(plot=Vln_GutCandidates_organ, filename='CandidateGenes_Vlnplot_organ_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)
ggsave(plot=Feature_GutCandidates, filename='CandidateGenes_FeaturePlot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)

ggsave(plot=Vln_GutCandidates, filename='CandidateGenes_Vlnplot_Gutenriched_NaiveTcells_childrenadults_byage.pdf', height = 16, width=18)
ggsave(plot=Vln_GutCandidates_organ, filename='CandidateGenes_Vlnplot_organ_Gutenriched_NaiveTcells_childrenadults_byage.pdf', height = 16, width=18)
ggsave(plot=Feature_GutCandidates, filename='CandidateGenes_FeaturePlot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)

ggsave(plot=Vln_GutCandidates, filename='CandidateGenes_Vlnplot_Gutenriched_NaiveTcells_childrenadults_byorgan.pdf', height = 16, width=18)


#### heatmap
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
load('ListDE_Naive.gutblood.merged.all.SCTonmtrb_res0.6_allmarkers_MASTRNAdata_childrenadults.Robj')
top75 <- lapply(ListDE_Naive.gutblood.merged.all_res0.6_allmarkers_MASTRNAdata, function(df) {df %>% subset(avg_log2FC>0) %>% top_n(75, avg_log2FC)})
top75 <- bind_rows(top75)

suppressPackageStartupMessages({
  library(rlang)
})
DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               additional.cols.use = NULL,  # Add this parameter for additional group colors
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by != i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
          # Use the colors for the main group (group.by)
          cols[[colname]] <- cols.use[[colname]] %||% c('black','lightblue','orange')  
        } else {
          colid = colname
          # Use the colors for the additional group (additional.group.by)
          cols[[colname]] <- additional.cols.use[[colname]] %||% c('blue','orange')  
        }
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

Naive.gutblood.merged.all <- ScaleData(Naive.gutblood.merged.all)
c <- DoMultiBarHeatmap(
  Naive.gutblood.merged.all,
  additional.group.by = "group.ident",
  additional.group.sort.by = "group.ident",
  features = top75$gene,
  cells = NULL,
  group.by = c("SCT_snn_res.0.6"),
  label=F
) + scale_color_manual(values=c('black','lightblue','orange'))
c

ggsave(filename='heatmap_childrenadults_res.0.6.pdf',plot=c, width=10,height=12)

#### Differential gene expression - children and adult ####
DefaultAssay(object = Naive.gutblood.merged.all) <- "RNA"
Naive.gutblood.merged.all <- JoinLayers(Naive.gutblood.merged.all)

Naive.gutblood.merged.all@active.ident <- Naive.gutblood.merged.all$SCT_snn_res.0.7

##Find Markers that are specific for each cluster
Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata=FindAllMarkers(Naive.gutblood.merged.all, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                                   min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)

## Create list
ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata<- split(Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata, 
                                                                   f=Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata$cluster)
## Filter on adj.P-value
ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata <-lapply(ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata<-lapply(ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

##save as Robj
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
save(ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata, file='ListDE_Naive.gutblood.merged.all.SCTonmtrb_res0.7_allmarkers_MASTRNAdata_childrenadults.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(ListDE_Naive.gutblood.merged.all_res0.7_allmarkers_MASTRNAdata, file='ListDE_Naive.gutblood.merged.all.SCTonmtrb_res0.7_allmarkers_MASTRNAdata_childrenadults.xlsx')

## compare 2 to the others
mapping <- c('0' = '0', '1' = '0', '2' = '1')
old_clusters <- as.character(Naive.gutblood.merged.all$SCT_snn_res.0.6)
new_clusters <- mapping[old_clusters]
names(new_clusters) <- colnames(Naive.gutblood.merged.all)
Naive.gutblood.merged.all$res.0.6_manual <- factor(new_clusters)

DefaultAssay(object = Naive.gutblood.merged.all) <- "RNA"
Naive.gutblood.merged.all <- JoinLayers(Naive.gutblood.merged.all)

Naive.gutblood.merged.all@active.ident <- Naive.gutblood.merged.all$res.0.6_manual

##Find Markers that are specific for each cluster
Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata=FindMarkers(Naive.gutblood.merged.all, ident.1='1', ident.2='0',
                                                                        test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                                        min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)

Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata$gene <- rownames(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata)
##optional: add pct.fold = how large is the absolute difference in percentage?
Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata$pct.fold <- Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata$pct.1/Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata$pct.2

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata <-dplyr::filter(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata, p_val_adj<0.05)
## Sort on logFC
Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata<-Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata[order(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata$avg_log2FC, decreasing=T),]

##save as Robj
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
save(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata, file='Naive.gutblood.merged.all.SCTonmtrb_res0.6_manual_markers_MASTRNAdata_childrenadults_C2vsC10.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata, file='Naive.gutblood.merged.all.SCTonmtrb_res0.6_manual_markers_MASTRNAdata_childrenadults_C2vsC10.xlsx')


#### 3. Pseudobulk Naive - gut vs blood - children ####
aggregatedCounts <- AggregateExpression(Naive.gutblood.merged, group.by=c('orig.ident'),
                                        assay='RNA', slot='counts', return.seurat = F)
aggregatedCounts <- as.data.frame(aggregatedCounts$RNA)
aggregatedCounts[1:8, ]

colData <- data.frame(row.names = colnames(aggregatedCounts), origID = colnames(aggregatedCounts), 
                      groupID=ifelse(grepl('EA024|EA026|EA029|EA030', colnames(aggregatedCounts)), 'Blood', 'Gut'))
colData$groupID <- as.factor(colData$groupID)

# colData$age <- c(31,31,27,27,41,41,41,41)
# colData$sex <- c('F','F','F','F','F','F','M','M')
# colData$seqrun <- c(1,1,1,1,1,1,2,2)
# colData$seqrun <- as.factor(colData$seqrun)
colData$donorID <- c('HINT139', 'HINT139','HINT140','HINT140', 'HINT148', 'HINT148', 'HINT149')
  
## run with and without without age
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts, colData = colData, design = ~ donorID+groupID)
keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

## first check on clustering/quality
## Transform counts for data visualization; blind=T means that design is not taken into account (good for QC check, for downstream things you might want to take all the info you have according to manual)
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=TRUE)

## Plot PCA
Fig_PCApseudo_Organ <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "groupID") +
  scale_color_manual(values=c('blue',  'orange'), name='')+
  ggtitle('Organ')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Organ

Fig_PCApseudo_Donor <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "donorID") +
  scale_color_manual(values=c('seagreen3','purple','skyblue','gold'), name='')+
  ggtitle('Donor variation')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Donor

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(plot=Fig_PCApseudo_Organ, filename = 'Naivegutblood_PCApseudobulk_Organ_children.pdf', height=6, width=6)
ggsave(plot=Fig_PCApseudo_Donor, filename = 'Naivegutblood_PCApseudobulk_Donor_children.pdf', height=6, width=6)

## DESeq2 analysis vs Control
dds$Organ <- relevel(dds$groupID, ref = "Blood")

dds <- DESeq(dds)

Fig_dispersionPseudo <- plotDispEsts(dds)
ggsave(plot=Fig_dispersionPseudo, filename = 'Naivegutblood_PCApseudobulk_dispersionplot_children.pdf')

resnames <- resultsNames(dds)[-1]
results <- list()
for(i in 1:length(resnames)){
  results[[i]] <- lfcShrink(dds, coef=resnames[i], type="apeglm")
}
names(results) <- resnames
results <- lapply(results, function(x){
  x <- as.data.frame(x)
  cbind(x, Gene=rownames(x))
})
results <- lapply(results, function(x){
  dplyr::filter(x, padj < 0.05) %>%
    dplyr::arrange(log2FoldChange)
})
results[[4]]

write_xlsx(results[[4]],
           path="T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels/DESeq2results_Naivegutblood_PCApseudobulk_paired.xlsx")

##visualizations
results_thres <- lapply(results, function(x) {
  x_filtered <- x[!is.na(x$padj),] %>%
    mutate(threshold = padj < 0.05)
  
  # Identify rows passing any slice_min/max selection
  selected_rows <- bind_rows(
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange)),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange))
  ) %>% distinct()
  extra_rows <-  x_filtered %>% filter(Gene %in% c("LEF1", "ITGA6", "GIMAP4", "KLF2", "SELL", "BTG1","BTG2","CXCR4","SOX4","IL7R","FOXP1","KLF6","ICOS")) # Corrected gene selection
  
  print(rownames(selected_rows)) # Debugging
  print(rownames(extra_rows))
  # Add a boolean column indicating selection status
  x_filtered <- x_filtered %>%
    mutate(selected = Gene %in% selected_rows$Gene, 
          selected_plus = (Gene %in% selected_rows$Gene)|(Gene %in% extra_rows$Gene),
           selected_names = selected_plus & !str_detect(Gene, "^(AC[0-9]*|RP|NCB|PIM1|LCN2|RACK1|TPT1)"))
  
  return(x_filtered)
})


## Generate plot
Fig_Volcano_FSIvsFSP <- ggplot(results_thres[[4]], aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = selected)) +
  ggtitle('Naive T cells - Pseudobulk Blood vs Gut',
          subtitle = 'Blood                                                                                                            Gut') +
  geom_text_repel(aes(label=ifelse(selected_names==T,as.character(Gene),''),
                      x = log2FoldChange, y = -log10(padj)), size=6.5, direction='both',force_pull=0.5,nudge_y=1,
                  max.overlaps = 18)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  labs(caption = 'Red = top25 p-adj and log2FC')+
  scale_color_manual(values = c("grey60", "red3")) +
  scale_x_continuous(limits = c(-5.6,10.8))+
  theme(legend.position = "none",
        axis.text = element_text(size=1.5*20), axis.title = element_text(size=1.5*22),
        plot.title = element_text(size=1.5*30, hjust=0.5), plot.subtitle = element_text(size=1.5*28, hjust=0.5),
        plot.caption = element_text(size=1.5*20),plot.caption.position = 'panel',
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        plot.margin = margin(20,10,10,10))
Fig_Volcano_FSIvsFSP

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(Fig_Volcano_FSIvsFSP,filename=('Naivegutblood_PCApseudobulk_paired_Volcanopseudobulk_2_children.pdf'), height=14, width=28)

#### 3.2 Pseudobulk Naive - gut vs blood - children - not the gut-cluster 2 ####
Naive.gutblood.merged <- SetIdent(Naive.gutblood.merged,value='SCT_snn_res.0.8')
Naive.gutblood.merged_not2 <- subset(Naive.gutblood.merged, idents=c('0','1'))
aggregatedCounts <- AggregateExpression(Naive.gutblood.merged_not2, group.by=c('orig.ident'),
                                        assay='RNA', slot='counts', return.seurat = F)
aggregatedCounts <- as.data.frame(aggregatedCounts$RNA)
aggregatedCounts[1:8, ]

colData <- data.frame(row.names = colnames(aggregatedCounts), origID = colnames(aggregatedCounts), 
                      groupID=ifelse(grepl('EA024|EA026|EA029|EA030', colnames(aggregatedCounts)), 'Blood', 'Gut'))
colData$groupID <- as.factor(colData$groupID)

# colData$age <- c(31,31,27,27,41,41,41,41)
# colData$sex <- c('F','F','F','F','F','F','M','M')
# colData$seqrun <- c(1,1,1,1,1,1,2,2)
# colData$seqrun <- as.factor(colData$seqrun)
colData$donorID <- c('HINT139', 'HINT139','HINT140','HINT140', 'HINT148', 'HINT148', 'HINT149')


dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts, colData = colData, design = ~ donorID+groupID)

keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

## first check on clustering/quality
## Transform counts for data visualization; blind=T means that design is not taken into account (good for QC check, for downstream things you might want to take all the info you have according to manual)
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=TRUE)

## Plot PCA
Fig_PCApseudo_Organ <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "groupID") +
  scale_color_manual(values=c('blue',  'orange'), name='')+
  ggtitle('Organ')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Organ

Fig_PCApseudo_Donor <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "donorID") +
  scale_color_manual(values=c('seagreen3','purple3','skyblue','gold','black','orange','cyan','magenta'), name='')+
  ggtitle('Donor variation')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Donor

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(plot=Fig_PCApseudo_Organ, filename = 'Naivegutblood_PCApseudobulk_Organ_children_excl.clust2.pdf', height=6, width=6)
ggsave(plot=Fig_PCApseudo_Donor, filename = 'Naivegutblood_PCApseudobulk_Donor_children_excl.clust2.pdf', height=6, width=6)

## DESeq2 analysis vs Control
dds$Organ <- relevel(dds$groupID, ref = "Blood")

dds <- DESeq(dds)

Fig_dispersionPseudo <- plotDispEsts(dds)
ggsave(plot=Fig_dispersionPseudo, filename = 'Naivegutblood_PCApseudobulk_dispersionplot_children.pdf')

resnames <- resultsNames(dds)[-1]
results <- list()
for(i in 1:length(resnames)){
  results[[i]] <- lfcShrink(dds, coef=resnames[i], type="apeglm")
}
names(results) <- resnames
results <- lapply(results, function(x){
  x <- as.data.frame(x)
  cbind(x, Gene=rownames(x))
})
results <- lapply(results, function(x){
  dplyr::filter(x, padj < 0.05) %>%
    dplyr::arrange(log2FoldChange)
})
results[[4]]

write_xlsx(results[[4]],
           path="T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels/DESeq2results_Naivegutblood_PCApseudobulk_paired_children_excl.clust2.xlsx")

##visualizations
results_thres <- lapply(results, function(x) {
  x_filtered <- x[!is.na(x$padj),] %>%
    mutate(threshold = padj < 0.05)
  
  # Identify rows passing any slice_min/max selection
  selected_rows <- bind_rows(
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange)),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange))
  ) %>% distinct()
  extra_rows <-  x_filtered %>% filter(Gene %in% c("LEF1", "ITGA6", "GIMAP4", "KLF2", "SELL", "BTG1","BTG2","CXCR4","CD38","SOX4","IL7R","FOXP1","KLF6","ICOS", "SOCS3")) # Corrected gene selection
  
  print(rownames(selected_rows)) # Debugging
  print(rownames(extra_rows))
  # Add a boolean column indicating selection status
  x_filtered <- x_filtered %>%
    mutate(selected = Gene %in% selected_rows$Gene, 
           selected_plus = (Gene %in% selected_rows$Gene)|(Gene %in% extra_rows$Gene),
           selected_names = selected_plus & !str_detect(Gene, "^(AC[0-9]*|RP|NCB|PIM1|LCN2|RACK1|TPT1|LINC|LRR|SNHG)"))
  
  return(x_filtered)
})


## Generate plot
Fig_Volcano_FSIvsFSP <- ggplot(results_thres[[4]], aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = selected)) +
  ggtitle('Naive T cells - Pseudobulk Blood vs Gut - Children - gut cluster 2 excluded',
          subtitle = 'Blood                                                                                                                     Gut') +
  geom_text_repel(aes(label=ifelse(selected_names==T,as.character(Gene),''),
                      x = log2FoldChange, y = -log10(padj)), size=7.5, direction='both',force_pull=0.5,nudge_y=1,
                  max.overlaps = 25)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  labs(caption = 'Red = top25 p-adj and log2FC')+
  scale_color_manual(values = c("grey60", "red3")) +
  theme(legend.position = "none",
        axis.text = element_text(size=1.5*20), axis.title = element_text(size=1.5*22),
        plot.title = element_text(size=1.5*30, hjust=0.5), plot.subtitle = element_text(size=1.5*28, hjust=0.5),
        plot.caption = element_text(size=1.5*20),plot.caption.position = 'panel',
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        plot.margin = margin(20,10,10,10))
Fig_Volcano_FSIvsFSP

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(Fig_Volcano_FSIvsFSP,filename=('Naivegutblood_PCApseudobulk_paired_Volcanopseudobulk_children_excl.clust2.pdf'), height=14, width=28)


#### 3.3 Pseudobulk Naive - gut vs blood - children and adults together ####
aggregatedCounts <- AggregateExpression(Naive.gutblood.merged.all, group.by=c('orig.ident'),
                                        assay='RNA', slot='counts', return.seurat = F)
aggregatedCounts <- as.data.frame(aggregatedCounts$RNA)
aggregatedCounts[1:8, ]

colData <- data.frame(row.names = colnames(aggregatedCounts), origID = colnames(aggregatedCounts), 
                      ageID=ifelse(grepl('EA012|EA013|EA014|EA015|EA016|EA017|EA021|EA022', colnames(aggregatedCounts)), 'Children', 'Adult'),
                      groupID=ifelse(grepl('EA013|EA015|EA017|EA022|EA024|EA026|EA029|EA030', colnames(aggregatedCounts)), 'Blood', 'Gut'))
colData$ageID <- as.factor(colData$ageID)
colData$groupID <- as.factor(colData$groupID)

# colData$age <- c(31,31,27,27,41,41,41,41)
# colData$sex <- c('F','F','F','F','F','F','M','M')
# colData$seqrun <- c(1,1,1,1,1,1,2,2)
# colData$seqrun <- as.factor(colData$seqrun)
colData$donorID <- c('HINT129','HINT129','HINT130','HINT130','HINT131','HINT131','HINT136','HINT136',
                     'HINT139', 'HINT139','HINT140','HINT140', 'HINT148', 'HINT148', 'HINT149')
## run with and without without age
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts[,1:8], colData = colData[1:8,], design = ~ donorID+groupID)
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts[,9:15], colData = colData[9:15,], design = ~ donorID+groupID)
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts, colData = colData, design = ~ donorID+groupID)

keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

## first check on clustering/quality
## Transform counts for data visualization; blind=T means that design is not taken into account (good for QC check, for downstream things you might want to take all the info you have according to manual)
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=TRUE)

## Plot PCA
Fig_PCApseudo_Organ <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "groupID") +
  scale_color_manual(values=c('blue',  'orange'), name='')+
  ggtitle('Organ')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Organ

Fig_PCApseudo_Donor <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "donorID") +
  scale_color_manual(values=c('seagreen3','purple','skyblue','gold'), name='')+
  ggtitle('Donor variation')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Donor

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(plot=Fig_PCApseudo_Organ, filename = 'Naivegutblood_PCApseudobulk_Organ_children.pdf', height=6, width=6)
ggsave(plot=Fig_PCApseudo_Donor, filename = 'Naivegutblood_PCApseudobulk_Donor_children.pdf', height=6, width=6)

## DESeq2 analysis vs Control
dds$Organ <- relevel(dds$groupID, ref = "Blood")

dds <- DESeq(dds)

Fig_dispersionPseudo <- plotDispEsts(dds)
ggsave(plot=Fig_dispersionPseudo, filename = 'Naivegutblood_PCApseudobulk_dispersionplot_children.pdf')

resnames <- resultsNames(dds)[-1]
results <- list()
for(i in 1:length(resnames)){
  results[[i]] <- lfcShrink(dds, coef=resnames[i], type="apeglm")
}
names(results) <- resnames
results <- lapply(results, function(x){
  x <- as.data.frame(x)
  cbind(x, Gene=rownames(x))
})
results <- lapply(results, function(x){
  dplyr::filter(x, padj < 0.05) %>%
    dplyr::arrange(log2FoldChange)
})
results[[4]]

write_xlsx(results[[4]],
           path="T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels/DESeq2results_Naivegutblood_PCApseudobulk_paired_children.xlsx")

##visualizations
results_thres <- lapply(results, function(x) {
  x_filtered <- x[!is.na(x$padj),] %>%
    mutate(threshold = padj < 0.05)
  
  # Identify rows passing any slice_min/max selection
  selected_rows <- bind_rows(
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange)),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange))
  ) %>% distinct()
  extra_rows <-  x_filtered %>% filter(Gene %in% c("LEF1", "ITGA6", "GIMAP4", "KLF2", "SELL", "BTG1","BTG2","CXCR4","CD38","SOX4","IL7R","FOXP1","KLF6","ICOS", "SOCS3")) # Corrected gene selection
  
  print(rownames(selected_rows)) # Debugging
  print(rownames(extra_rows))
  # Add a boolean column indicating selection status
  x_filtered <- x_filtered %>%
    mutate(selected = Gene %in% selected_rows$Gene, 
           selected_plus = (Gene %in% selected_rows$Gene)|(Gene %in% extra_rows$Gene),
           selected_names = selected_plus & !str_detect(Gene, "^(AC[0-9]*|RP|NCB|PIM1|LCN2|RACK1|TPT1|LINC|LRR|SNHG)"))
  
  return(x_filtered)
})


## Generate plot
Fig_Volcano_FSIvsFSP <- ggplot(results_thres[[8]], aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = selected)) +
  ggtitle('Naive T cells - Pseudobulk Blood vs Gut - Children&Adults',
          subtitle = 'Blood                                                                                                                     Gut') +
  geom_text_repel(aes(label=ifelse(selected_names==T,as.character(Gene),''),
                      x = log2FoldChange, y = -log10(padj)), size=7.5, direction='both',force_pull=0.5,nudge_y=1,
                  max.overlaps = 25)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  labs(caption = 'Red = top25 p-adj and log2FC')+
  scale_color_manual(values = c("grey60", "red3")) +
  theme(legend.position = "none",
        axis.text = element_text(size=1.5*20), axis.title = element_text(size=1.5*22),
        plot.title = element_text(size=1.5*30, hjust=0.5), plot.subtitle = element_text(size=1.5*28, hjust=0.5),
        plot.caption = element_text(size=1.5*20),plot.caption.position = 'panel',
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        plot.margin = margin(20,10,10,10))
Fig_Volcano_FSIvsFSP

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(Fig_Volcano_FSIvsFSP,filename=('Naivegutblood_PCApseudobulk_paired_Volcanopseudobulk_childrenadults.pdf'), height=14, width=28)

#### 3.4 Pseudobulk Naive - gut vs blood - children and adults together - not the gut-cluster 2 ####
Naive.gutblood.merged.all <- SetIdent(Naive.gutblood.merged.all,value='SCT_snn_res.0.6')
Naive.gutblood.merged.all_not2 <- subset(Naive.gutblood.merged.all, idents=c('0','1'))
aggregatedCounts <- AggregateExpression(Naive.gutblood.merged.all_not2, group.by=c('orig.ident'),
                                        assay='RNA', slot='counts', return.seurat = F)
aggregatedCounts <- as.data.frame(aggregatedCounts$RNA)
aggregatedCounts[1:8, ]

colData <- data.frame(row.names = colnames(aggregatedCounts), origID = colnames(aggregatedCounts), 
                      ageID=ifelse(grepl('EA012|EA013|EA014|EA015|EA016|EA017|EA021|EA022', colnames(aggregatedCounts)), 'Children', 'Adult'),
                      groupID=ifelse(grepl('EA013|EA015|EA017|EA022|EA024|EA026|EA029|EA030', colnames(aggregatedCounts)), 'Blood', 'Gut'))
colData$ageID <- as.factor(colData$ageID)
colData$groupID <- as.factor(colData$groupID)

# colData$age <- c(31,31,27,27,41,41,41,41)
# colData$sex <- c('F','F','F','F','F','F','M','M')
# colData$seqrun <- c(1,1,1,1,1,1,2,2)
# colData$seqrun <- as.factor(colData$seqrun)
colData$donorID <- c('HINT129','HINT129','HINT130','HINT130','HINT131','HINT131','HINT136','HINT136',
                     'HINT139', 'HINT139','HINT140','HINT140', 'HINT148', 'HINT148', 'HINT149')
## run with and without without age
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts[,1:8], colData = colData[1:8,], design = ~ donorID+groupID)
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts[,9:15], colData = colData[9:15,], design = ~ donorID+groupID)

dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts, colData = colData, design = ~ donorID+groupID)

keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

## first check on clustering/quality
## Transform counts for data visualization; blind=T means that design is not taken into account (good for QC check, for downstream things you might want to take all the info you have according to manual)
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=TRUE)

## Plot PCA
Fig_PCApseudo_Organ <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "groupID") +
  scale_color_manual(values=c('blue',  'orange'), name='')+
  ggtitle('Organ')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Organ

Fig_PCApseudo_Donor <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "donorID") +
  scale_color_manual(values=c('seagreen3','purple3','skyblue','gold','black','orange','cyan','magenta'), name='')+
  ggtitle('Donor variation')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Donor

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(plot=Fig_PCApseudo_Organ, filename = 'Naivegutblood_PCApseudobulk_Organ_childrenadults_excl.clust2.pdf', height=6, width=6)
ggsave(plot=Fig_PCApseudo_Donor, filename = 'Naivegutblood_PCApseudobulk_Donor_childrenadults_excl.clust2.pdf', height=6, width=6)

## DESeq2 analysis vs Control
dds$Organ <- relevel(dds$groupID, ref = "Blood")

dds <- DESeq(dds)

Fig_dispersionPseudo <- plotDispEsts(dds)
ggsave(plot=Fig_dispersionPseudo, filename = 'Naivegutblood_PCApseudobulk_dispersionplot_children.pdf')

resnames <- resultsNames(dds)[-1]
results <- list()
for(i in 1:length(resnames)){
  results[[i]] <- lfcShrink(dds, coef=resnames[i], type="apeglm")
}
names(results) <- resnames
results <- lapply(results, function(x){
  x <- as.data.frame(x)
  cbind(x, Gene=rownames(x))
})
results <- lapply(results, function(x){
  dplyr::filter(x, padj < 0.05) %>%
    dplyr::arrange(log2FoldChange)
})
results[[8]]

write_xlsx(results[[8]],
           path="T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels/DESeq2results_Naivegutblood_PCApseudobulk_paired_childrenadults_excl.clust2.xlsx")

##visualizations
results_thres <- lapply(results, function(x) {
  x_filtered <- x[!is.na(x$padj),] %>%
    mutate(threshold = padj < 0.05)
  
  # Identify rows passing any slice_min/max selection
  selected_rows <- bind_rows(
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange)),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_min(n = 25, order_by = padj),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_max(n = 25, order_by = abs(log2FoldChange))
  ) %>% distinct()
  extra_rows <-  x_filtered %>% filter(Gene %in% c("LEF1", "ITGA6", "GIMAP4", "KLF2", "SELL", "BTG1","BTG2","CXCR4","CD38","SOX4","IL7R","FOXP1","KLF6","ICOS", "SOCS3")) # Corrected gene selection
  
  print(rownames(selected_rows)) # Debugging
  print(rownames(extra_rows))
  # Add a boolean column indicating selection status
  x_filtered <- x_filtered %>%
    mutate(selected = Gene %in% selected_rows$Gene, 
           selected_plus = (Gene %in% selected_rows$Gene)|(Gene %in% extra_rows$Gene),
           selected_names = selected_plus & !str_detect(Gene, "^(AC[0-9]*|RP|NCB|PIM1|LCN2|RACK1|TPT1|LINC|LRR|SNHG)"))
  
  return(x_filtered)
})


## Generate plot
Fig_Volcano_FSIvsFSP <- ggplot(results_thres[[8]], aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = selected)) +
  ggtitle('Naive T cells - Pseudobulk Blood vs Gut - Children&Adults - gut cluster 2 excluded',
          subtitle = 'Blood                                                                                                                     Gut') +
  geom_text_repel(aes(label=ifelse(selected_names==T,as.character(Gene),''),
                      x = log2FoldChange, y = -log10(padj)), size=7.5, direction='both',force_pull=0.5,nudge_y=1,
                  max.overlaps = 25)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  labs(caption = 'Red = top25 p-adj and log2FC')+
  scale_color_manual(values = c("grey60", "red3")) +
  theme(legend.position = "none",
        axis.text = element_text(size=1.5*20), axis.title = element_text(size=1.5*22),
        plot.title = element_text(size=1.5*30, hjust=0.5), plot.subtitle = element_text(size=1.5*28, hjust=0.5),
        plot.caption = element_text(size=1.5*20),plot.caption.position = 'panel',
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        plot.margin = margin(20,10,10,10))
Fig_Volcano_FSIvsFSP

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(Fig_Volcano_FSIvsFSP,filename=('Naivegutblood_PCApseudobulk_paired_Volcanopseudobulk_childrenadults_excl.clust2.pdf'), height=14, width=28)


#### 3.5 Pseudobulk Naive - adult vs children ####
aggregatedCounts <- AggregateExpression(Naive.gutblood.merged.all, group.by=c('orig.ident'),
                                        assay='RNA', slot='counts', return.seurat = F)
aggregatedCounts <- as.data.frame(aggregatedCounts$RNA)
aggregatedCounts[1:8, ]

colData <- data.frame(row.names = colnames(aggregatedCounts), origID = colnames(aggregatedCounts), 
                      ageID=ifelse(grepl('EA012|EA013|EA014|EA015|EA016|EA017|EA021|EA022', colnames(aggregatedCounts)), 'Children', 'Adult'),
                      groupID=ifelse(grepl('EA013|EA015|EA017|EA022|EA024|EA026|EA029|EA030', colnames(aggregatedCounts)), 'Blood', 'Gut'))
colData$ageID <- as.factor(colData$ageID)
colData$groupID <- as.factor(colData$groupID)

colData$age <- c(31,31,27,27,41,41,41,41)
colData$sex <- c('F','F','F','F','F','F','M','M')
colData$seqrun <- c(1,1,1,1,1,1,2,2)
colData$seqrun <- as.factor(colData$seqrun)
colData$donorID <- c('HINT129','HINT129','HINT130','HINT130','HINT131','HINT131','HINT136','HINT136',
                     'HINT139', 'HINT139','HINT140','HINT140', 'HINT148', 'HINT148', 'HINT149')

## run with and without tissue
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts[,c(2,4,6,8,10,12,14,15)], colData = colData[c(2,4,6,8,10,12,14,15),], design = ~ ageID)
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts[,c(1,3,5,7,9,11,13)], colData = colData[c(1,3,5,7,9,11,13),], design = ~ ageID)
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts, colData = colData, design = ~ ageID)
dds <- DESeqDataSetFromMatrix(countData = aggregatedCounts, colData = colData, design = ~ groupID+ageID)

keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

## first check on clustering/quality
## Transform counts for data visualization; blind=T means that design is not taken into account (good for QC check, for downstream things you might want to take all the info you have according to manual)
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=TRUE)

## Plot PCA
Fig_PCApseudo_age <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "ageID") +
  scale_color_manual(values=c('pink',  'black'), name='')+
  ggtitle('Age')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_age

Fig_PCApseudo_Organ <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "groupID") +
  scale_color_manual(values=c('blue',  'orange'), name='')+
  ggtitle('Organ')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Organ

Fig_PCApseudo_Donor <- DESeq2::plotPCA(vsd, ntop = 500, intgroup = "donorID") +
  scale_color_manual(values=c('seagreen3','purple3','skyblue','gold','black','orange','cyan','magenta'), name='')+
  ggtitle('Donor variation')+
  theme(plot.title = element_text(hjust=0.5, vjust=2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=.5))
Fig_PCApseudo_Donor

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(plot=Fig_PCApseudo_age, filename = 'Naivegutblood_PCApseudobulk_Age_childrenadult.pdf', height=6, width=6)
ggsave(plot=Fig_PCApseudo_Organ, filename = 'Naivegutblood_PCApseudobulk_Organ_childrenadult.pdf', height=6, width=6)
ggsave(plot=Fig_PCApseudo_Donor, filename = 'Naivegutblood_PCApseudobulk_Donor_childrenadult.pdf', height=6, width=6)

## DESeq2 analysis vs Control
dds$Age <- relevel(dds$ageID, ref = "Adult")

dds <- DESeq(dds)

Fig_dispersionPseudo <- plotDispEsts(dds)
ggsave(plot=Fig_dispersionPseudo, filename = 'Naivegutblood_PCApseudobulk_dispersionplot_children.pdf')

resnames <- resultsNames(dds)[-1]
results <- list()
for(i in 1:length(resnames)){
  results[[i]] <- lfcShrink(dds, coef=resnames[i], type="apeglm")
}
names(results) <- resnames
results <- lapply(results, function(x){
  x <- as.data.frame(x)
  cbind(x, Gene=rownames(x))
})
results <- lapply(results, function(x){
  dplyr::filter(x, padj < 0.05) %>%
    dplyr::arrange(log2FoldChange)
})
results[[1]]

write_xlsx(results[[1]],
           path="T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels/DESeq2results_NaiveGutBloodtogether_ChildrenvsAdult.xlsx")

##visualizations
results_thres <- lapply(results, function(x) {
  x_filtered <- x[!is.na(x$padj),] %>%
    mutate(threshold = padj < 0.05)
  
  # Identify rows passing any slice_min/max selection
  selected_rows <- bind_rows(
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_min(n = 35, order_by = padj),
    x_filtered %>% filter(log2FoldChange > 0) %>% slice_max(n = 35, order_by = abs(log2FoldChange)),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_min(n = 35, order_by = padj),
    x_filtered %>% filter(log2FoldChange < 0) %>% slice_max(n = 35, order_by = abs(log2FoldChange))
  ) %>% distinct()
  extra_rows <-  x_filtered %>% filter(Gene %in% c("LEF1", "ITGA6","STAT4", "GIMAP4", "KLF2", "SELL", "BTG1","BTG2","CXCR4","SOX4","IL7R","FOXP1","KLF6","ICOS")) # Corrected gene selection
  
  print(rownames(selected_rows)) # Debugging
  print(rownames(extra_rows))
  # Add a boolean column indicating selection status
  x_filtered <- x_filtered %>%
    mutate(selected = Gene %in% selected_rows$Gene, 
           selected_plus = (Gene %in% selected_rows$Gene)|(Gene %in% extra_rows$Gene),
           selected_names = selected_plus & !str_detect(Gene, "^(AC[0-9]*|RP|NCB|PIM1|LCN2|RACK1|TPT1)"))
  
  return(x_filtered)
})


## Generate plot
Fig_Volcano_FSIvsFSP <- ggplot(results_thres[[1]], aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = selected)) +
  ggtitle('Naive T cells - Pseudobulk Adult vs Children - Gut', 
          subtitle = 'Adults                                                                                                     Children') +
  geom_text_repel(aes(label=ifelse(selected_names==T,as.character(Gene),''),
                      x = log2FoldChange, y = -log10(padj)), size=8.5, direction='both',force_pull=0.5,nudge_y=1,
                  max.overlaps = 10)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  labs(caption = 'Red = top25 p-adj and log2FC')+
  scale_color_manual(values = c("grey60", "red3")) +
  theme(legend.position = "none",
        axis.text = element_text(size=1.5*20), axis.title = element_text(size=1.5*22),
        plot.title = element_text(size=1.5*30, hjust=0.5), plot.subtitle = element_text(size=1.5*28, hjust=0.5),
        plot.caption = element_text(size=1.5*20),plot.caption.position = 'panel',
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        plot.margin = margin(20,10,10,10))
Fig_Volcano_FSIvsFSP

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
ggsave(Fig_Volcano_FSIvsFSP,filename=('Naivegutblood_PCApseudobulk_paired_Volcanopseudobulk_gut_ChildrenvsAdult.pdf'), height=14, width=28)


##### 4. Gut - adult children combined ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/")
load('EA012SCT_pctmtribo_nocellcylce.Robj')
load('EA014SCT_pctmtribo_nocellcylce.Robj')
load('EA016SCT_pctmtribo_nocellcylce.Robj')

load('EA021SCT_pctmtribo_nocellcylce.Robj')

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/")
load('EA013SCT_pctmtribo_nocellcylce.Robj')
load('EA015SCT_pctmtribo_nocellcylce.Robj')
load('EA017SCT_pctmtribo_nocellcylce.Robj')
load('EA022SCT_pctmtribo_nocellcylce.Robj')

### 
Naive.gut.merged.all <- merge(EA012_mtrb, c(EA014_mtrb,  EA016_mtrb,  EA021_mtrb, 
                                                 EA023_mtrb,  EA025_mtrb, EA028_mtrb))


setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
save(Naive.gut.merged.all, file='Naive.gut.merged.all_SCTmtrb_childrenadults.Robj')
load(file = "T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Naive.gut.merged.all_SCTmtrb_childrenadults.Robj")

## after merge, you have to set features again: variablefeatures are just all SCT features
Naive.gut.merged.all <- RunPCA(Naive.gut.merged.all, features = rownames(Naive.gut.merged.all@assays[["SCT"]]@scale.data))
Naive.gut.merged.all <- RunUMAP(Naive.gut.merged.all, dims = 1:30)

##clustering
Naive.gut.merged.all <- FindNeighbors(object = Naive.gut.merged.all, reduction = "pca", dims = 1:30)
Naive.gut.merged.all <- FindClusters(object = Naive.gut.merged.all, resolution = seq(0,1.3,0.1))

##add group
Naive.gut.merged.all$group.ident <- Naive.gut.merged.all$orig.ident
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$orig.ident=='EA013'|Naive.gut.merged.all$ orig.ident=='EA015'|
                                      Naive.gut.merged.all$ orig.ident=='EA017'|Naive.gut.merged.all$ orig.ident=='EA022'|
                                      Naive.gut.merged.all$orig.ident=='EA024'|Naive.gut.merged.all$orig.ident=='EA026'|
                                      Naive.gut.merged.all$orig.ident=='EA029'|Naive.gut.merged.all$orig.ident=='EA030',
                                    'group.ident'] <- 'Blood'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$ orig.ident=='EA012'|Naive.gut.merged.all$ orig.ident=='EA014'|
                                      Naive.gut.merged.all$ orig.ident=='EA016'|Naive.gut.merged.all$ orig.ident=='EA021'|
                                      Naive.gut.merged.all$orig.ident=='EA023'|Naive.gut.merged.all$orig.ident=='EA025'|
                                      Naive.gut.merged.all$orig.ident=='EA028',
                                    'group.ident'] <- 'Gut'
Naive.gut.merged.all$donor.ident <- Naive.gut.merged.all$orig.ident
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$orig.ident=='EA023'|Naive.gut.merged.all$orig.ident=='EA024',
                                    'donor.ident'] <- 'HINT139'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$orig.ident=='EA025'|Naive.gut.merged.all$orig.ident=='EA026',
                                    'donor.ident'] <- 'HINT140'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$orig.ident=='EA028'|Naive.gut.merged.all$orig.ident=='EA029',
                                    'donor.ident'] <- 'HINT148'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$orig.ident=='EA030',
                                    'donor.ident'] <- 'HINT149'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$ orig.ident=='EA012'|Naive.gut.merged.all$ orig.ident=='EA013',
                                    'donor.ident'] <- 'HINT129'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$ orig.ident=='EA014'|Naive.gut.merged.all$ orig.ident=='EA015',
                                    'donor.ident'] <- 'HINT130'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$ orig.ident=='EA016'|Naive.gut.merged.all$ orig.ident=='EA017',
                                    'donor.ident'] <- 'HINT131'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$ orig.ident=='EA021'|Naive.gut.merged.all$ orig.ident=='EA022',
                                    'donor.ident'] <- 'HINT136'
Naive.gut.merged.all$age.ident <- Naive.gut.merged.all$orig.ident
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$orig.ident=='EA012'|Naive.gut.merged.all$ orig.ident=='EA013'|
                                      Naive.gut.merged.all$ orig.ident=='EA014'|Naive.gut.merged.all$ orig.ident=='EA015'|
                                      Naive.gut.merged.all$orig.ident=='EA016'|Naive.gut.merged.all$orig.ident=='EA017'|
                                      Naive.gut.merged.all$orig.ident=='EA021'|Naive.gut.merged.all$orig.ident=='EA022',
                                    'age.ident'] <- 'Children'
Naive.gut.merged.all@meta.data[Naive.gut.merged.all$ orig.ident=='EA023'|Naive.gut.merged.all$ orig.ident=='EA024'|
                                      Naive.gut.merged.all$ orig.ident=='EA025'|Naive.gut.merged.all$ orig.ident=='EA026'|
                                      Naive.gut.merged.all$orig.ident=='EA028'|Naive.gut.merged.all$orig.ident=='EA029'|
                                      Naive.gut.merged.all$orig.ident=='EA030',
                                    'age.ident'] <- 'Adults'

##visualization
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
DefaultAssay(Naive.gut.merged.all) <- 'RNA'

Fig_clustree.merged <- clustree(Naive.gut.merged.all, prefix='SCT_snn_res.')
ggsave(plot=Fig_clustree.merged, file='gut_Merged_Clustree_SCTonmtrb_childrenadults.pdf', width=8, height=12)

Fig_Dimplot.merged_OrigID <- DimPlot(object = Naive.gut.merged.all, reduction = "umap", group.by = 'orig.ident', pt.size=0.9)+
  scale_color_manual(values=c('pink','green', 'violet',  'seagreen3','purple4','seagreen1','lightblue'))+
  ggtitle('Sample ID', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_groupID <- DimPlot(object = Naive.gut.merged.all, reduction = "umap", group.by = 'group.ident', pt.size=0.9)+
  scale_color_manual(values=c('blue3',  'orange'))+
  ggtitle('Gut vs Blood', sub= 'Children&Adult Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_ageID <- DimPlot(object = Naive.gut.merged.all, reduction = "umap", group.by = 'age.ident', pt.size=0.9)+
  scale_color_manual(values=c('violet',  'black'))+
  ggtitle('Adult vs Children', sub= 'Children&Adult Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_DonorID <- DimPlot(object = Naive.gut.merged.all, reduction = "umap", group.by = 'donor.ident', pt.size=0.9)+
  scale_color_manual(values=c( 'seagreen3','purple3','skyblue','gold','black','orange','cyan','magenta'))+
  ggtitle('Donor ID', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_Dimplot.merged_res0.6 <- DimPlot(object = Naive.gut.merged.all, group.by = 'SCT_snn_res.0.6', reduction = "umap", pt.size=0.9,
                                     cols=c('black','skyblue2'))+
  ggtitle('Clustering resolution 0.6', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_res0.9 <- DimPlot(object = Naive.gut.merged.all,group.by = 'SCT_snn_res.0.9', reduction = "umap", pt.size=0.9,
                                     cols=c('black','skyblue2','orange'))+
  ggtitle('Clustering resolution 0.9', sub= 'SCT on pct.mtrb')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5))

Fig_Dimplot.merged_Phase <- DimPlot(object = Naive.gut.merged.all, reduction = "umap", group.by = 'Phase', pt.size=0.9)+
  scale_color_colorblind()
Fig_Featureplot.merged_nFeature <- FeaturePlot(Naive.gut.merged.all, features='nFeature_RNA', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_mtrb <- FeaturePlot(Naive.gut.merged.all, features='percent.mt', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_rb <- FeaturePlot(Naive.gut.merged.all, features='percent.ribo', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Naive.gut.merged.all$flowCD31_bin <- Naive.gut.merged.all$flowCD31 > 300
Fig_Featureplot.merged_CD31 <- DimPlot(Naive.gut.merged.all, group.by ='flowCD31_bin', pt.size=0.9, order=T) + 
  scale_color_manual(values=c('orange2', 'skyblue')) + ggtitle('surface CD31+')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_clustering_ext <-  Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.9+Fig_Featureplot.merged_CD31+
  Fig_Featureplot.merged_mtrb+ Fig_Featureplot.merged_rb + Fig_Featureplot.merged_nFeature+Fig_Dimplot.merged_Phase+
  plot_layout(ncol=5)

Fig_clustering <- Fig_Dimplot.merged_DonorID + Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.9+
  plot_layout(ncol=4)

ggsave('gut_Merged_ClusteringOverview_SCTonmtrb_childrenadults.pdf',Fig_clustering_ext, width=24, height=10)
ggsave('gut_Merged_ClusteringOverview_simple_SCTonmtrb_childrenadults.pdf',Fig_clustering, width=24, height=6)

Fig_Vlnplot.merged_CD31 <- VlnPlot(Naive.gut.merged.all, group.by = 'SCT_snn_res.0.6', feature='flowCD31')+
  VlnPlot(Naive.gut.merged.all, group.by = 'SCT_snn_res.0.8', feature='flowCD31')+
  VlnPlot(Naive.gut.merged.all, group.by = 'SCT_snn_res.0.9', feature='flowCD31')
ggsave('gut_Merged_CD31Vlns_SCTonmtrb_childrenadults.pdf',Fig_Vlnplot.merged_CD31, width=14, height=20)

Fig_barplot_res0.6 <- ggplot(Naive.gut.merged.all@meta.data, aes(x=age.ident, fill=SCT_snn_res.0.6)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.6 - SCT on pct.mtrb')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))

Fig_barplot_res0.9 <- ggplot(Naive.gut.merged.all@meta.data, aes(x=age.ident, fill=SCT_snn_res.0.9)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.9 - SCT on pct.mtrb  ')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))


Fig_barplots <- Fig_barplot_res0.6 + Fig_barplot_res0.9 + plot_layout(ncol=2)
ggsave(plot=Fig_barplots, filename='gut_Merged_Barplot_groupsperclusters_SCTonmtrb_childrenadults.pdf',width=16, height=5)


#### visualize potential gut-enriched naive T cell genes
DefaultAssay(object = Naive.gut.merged.all) <- "RNA"
Naive.gut.merged.all@active.ident <- Naive.gut.merged.all$SCT_snn_res.0.6
candidate_genes <- c('CXCR4','TNFAIP3','CTLA4','SOCS3', 
                     'KLF6', 'IL2RA', 'CCL5', 'ICOS',
                     'BTG1', 'LEF1','KLF2','SOX4',
                     'CREM','SRGN', 'SLC2A3',  'PTGER4',
                     'GPR183', 'CD69', 'LDLR', 'HIF1A')
FeaturePlot(Naive.gut.merged.all, features=candidate_genes)
VlnPlot(Naive.gut.merged.all, features=candidate_genes)

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
Vln_GutCandidates <- VlnPlot(Naive.gut.merged.all, features=candidate_genes, pt.size=0.00001) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
Feature_GutCandidates <- FeaturePlot(Naive.gut.merged.all, features=candidate_genes, pt.size=1) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
ggsave(plot=Vln_GutCandidates, filename='CandidateGenes_Vlnplot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)
ggsave(plot=Feature_GutCandidates, filename='CandidateGenes_FeaturePlot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)

##Differential gene expression
DefaultAssay(object = Naive.gut.merged.all) <- "RNA"
Naive.gut.merged.all <- JoinLayers(Naive.gut.merged.all)

Naive.gut.merged.all@active.ident <- Naive.gut.merged.all$SCT_snn_res.0.9

##Find Markers that are specific for each cluster
Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata=FindAllMarkers(Naive.gut.merged.all, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                                       min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)

## Create list
ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata<- split(Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata, 
                                                                       f=Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata$cluster)
## Filter on adj.P-value
ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata <-lapply(ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata<-lapply(ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

##save as Robj
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
save(ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata, file='ListDE_Naive.gut.merged.all.SCTonmtrb_res0.9_allmarkers_MASTRNAdata_childrenadults.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(ListDE_Naive.gut.merged.all_res0.9_allmarkers_MASTRNAdata, file='ListDE_Naive.gut.merged.all.SCTonmtrb_res0.9_allmarkers_MASTRNAdata_childrenadults.xlsx')

##### 10. Blood - adult children combined ####
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/")
load('EA012SCT_pctmtribo_nocellcylce.Robj')
load('EA014SCT_pctmtribo_nocellcylce.Robj')
load('EA016SCT_pctmtribo_nocellcylce.Robj')

load('EA021SCT_pctmtribo_nocellcylce.Robj')

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/")
load('EA013SCT_pctmtribo_nocellcylce.Robj')
load('EA015SCT_pctmtribo_nocellcylce.Robj')
load('EA017SCT_pctmtribo_nocellcylce.Robj')
load('EA022SCT_pctmtribo_nocellcylce.Robj')

### 
Naive.blood.merged.all <- merge(EA013_mtrb, c(EA015_mtrb,  EA017_mtrb,  EA022_mtrb, 
                                            EA024_mtrb,  EA026_mtrb, EA029_mtrb, EA030_mtrb))


setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
save(Naive.blood.merged.all, file='Naive.blood.merged.all_SCTmtrb_childrenadults.Robj')
load(file = "T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Naive.blood.merged.all_SCTmtrb_childrenadults.Robj")

## after merge, you have to set features again: variablefeatures are just all SCT features
Naive.blood.merged.all <- RunPCA(Naive.blood.merged.all, features = rownames(Naive.blood.merged.all@assays[["SCT"]]@scale.data))
Naive.blood.merged.all <- RunUMAP(Naive.blood.merged.all, dims = 1:30)

##clustering
Naive.blood.merged.all <- FindNeighbors(object = Naive.blood.merged.all, reduction = "pca", dims = 1:30)
Naive.blood.merged.all <- FindClusters(object = Naive.blood.merged.all, resolution = seq(0,1.3,0.1))

##add group
Naive.blood.merged.all$group.ident <- Naive.blood.merged.all$orig.ident
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$orig.ident=='EA013'|Naive.blood.merged.all$ orig.ident=='EA015'|
                                 Naive.blood.merged.all$ orig.ident=='EA017'|Naive.blood.merged.all$ orig.ident=='EA022'|
                                 Naive.blood.merged.all$orig.ident=='EA024'|Naive.blood.merged.all$orig.ident=='EA026'|
                                 Naive.blood.merged.all$orig.ident=='EA029'|Naive.blood.merged.all$orig.ident=='EA030',
                               'group.ident'] <- 'Blood'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$ orig.ident=='EA012'|Naive.blood.merged.all$ orig.ident=='EA014'|
                                 Naive.blood.merged.all$ orig.ident=='EA016'|Naive.blood.merged.all$ orig.ident=='EA021'|
                                 Naive.blood.merged.all$orig.ident=='EA023'|Naive.blood.merged.all$orig.ident=='EA025'|
                                 Naive.blood.merged.all$orig.ident=='EA028',
                               'group.ident'] <- 'Gut'
Naive.blood.merged.all$donor.ident <- Naive.blood.merged.all$orig.ident
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$orig.ident=='EA023'|Naive.blood.merged.all$orig.ident=='EA024',
                               'donor.ident'] <- 'HINT139'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$orig.ident=='EA025'|Naive.blood.merged.all$orig.ident=='EA026',
                               'donor.ident'] <- 'HINT140'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$orig.ident=='EA028'|Naive.blood.merged.all$orig.ident=='EA029',
                               'donor.ident'] <- 'HINT148'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$orig.ident=='EA030',
                               'donor.ident'] <- 'HINT149'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$ orig.ident=='EA012'|Naive.blood.merged.all$ orig.ident=='EA013',
                               'donor.ident'] <- 'HINT129'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$ orig.ident=='EA014'|Naive.blood.merged.all$ orig.ident=='EA015',
                               'donor.ident'] <- 'HINT130'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$ orig.ident=='EA016'|Naive.blood.merged.all$ orig.ident=='EA017',
                               'donor.ident'] <- 'HINT131'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$ orig.ident=='EA021'|Naive.blood.merged.all$ orig.ident=='EA022',
                               'donor.ident'] <- 'HINT136'
Naive.blood.merged.all$age.ident <- Naive.blood.merged.all$orig.ident
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$orig.ident=='EA012'|Naive.blood.merged.all$ orig.ident=='EA013'|
                                 Naive.blood.merged.all$ orig.ident=='EA014'|Naive.blood.merged.all$ orig.ident=='EA015'|
                                 Naive.blood.merged.all$orig.ident=='EA016'|Naive.blood.merged.all$orig.ident=='EA017'|
                                 Naive.blood.merged.all$orig.ident=='EA021'|Naive.blood.merged.all$orig.ident=='EA022',
                               'age.ident'] <- 'Children'
Naive.blood.merged.all@meta.data[Naive.blood.merged.all$ orig.ident=='EA023'|Naive.blood.merged.all$ orig.ident=='EA024'|
                                 Naive.blood.merged.all$ orig.ident=='EA025'|Naive.blood.merged.all$ orig.ident=='EA026'|
                                 Naive.blood.merged.all$orig.ident=='EA028'|Naive.blood.merged.all$orig.ident=='EA029'|
                                 Naive.blood.merged.all$orig.ident=='EA030',
                               'age.ident'] <- 'Adults'

##visualization
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
DefaultAssay(Naive.blood.merged.all) <- 'RNA'

Fig_clustree.merged <- clustree(Naive.blood.merged.all, prefix='SCT_snn_res.')
ggsave(plot=Fig_clustree.merged, file='gut_Merged_Clustree_SCTonmtrb_childrenadults.pdf', width=8, height=12)

Fig_Dimplot.merged_OrigID <- DimPlot(object = Naive.blood.merged.all, reduction = "umap", group.by = 'orig.ident', pt.size=0.9)+
  scale_color_manual(values=c('pink','green', 'violet',  'seagreen3','purple4','seagreen1','lightblue'))+
  ggtitle('Sample ID', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_groupID <- DimPlot(object = Naive.blood.merged.all, reduction = "umap", group.by = 'group.ident', pt.size=0.9)+
  scale_color_manual(values=c('blue3',  'orange'))+
  ggtitle('Gut vs Blood', sub= 'Children&Adult Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_ageID <- DimPlot(object = Naive.blood.merged.all, reduction = "umap", group.by = 'age.ident', pt.size=0.9)+
  scale_color_manual(values=c('violet',  'black'))+
  ggtitle('Adult vs Children', sub= 'Children&Adult Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_DonorID <- DimPlot(object = Naive.blood.merged.all, reduction = "umap", group.by = 'donor.ident', pt.size=0.9)+
  scale_color_manual(values=c( 'seagreen3','purple3','skyblue','gold','black','orange','cyan','magenta'))+
  ggtitle('Donor ID', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_Dimplot.merged_res0.5 <- DimPlot(object = Naive.blood.merged.all, group.by = 'SCT_snn_res.0.5', reduction = "umap", pt.size=0.9,
                                     cols=c('black','skyblue2'))+
  ggtitle('Clustering resolution 0.5', sub= 'childrenadults Naive T cells')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))

Fig_Dimplot.merged_res0.9 <- DimPlot(object = Naive.blood.merged.all,group.by = 'SCT_snn_res.0.9', reduction = "umap", pt.size=0.9,
                                     cols=c('black','skyblue2','orange'))+
  ggtitle('Clustering resolution 0.9', sub= 'SCT on pct.mtrb')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5))

Fig_Dimplot.merged_Phase <- DimPlot(object = Naive.blood.merged.all, reduction = "umap", group.by = 'Phase', pt.size=0.9)+
  scale_color_colorblind()
Fig_Featureplot.merged_nFeature <- FeaturePlot(Naive.blood.merged.all, features='nFeature_RNA', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_mtrb <- FeaturePlot(Naive.blood.merged.all, features='percent.mt', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Fig_Featureplot.merged_rb <- FeaturePlot(Naive.blood.merged.all, features='percent.ribo', pt.size=0.9)+
  scale_color_gradient(low = 'yellow', high='red')
Naive.blood.merged.all$flowCD31_bin <- Naive.blood.merged.all$flowCD31 > 300
Fig_Featureplot.merged_CD31 <- DimPlot(Naive.blood.merged.all, group.by ='flowCD31_bin', pt.size=0.9, order=T) + 
  scale_color_manual(values=c('orange2', 'skyblue')) + ggtitle('surface CD31+')+
  theme(plot.subtitle = element_text(hjust=0.5), plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=16))


Fig_clustering_ext <-  Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.9+Fig_Featureplot.merged_CD31+
  Fig_Featureplot.merged_mtrb+ Fig_Featureplot.merged_rb + Fig_Featureplot.merged_nFeature+Fig_Dimplot.merged_Phase+
  plot_layout(ncol=5)

Fig_clustering <- Fig_Dimplot.merged_DonorID + Fig_Dimplot.merged_groupID +Fig_Dimplot.merged_res0.6+Fig_Dimplot.merged_res0.9+
  plot_layout(ncol=4)

ggsave('blood_Merged_ClusteringOverview_SCTonmtrb_childrenadults.pdf',Fig_clustering_ext, width=24, height=10)
ggsave('blood_Merged_ClusteringOverview_simple_SCTonmtrb_childrenadults.pdf',Fig_clustering, width=24, height=6)

Fig_Vlnplot.merged_CD31 <- VlnPlot(Naive.blood.merged.all, group.by = 'SCT_snn_res.0.6', feature='flowCD31')+
  VlnPlot(Naive.blood.merged.all, group.by = 'SCT_snn_res.0.8', feature='flowCD31')+
  VlnPlot(Naive.blood.merged.all, group.by = 'SCT_snn_res.0.9', feature='flowCD31')
ggsave('blood_Merged_CD31Vlns_SCTonmtrb_childrenadults.pdf',Fig_Vlnplot.merged_CD31, width=14, height=20)

Fig_barplot_res0.5 <- ggplot(Naive.blood.merged.all@meta.data, aes(x=age.ident, fill=SCT_snn_res.0.5)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.5 - SCT on pct.mtrb')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))

Fig_barplot_res0.9 <- ggplot(Naive.blood.merged.all@meta.data, aes(x=age.ident, fill=SCT_snn_res.0.9)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab("Fraction") + 
  labs(fill = "Cluster", title='Proportionplot clusters/group', subtitle='res 0.9 - SCT on pct.mtrb  ')+
  theme(plot.title = element_text(hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5),
        axis.title.y = element_text(vjust=2.5), axis.text.x = element_text(size=10))+
  scale_y_continuous(expand = c(0,0))


Fig_barplots <- Fig_barplot_res0.5 + Fig_barplot_res0.9 + plot_layout(ncol=2)
ggsave(plot=Fig_barplots, filename='blood_Merged_Barplot_groupsperclusters_SCTonmtrb_childrenadults.pdf',width=16, height=5)


#### visualize potential gut-enriched naive T cell genes
DefaultAssay(object = Naive.blood.merged.all) <- "RNA"
Naive.blood.merged.all@active.ident <- Naive.blood.merged.all$SCT_snn_res.0.5
candidate_genes <- c('CXCR4','TNFAIP3','CTLA4','SOCS3', 
                     'KLF6', 'IL2RA', 'CCL5', 'ICOS',
                     'BTG1', 'LEF1','KLF2','SOX4',
                     'CREM','SRGN', 'SLC2A3',  'PTGER4',
                     'GPR183', 'CD69', 'LDLR', 'HIF1A')
FeaturePlot(Naive.blood.merged.all, features=candidate_genes)
VlnPlot(Naive.blood.merged.all, features=candidate_genes)

setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
Vln_GutCandidates <- VlnPlot(Naive.blood.merged.all, features=candidate_genes, pt.size=0.00001) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
Feature_GutCandidates <- FeaturePlot(Naive.blood.merged.all, features=candidate_genes, pt.size=1) & 
  theme(axis.title.x = element_blank(), axis.text = element_text(size=20), axis.title.y = element_text(size=22),
        plot.title = element_text(size=26))
ggsave(plot=Vln_GutCandidates, filename='CandidateGenes_Vlnplot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)
ggsave(plot=Feature_GutCandidates, filename='CandidateGenes_FeaturePlot_Gutenriched_NaiveTcells_childrenadults.pdf', height = 16, width=18)

##Differential gene expression
DefaultAssay(object = Naive.blood.merged.all) <- "RNA"
Naive.blood.merged.all <- JoinLayers(Naive.blood.merged.all)

Naive.blood.merged.all@active.ident <- Naive.blood.merged.all$SCT_snn_res.0.5

##Find Markers that are specific for each cluster
Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata=FindAllMarkers(Naive.blood.merged.all, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)

## Create list
ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata<- split(Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata, 
                                                                  f=Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata$cluster)
## Filter on adj.P-value
ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata <-lapply(ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata<-lapply(ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})

##save as Robj
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/Figures&Excels")
save(ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata, file='ListDE_Naive.blood.merged.all.SCTonmtrb_res0.5_allmarkers_MASTRNAdata_childrenadults.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(ListDE_Naive.blood.merged.all_res0.5_allmarkers_MASTRNAdata, file='ListDE_Naive.blood.merged.all.SCTonmtrb_res0.5_allmarkers_MASTRNAdata_childrenadults.xlsx')


