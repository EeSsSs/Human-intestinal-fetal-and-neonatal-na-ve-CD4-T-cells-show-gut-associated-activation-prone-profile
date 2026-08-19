#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
library(ggbreak)
library(Seurat)
library(sctransform)
library(biomaRt)
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
library(monocle3)
library(SeuratWrappers)
library(tradeSeq)
library(scatterpie)
# library(batchelor)


### 1. Setting working directory and loading data / saving data ----------------------------------------------------------------
setwd("~/")
load('FSI_CD4subset_wTreg.Robj')
load('FSP_CD4subset_wTreg.Robj')
load('FLN_CD4subset_wTreg.Robj')

#### 3. integrate ---------------------------------------------------------------------------------

FLN_CD4 <- SetIdent(FLN_CD4, value = "orig.ident")
FSP_CD4 <- SetIdent(FSP_CD4, value = "orig.ident")
FSI_CD4 <- SetIdent(FSI_CD4, value = "orig.ident")

colnames(FSP_CD4@meta.data) 
colnames(FSP_CD4@meta.data) [27] <- 'TotalTcells_clusterID'
FSP_CD4$TotalTcells_clusterID <- paste0('FSP - ', FSP_CD4$TotalTcells_clusterID)

colnames(FLN_CD4@meta.data)
colnames(FLN_CD4@meta.data) [27] <- 'TotalTcells_clusterID'
FLN_CD4$TotalTcells_clusterID <- paste0('FLN - ', FLN_CD4$TotalTcells_clusterID)

colnames(FSI_CD4@meta.data)
colnames(FSI_CD4@meta.data) [25] <- 'TotalTcells_clusterID'
FSI_CD4$TotalTcells_clusterID <- paste0('FSI - ', FSI_CD4$TotalTcells_clusterID)

### merge to get the layers
tissues.merged <- merge(FSP_CD4, c(FLN_CD4,FSI_CD4), merge.data=T)

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
        group.by = 'harmony_clusters_res.1.0')

#### merge all classical naive clusters --------
DimPlot(CD4_tissues_wTreg_2,reduction = "umap.harmony" ,
        group.by = 'harmony_clusters_res.0.7', pt.size=1)+
  scale_color_manual(values=c(palette.colors(palette = "R4"),
                              'darkgreen','maroon2', 'orange2', 'violet', 'gold','skyblue'))

mapping <- c('0' = '1', '1' = '0', '2' = '0', '3' = '0', '4' = '4',
             '5' = '7', '6' = '2', '7' = '8', '8' = '3', '9' = '5', '10' = '6')
old_clusters <- as.character(CD4_tissues_wTreg_2$harmony_clusters_res.0.7)
new_clusters <- mapping[old_clusters]
names(new_clusters) <- colnames(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2$tissues_res.0.7_manual <- factor(new_clusters)

DimPlot(CD4_tissues_wTreg_2,reduction = "umap.harmony" ,
        group.by = 'tissues_res.0.7_manual', pt.size=1)+
  scale_color_manual(values=c(palette.colors(palette = "R4"),
                              'darkgreen','maroon2', 'orange2', 'violet', 'gold','skyblue'))

mapping_tissue <- c('FLN' = 'MLN', 'FSI' = 'SI', 'FSP' = 'SP')
old_tissues <- as.character(CD4_tissues_wTreg_2$orig.ident)
new_tissues <- mapping_tissue[old_tissues]
names(new_tissues) <- colnames(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2$orig.ident <- factor(new_tissues)

rm(CD4_tissues_wTreg)
rm(tissues.merged)
rm(tissues.integrated)
rm(FSI_CD4)

rm(FLN_CD4)
rm(FSP_CD4)

## save UMAP coordinates
# setwd("")
# umap_coords <- Embeddings(CD4_tissues_wTreg_2, "umap.harmony")
# saveRDS(umap_coords, "umap_coords.rds")

#### 5. Figures Paper -------------------------------------------------------------
setwd("")
cluster_colors <- c('black','coral4', 'lightblue3','maroon','orange','gold', 'seagreen','turquoise','deepskyblue2')
tissue_colors <- c('turquoise', 'gold', 'maroon')
donor_colors <- c(palette.colors(palette = "R4")[1:4],palette.colors(palette = "R4")[6:7])

DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)

CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_clusters <- as.factor(CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_manual)
levels(CD4_tissues_wTreg_2@meta.data$tissues_res.0.7_clusters) <- c( '0: Classical Naive', '1: Tissue Naive-like',
                                                                     '2: Early activated/stressed',
                                                                     '3: Circulating Memory',
                                                                     '4: Th1', '5: Th17', 
                                                                     '6: Activated/regulatory', '7: Treg 1', '8: Treg 2')
#### Fig 1 ####
pdf('CD4_tissues_wTreg_noD010_clustergenes_fig1_res.0.7.pdf', width=18, height=6)
FeaturePlot(SetIdent(CD4_tissues_wTreg_2, value = "tissues_res.0.7_manual"), 
            c('adt_CD45RA','FOXP3', 'KLRB1','HOPX', 'ZBTB16', 'IL2', 
              'adt_CD69.1', 'RORC', 'CCR6','STAT4','CXCR3','IFNG'), 
            label=T, pt.size=.001, label.size = 4, ncol=6, order=T, cols=c('gold2','purple'))&
  ylab('UMAP')&xlab('UMAP')&
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_text(size=16), 
        legend.text = element_text(size=10),legend.key.size = unit(10,'points'),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

fig1a <- DimPlot(object = CD4_tissues_wTreg_2, reduction = "umap.harmony", pt.size=1, group.by = 'orig.ident', 
                 cols=tissue_colors) +
  ylab('UMAP')+xlab('UMAP')+
  theme(axis.text = element_blank(), axis.title = element_text(size=26), axis.ticks = element_blank(),
        legend.text = element_text(size=32),legend.title = element_text(size=30),
        plot.title = element_text(size=34, face='bold', hjust=0.5),plot.subtitle = element_text(size=32, hjust=0.5),
        plot.margin = margin(20,40,20,40),
        text=element_text(size=8))+
  labs(title=expression(bold('Tissue origin of fetal CD4'^+''*' T cells')), color='')+
  guides(color=guide_legend(override.aes=list(size=8)))+
  scale_color_manual(values=tissue_colors, labels=c('MLN', 'Small Intestine', 'Spleen'))

fig1b <- DimPlot(object = CD4_tissues_wTreg_2, reduction = "umap.harmony", pt.size=1, group.by = 'tissues_res.0.7_manual', 
        label = T, label.size = 12,label.color = c('grey80','grey80','grey30','grey30','grey80',
                                                   'grey30','grey30','grey80','grey30'), 
        repel=T,label.box = T,cols=cluster_colors) +
  ylab('UMAP')+xlab('UMAP')+
  theme(axis.text = element_blank(), axis.title = element_text(size=26), axis.ticks = element_blank(),
        legend.text = element_text(size=32),legend.title = element_text(size=30),
        plot.title = element_text(size=34, face='bold', hjust=0.5),plot.subtitle = element_text(size=32, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))+
  labs(title=expression(bold('Cluster ID of fetal CD4'^+''*' T cells')), color='')+
  guides(color=guide_legend(override.aes=list(size=7)))+
  scale_color_manual(values=cluster_colors, 
                     labels= c( '0: Classical Naive', '1: Tissue Naive-like',
                                '2: Early activated/stressed',
                                '3: Circulating Memory',
                                '4: Th1', '5: Th17', 
                                '6: Activated/regulatory', '7: Treg 1', '8: Treg 2'))

percentages_0.5 <- CD4_tissues_wTreg_2@meta.data %>% group_by(orig.ident, tissues_res.0.7_clusters) %>%
  summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ungroup()

percentages_0.5$text_color <- c('grey60',rep('grey20',8),'grey60',rep('grey20',8),'grey60',rep('grey20',7))

fig1c <- ggplot(CD4_tissues_wTreg_2@meta.data, aes(x=orig.ident, fill=tissues_res.0.7_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  geom_text_repel( data = percentages_0.5, force=0.01,
                   min.segment.length = unit(0.2, 'lines'),
                   aes(x = orig.ident, y = freq, 
                       label = scales::percent(freq, accuracy = 0.1), color=text_color, 
                       group = tissues_res.0.7_clusters),
                   position = position_fill(vjust = 0.5), 
                   size = 6) +
  labs(title=expression(bold('Tissue distribution of fetal CD4'^+''*' T cell clusters')), subtitle='')+
  theme(axis.text.x = element_text(size=34), 
        axis.text.y = element_text(size=28),axis.title.y = element_text(size=30,vjust=1.5),
              legend.text = element_text(size=32),legend.title = element_text(size=30),
              plot.title = element_text(size=34, face='bold', hjust=0.5),plot.subtitle = element_text(size=32, hjust=0.5),
              plot.margin = margin(20,150,20,40),
              text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=cluster_colors)+
  scale_color_manual(values=c('black', 'grey60'))+
  guides(color='none', fill='none')

fig1_1 <- fig1a + fig1b + fig1c + plot_layout(ncol = 3)
ggsave("CD4_tissues_wTreg_noD010_clusteringoverview_fig1_res.0.7_NEW.pdf", fig1_1, width = 32, height = 9)

#### Fig 2 ####
##horizontal legend
leg2 <- ggplot(CD4_tissues_wTreg_2@meta.data, aes(x=orig.ident, fill=tissues_res.0.7_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(title=expression(bold('Tissue distribution of fetal CD4'^+''*' T cell clusters')), subtitle='')+
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=22),axis.title.y = element_text(size=24,vjust=1.5),
        legend.text = element_text(size=22),legend.title = element_blank(),
        legend.position = 'bottom', legend.byrow = T,
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,150,20,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=cluster_colors)+
  scale_color_manual(values=c('black', 'grey60'))

fig2a <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped),],
       aes(x=orig.ident, fill=Clone_Frequency_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*alpha*beta*" clone size"),
title=expression(bold("TCR"*alpha*beta*" clonality in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=36), 
        axis.text.y = element_text(size=30),axis.title.y = element_text(size=34,vjust=1.5),
        legend.text = element_text(size=34),legend.title = element_text(size=36),
        plot.title = element_text(size=40, face='bold', hjust=0.5, vjust=2),plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(20,20,20,60),
        legend.key.size = unit(1, "cm"),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito"))[c(1,2,3)])+
  guides(fill = guide_legend(override.aes = list(size = 7)))

fig2b <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped),],
       aes(x=tissues_res.0.7_clusters, fill=Clone_Frequency_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  ylab(expression("% of CD4"^+''*" T cells/tissue cluster")) + 
  labs(fill = expression("TCR"*alpha*beta*" clone size"),
       title=expression(bold("TCR"*alpha*beta*" clonality in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=36/1.037,  color=cluster_colors),
        axis.text.y = element_text(size=30/1.037),axis.title.y = element_text(size=34/1.037,vjust=1.5),
        axis.title.x = element_text(size=36/1.037),
        legend.text = element_text(size=34/1.037),legend.title = element_text(size=36/1.037),
        plot.title = element_text(size=40/1.037, face='bold', hjust=0.5),plot.subtitle = element_text(size=36/1.037, hjust=0.5),
        plot.margin = margin(20,80,20,20),panel.spacing = unit(2, "lines"),
        text=element_text(size=8), strip.text = element_text(size=34/1.037),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(expand = c(0,0), limits=c(0,0.1),labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito"))[c(1,2,3)])+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))

heatmap_data_res.0.7_tissue <- overlapping_clones_res.0.7_tissue %>%
  group_by(tissues_res.0.7_manual_group1, orig.ident_group1, tissues_res.0.7_manual_group2, orig.ident_group2) %>%
  summarise(Mean_Jaccard = mean(JaccardIndex, na.rm = TRUE), summedClones = sum(SharedClones),
            meanClones=mean(SharedClones),.groups = "drop")%>%
  mutate(Mean_Jaccard_color=ifelse(Mean_Jaccard<0.001,T,F), 
         summedClones_color=ifelse(summedClones<3,T,F),
         meanClones_color=ifelse(meanClones<0.3,T,F),
         group1=paste0(tissues_res.0.7_manual_group1, '_',orig.ident_group1),
         group2=paste0(tissues_res.0.7_manual_group2,'_',orig.ident_group2),
         comparison=paste0(pmin(group1,group2), '-', pmax(group1,group2))) %>%
  distinct(comparison, .keep_all = TRUE)
heatmap_data_res.0.7_tissue <- heatmap_data_res.0.7_tissue[heatmap_data_res.0.7_tissue$group1!=heatmap_data_res.0.7_tissue$group2,]

fig2c <- ggplot(heatmap_data_res.0.7_tissue, aes(x = group1, y = group2, fill = Mean_Jaccard)) +
  geom_tile() +
  scale_color_manual(values=c('black','grey40'))+
  guides(color='none')+ scale_fill_gradient(low='#310062', high='gold', limits=c(0,0.025)) + 
  theme_classic()+ 
  labs(title = expression(bold("TCR"*alpha*beta*" clonal overlap")),
       fill = "Mean Jaccard Index")+
  xlab('Cluster ID_Tissue')+ylab('Cluster ID_Tissue')+
  theme(axis.text.x = element_text(color=c(rep(cluster_colors[1:5],each=3), rep(cluster_colors[6],1), rep(cluster_colors[7:9],each=3)),
                                   angle=90, size=36), 
        axis.text.y = element_text(size=32, color=c(rep(cluster_colors[1],each=2), rep(cluster_colors[2:5],each=3),rep(cluster_colors[6],1), rep(cluster_colors[7:9],each=3))), 
        axis.title.y = element_text(size=34*1.2, vjust=3),   axis.title.x = element_text(size=34*1.2, vjust=-1), 
        legend.text = element_text(size=34*1.2),legend.title = element_text(size=36*1.2),
        plot.title = element_text(size=40*1.2, face='bold', hjust=0.5),plot.subtitle = element_text(size=36*1.2, hjust=0.5),
        legend.key.size = unit(1.3, "cm"),plot.margin = margin(5,20,20,20))


fig2d <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped),],
                aes(x=orig.ident, fill=Clone_Frequency_TRA_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*alpha*" clonotype size"), 
       title=expression(bold("TCR"*alpha*" frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=36*1.06), 
        axis.text.y = element_text(size=30*1.06),axis.title.y = element_text(size=34*1.06,vjust=1.5),
        legend.text = element_text(size=34*1.06),legend.title = element_text(size=36*1.06),
        plot.title = element_text(size=40*1.06, face='bold', hjust=0.5, vjust=2),plot.subtitle = element_text(size=36*1.06, hjust=0.5),
        plot.margin = margin(20,20,20,60),
        legend.key.size = unit(1.06, "cm"),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))

fig2e <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped),],
                aes(x=tissues_res.0.7_clusters, fill=Clone_Frequency_TRA_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  ylab(expression("% of CD4"^+''*" T cells/tissue cluster")) + 
  labs(fill = expression("TCR"*alpha*" clonotype size"),
       title=expression(bold("TCR"*alpha*" clonotype frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=36,  color=cluster_colors),
        axis.text.y = element_text(size=30),axis.title.y = element_text(size=34,vjust=1.5),
        axis.title.x = element_text(size=36),
        legend.text = element_text(size=34),legend.title = element_text(size=36),
        plot.title = element_text(size=40, face='bold', hjust=0.5),plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(20,80,20,20),panel.spacing = unit(2, "lines"),
        text=element_text(size=8), strip.text = element_text(size=34),
        legend.key.size = unit(1, "cm"))+
  scale_y_continuous(expand = c(0,0), limits=c(0,0.5),labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))

heatmap_data_res.0.7_tissue_TRA <- overlapping_clones_res.0.7_tissue_TRA %>%
  group_by(tissues_res.0.7_manual_group1, orig.ident_group1, tissues_res.0.7_manual_group2, orig.ident_group2) %>%
  summarise(Mean_Jaccard = mean(JaccardIndex, na.rm = TRUE), summedClones = sum(SharedClones),
            meanClones=mean(SharedClones),.groups = "drop")%>%
  mutate(Mean_Jaccard_color=ifelse(Mean_Jaccard<0.005,T,F), 
         summedClones_color=ifelse(summedClones<8,T,F),
         meanClones_color=ifelse(meanClones<1.3,T,F),
         group1=paste0(tissues_res.0.7_manual_group1, '_',orig.ident_group1),
         group2=paste0(tissues_res.0.7_manual_group2,'_',orig.ident_group2),
         comparison=paste0(pmin(group1,group2), '-', pmax(group1,group2))) %>%
  distinct(comparison, .keep_all = TRUE)
heatmap_data_res.0.7_tissue_TRA <- heatmap_data_res.0.7_tissue_TRA[heatmap_data_res.0.7_tissue_TRA$group1!=heatmap_data_res.0.7_tissue_TRA$group2,]

fig2f <- ggplot(heatmap_data_res.0.7_tissue_TRA, aes(x = group1, y = group2, fill = Mean_Jaccard)) +
  geom_tile() + 
  scale_color_manual(values=c('black','grey40'))+
  guides(color='none')+ scale_fill_gradient(low='#310062', high='gold') + 
  theme_classic()+ 
  labs(title = expression(bold("TCR"*alpha*" clonotype overlap")),
       fill = "Mean Jaccard Index")+
  xlab('Cluster ID_Tissue')+ylab('Cluster ID_Tissue')+
  theme( axis.text.x = element_text(size=36, color=c(rep(cluster_colors[1:5],each=3), rep(cluster_colors[6],1), rep(cluster_colors[7:9],each=3)),
                                   angle=90), 
        axis.text.y = element_text(size=32,color=c(rep(cluster_colors[1],each=2), rep(cluster_colors[2:5],each=3),rep(cluster_colors[6],1), rep(cluster_colors[7:9],each=3))), 
        axis.title.y = element_text(size=34*1.2, vjust=3),   axis.title.x = element_text(size=34*1.2, vjust=-1), 
        legend.text = element_text(size=34*1.2),legend.title = element_text(size=36*1.2),
        plot.title = element_text(size=40*1.2, face='bold', hjust=0.5),plot.subtitle = element_text(size=36*1.2, hjust=0.5),
        legend.key.size = unit(1.3, "cm"),
        plot.margin = margin(5,20,20,20))

CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data[sample(nrow(CD4_tissues_wTreg_2_TCR@meta.data)), ]
fig2h <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Nadditions_TRA)&CD4_tissues_wTreg_2_TCR$Nadditions_TRA<50&
                                                    !is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA),], 
                aes(x=as.factor(Clone_Frequency_TRA), y=Nadditions_TRA)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.5, size=2, aes(color=orig.ident, fill=orig.ident), shape=21,stroke=0.2) +  # Add jitter for visibility
  labs(title=expression(bold("N-additions vs TCR"*alpha*" clonotype size")),
       x=expression("TCR"*alpha*" clonotype size"),
       y=expression("# of N-additions (TCR"*alpha*")"),
       color='Fetal tissue',
       fill='Fetal tissue') +
  scale_color_manual(values=tissue_colors)+
  scale_fill_manual(values=tissue_colors)+
  theme_classic() +
  theme(axis.text = element_text(size=30), axis.title = element_text(size=34),
        axis.text.x.top  = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top  = element_blank(),
        legend.text = element_text(size=34),legend.title = element_text(size=36),
        plot.title = element_text(size=40, face='bold', hjust=0.5),plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(20,20,20,20), 
        text=element_text(size=30), plot.caption.position = 'plot')+
  guides(color = guide_legend(override.aes = list(size = 10)))


fig2i <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Nadditions_TRA),], 
       aes(x=tissues_res.0.7_clusters, y=Nadditions_TRA)) +
  geom_boxplot() + theme_classic() +
  scale_y_sqrt(limits=c(0,35), expand=c(0.01,0), breaks=c(0,1,3,5,10,20,30))+
  labs(title=expression(bold("N-additions in TCR"*alpha*" chain of fetal CD4"^+''*" T cells")),
       x="Cluster ID",
       y=expression("# of N-additions (TCR"*alpha*")"))+
  theme(axis.text.x = element_text(size=36,  color=cluster_colors),
        axis.text.y = element_text(size=28),axis.title.y = element_text(size=34),
        axis.title.x = element_text(size=36),
        legend.text = element_text(size=34),legend.title = element_text(size=36),
        plot.title = element_text(size=40, face='bold', hjust=0.5),plot.subtitle = element_text(size=36, hjust=0.5),
        plot.margin = margin(20,20,20,20),panel.spacing = unit(2, "lines"),
        text=element_text(size=30), strip.text = element_text(size=34))+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))

CD4_tissues_wTreg_2_TCR$Publicity_TRA <- as.numeric(CD4_tissues_wTreg_2_TCR$Publicity_TRA)
fig2g_data <- CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity_TRA),c('Publicity_TRA','Clone_Frequency_TRA', 'orig.ident')]
fig2g_data <- fig2g_data %>% group_by(Publicity_TRA,Clone_Frequency_TRA,orig.ident) %>% summarize(freq=n())
fig2g_data <- pivot_wider(fig2g_data,names_from = orig.ident, values_from = freq)
fig2g_data$size <- rowSums(fig2g_data[,3:5], na.rm=T)
fig2g_data[4,3] <- 0
fig2g_data[4,5] <- 0
fig2g_data$Clone_Frequency_TRA <- as.factor(fig2g_data$Clone_Frequency_TRA)
fig2g <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity_TRA),], 
                aes(x=as.factor(Clone_Frequency_TRA), y=Publicity_TRA)) +
  geom_boxplot(median.linewidth = 2)+
  labs(title=expression(bold("TCR"*alpha*" publicity vs TCR"*alpha*" clonotype size")),
       x=expression("TCR"*alpha*" clonotype size"),
       y=expression("# of donors sharing TCR"*alpha*" clonotype")) +
  theme_classic() +
  theme(axis.text = element_text(size=30*1.07), axis.title = element_text(size=32),
        axis.text.x.top  = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top  = element_blank(),
        legend.text = element_text(size=34*1.07),legend.title = element_text(size=36*1.07),
        plot.title = element_text(size=40*1.07, face='bold', hjust=0.5),plot.subtitle = element_text(size=36*1.07, hjust=0.5),
        plot.margin = margin(20,20,20,20), 
        text=element_text(size=22), plot.caption.position = 'plot')

ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2a_res.0.7_NEW.pdf", fig2a, width=12,height=7.2)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2b_res.0.7_NEW.pdf", fig2b, width=19.2,height=7.2)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2c_res.0.7_NEW.pdf", fig2c, width=19,height=11)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2d_res.0.7_NEW.pdf", fig2d, width=13,height=7.2)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2e_res.0.7_NEW.pdf", fig2e, width=20,height=7.2)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2f_res.0.7_NEW.pdf", fig2f, width=19,height=11)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2g_res.0.7_NEW.pdf", fig2g, width=12,height=7.2)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2h_res.0.7_NEW.pdf", fig2h, width=13,height=7.2)
ggsave("CD4_tissues_wTreg_noD010_TCRrepertoire_fig2i_res.0.7_NEW.pdf", fig2i, width=15.6,height=7.2)
ggsave('horizontallegend.pdf', leg2, width=16,height=13)

#### Fig 3 ####
DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)

pdf('CD4_tissues_wTreg_noD010_naiveclustergenes_fig3a_res.0.7_NEW.pdf', width=32.5, height=10)
VlnPlot(SetIdent(CD4_tissues_wTreg_2, value = "tissues_res.0.7_clusters"), 
            c('CD27','CCR7', 'SELL', 'LEF1', 
              'TCF7', 'KLF2', 'CD44',
              'BACH2','CD38','HES4'), 
          ncol=5, pt.size=0, assay='RNA', slot='data', same.y.lims = T,
        cols=cluster_colors)&
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8')) &
  scale_y_continuous(expand=c(0,0))&
  theme(axis.text = element_text(size=22), axis.title = element_text(size=26), legend.text = element_text(size=22),
        plot.title = element_text(size=28, hjust=0.5), plot.subtitle = element_text(size=26, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=22))
dev.off()

pdf('CD4_tissues_wTreg_noD010_VlnCD45RA_fig3a_res.0.7.pdf', width=7, height=4)
VlnPlot(SetIdent(CD4_tissues_wTreg_2, value = "tissues_res.0.7_clusters"), 
        'CD45RA', pt.size=0, assay='RNA', slot='data', same.y.lims = T,
        cols=cluster_colors[c(1,2,7,6,5,8,9,3,4)])&
  scale_x_discrete(labels=c('0','1','6','5','4','7','8','2','3')) &
  scale_y_continuous(expand=c(0,0))&
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5), plot.subtitle = element_text(size=18, hjust=0.5),
        plot.margin = margin(10,10,10,10),
        text=element_text(size=14))
dev.off()

pdf('CD4_tissues_wTreg_noD010_res.0.7_manual_proportionplot_fig3a.pdf', width=12, height=7)
ggplot(CD4_tissues_wTreg_2@meta.data, aes(x=orig.ident, fill=tissues_res.0.7_clusters)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  geom_text_repel( data = percentages_0.5, force=0.01,
                   min.segment.length = unit(0.2, 'lines'),
                   aes(x = orig.ident, y = freq, 
                       label = scales::percent(freq, accuracy = 0.1), color=text_color, 
                       group = tissues_res.0.7_clusters),
                   position = position_fill(vjust = 0.5), 
                   size = 6) +
  labs(fill = "Cluster ID", title=expression(bold('Tissue distribution of fetal CD4'^+''*' T cell clusters')), subtitle='')+
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=22),axis.title.y = element_text(size=24,vjust=1.5),
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        plot.title = element_text(size=30, face='bold'),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=cluster_colors[c(1,2,7,6,5,8,9,3,4)])+
  scale_color_manual(values=c('black', 'grey60'))+
  guides(color='none')
dev.off()

DEG_C1vsC0 <- read_xlsx('~/PhD/Fetal 10X/DEG/CD4_Harmony_FetalTissues/RNA_CD4_wTreg_Harmony_RNA_MAST_0.7_manual_noD010_Naiveclusters_1vfs0.xlsx')
Fig_Volcano_C1vsC0 <- ggplot(DEG_C1vsC0, 
                             aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(colour = abs(avg_log2FC)), size=4) +
  ggtitle(expression(bold('Differential gene expression between CD45RA'^+''*' naive-like clusters')), 
          subtitle = 'Cluster 0                                                                                          Cluster 1') +
  geom_text_repel(aes(label=gene,x = avg_log2FC, y = -log10(p_val_adj)), 
                  size=7, direction='both', nudge_y = 0.25,
                  max.overlaps = 15)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_gradient(low = "gold", high = "blue") +
  scale_y_continuous(limits=c(0,300),expand=c(0,0))+
  scale_x_continuous(limits=c(-4,4))+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        axis.text = element_text(size=32),
        axis.title.x = element_text(size=36, vjust=-2),axis.title.y = element_text(size=36, vjust=3),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=46, face='bold', hjust=0.5, vjust=2),
        plot.subtitle = element_text(size=48, hjust=0.5, vjust=1),
        plot.margin = margin(40,40,40,40))
Fig_Volcano_C1vsC0
ggsave(Fig_Volcano_C1vsC0,filename=('C1vsC0_naiveclusters_res.0.7_manual_CD45RAhighsubsets_Volcano_fig3c.pdf'), height=14, width=26)

#### Fig S2 ####
donor_colors_1 <- c(palette.colors(palette = "R4")[1:4],palette.colors(palette = "R4")[6:8])
cluster_colors_1 <-  c('black','coral4', 'turquoise','blue3','grey','orange','maroon','lightblue3','gold', 'seagreen')
cluster_colors_2 <- c('coral4','black','blue3','grey50', 'orange', 'turquoise','lightblue3','deepskyblue2','maroon','gold', 'seagreen')

CD4_tissues_wTreg <- FindVariableFeatures(CD4_tissues_wTreg)
VariableFeatures(CD4_tissues_wTreg) <- VariableFeatures(CD4_tissues_wTreg)[-grep('TRBV|TRGV|TRAV|TRDV|TRBC|TRGC|TRAC|TRDC|TRBJ|TRGJ|TRAJ|TRDC',VariableFeatures(CD4_tissues_wTreg))]
CD4_tissues_wTreg <- ScaleData(CD4_tissues_wTreg)
CD4_tissues_wTreg <- RunPCA(CD4_tissues_wTreg)
CD4_tissues_wTreg <- RunUMAP(CD4_tissues_wTreg, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
CD4_tissues_wTreg <- FindNeighbors(object = CD4_tissues_wTreg, reduction = "harmony", dims = 1:30)
CD4_tissues_wTreg <- FindClusters(CD4_tissues_wTreg, resolution = seq(0.1:1, by=0.1),
                                    cluster.name = c('harmony_clusters_res.0.1', 'harmony_clusters_res.0.2',
                                                     'harmony_clusters_res.0.3','harmony_clusters_res.0.4','harmony_clusters_res.0.5',
                                                     'harmony_clusters_res.0.6','harmony_clusters_res.0.7','harmony_clusters_res.0.8',
                                                     'harmony_clusters_res.0.9','harmony_clusters_res.1.0'))

FigS2a <- DimPlot(object = CD4_tissues_wTreg, reduction = "umap.harmony", pt.size=1, group.by = 'harmony_clusters_res.0.5', 
                 label = T, label.size = 12, label.color = 'grey30',
                 repel=T,label.box = T,cols=cluster_colors_1) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))+
  labs(title=expression(bold('Cluster ID - before selection')), color='')+
  guides(color=guide_legend(override.aes=list(size=7)))

FigS2b <- FeaturePlot(CD4_tissues_wTreg, 'HSPA6', pt.size=1) + 
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        legend.key.size = unit(24, 'pt'),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))

FigS2c <- DimPlot(object = CD4_tissues_wTreg, reduction = "umap.harmony", pt.size=1, group.by = 'Donor', 
                  cols=donor_colors_1) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))+
  labs(title=expression(bold('Donor ID - before selection')), color='')+
  guides(color=guide_legend(override.aes=list(size=7)))

FigS2d <- DimPlot(object = CD4_tissues_wTreg_2, reduction = "umap.harmony", pt.size=1, group.by = 'harmony_clusters_res.0.7', 
                  label = T, label.size = 12, label.color = 'grey30',
                  repel=T,label.box = T,cols=cluster_colors_2) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))+
  labs(title=expression(bold('Cluster ID - after selection')), color='')+
  guides(color=guide_legend(override.aes=list(size=7)))

FigS2e <- DimPlot(object = CD4_tissues_wTreg_2, reduction = "umap.harmony", pt.size=1, group.by = 'Donor', 
                  cols=donor_colors) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))+
  labs(title=expression(bold('Donor ID - after selection')), color='')+
  guides(color=guide_legend(override.aes=list(size=7)))

DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)

Figs2_2a <- FeaturePlot(CD4_tissues_wTreg_2, 'nCount_RNA', pt.size=1) + 
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        legend.key.size = unit(24, 'pt'),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))
Figs2_2b <- FeaturePlot(CD4_tissues_wTreg_2, 'nFeature_RNA', pt.size=1) + 
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        legend.key.size = unit(24, 'pt'),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))
Figs2_2c <- FeaturePlot(CD4_tissues_wTreg_2, 'percent.mt', pt.size=1) + 
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        legend.key.size = unit(24, 'pt'),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))
Figs2_2d <- FeaturePlot(CD4_tissues_wTreg_2, 'percent.ribo', pt.size=1) + 
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        legend.key.size = unit(24, 'pt'),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))
Figs2_2e <-   DimPlot(CD4_tissues_wTreg_2, group.by = 'Phase', reduction = 'umap.harmony', pt.size=1,
                      cols=c('grey20', 'violet', 'seagreen1'))+
        guides(color=guide_legend(override.aes=list(size=7)))+ 
          theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        legend.key.size = unit(24, 'pt'),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,40,20,20),
        text=element_text(size=8))
                
FigS2_3 <- VlnPlot(SetIdent(CD4_tissues_wTreg_2, value = "tissues_res.0.7_manual"), 
                   features=CD4_tissues_wTreg_2@assays[["ADT"]]@counts@Dimnames[[1]], 
                   ncol=5, pt.size=0, assay='RNA', slot='data', same.y.lims = T,
                   cols=cluster_colors[c(1,2,7,6,5,8,9,3,4)])&
  scale_x_discrete(labels=c('0','1','6','5','4','7','8','2','3')) &
  scale_y_continuous(expand=c(0,0))&
  theme(axis.text = element_text(size=26), axis.title = element_text(size=28), 
        legend.text = element_text(size=26),
        plot.title = element_text(size=32, hjust=0.5), 
        plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,20),
        text=element_text(size=20))

FigS2 <- plot_grid(FigS2a, FigS2b,FigS2c,FigS2d,FigS2e, 
                   Figs2_2a, Figs2_2b,Figs2_2c,Figs2_2d,Figs2_2e, nrow=2,
                   align='h', axis='l')

ggsave(plot=FigS2,filename=('CD4_tissues_wTreg_clusteringoverview_figS2_QC_Donorselection.pdf'), height=14, width=42)
ggsave(plot=FigS2_3,filename=('CD4_tissues_wTreg_clusteringoverview_figS2_surfacemarkers.pdf'), height=14, width=42)

## use script from ES_fetal10X_blood.R to make figure S3A
## use script from ES_fetal10X_tissuescombined_CD4_wTreg_Harmony_wFBL.R to make the rest of Figure S3

#### Fig S5 ####
FigS5a <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data,
                          aes(x=orig.ident, fill=Clone_Frequency_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*alpha*beta*" clone size"),
       title=expression(bold("TCR"*alpha*beta*" clonality in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito"))[c(1:3)], na.value='black')

FigS5b <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped),],
                 aes(x=orig.ident, fill=Clone_Frequency_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*alpha*beta*" clone size"),
       title=expression(bold("TCR"*alpha*beta*" clonality in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(limits=c(0,0.05), expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito"))[c(1,2,3)])

FigS5c <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data,
                 aes(x=orig.ident, fill=Clone_Frequency_TRA_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*alpha*" clonotype size"),
       title=expression(bold("TCR"*alpha*" clonotype frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")), na.value='black')

FigS5d <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped),],
                 aes(x=orig.ident, fill=Clone_Frequency_TRA_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*alpha*" clonotype size"),
       title=expression(bold("TCR"*alpha*" clonotype frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(limits=c(0,0.31), expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))

FigS5e <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data,
                 aes(x=orig.ident, fill=Clone_Frequency_TRB_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*beta*" clonotype size"),
       title=expression(bold("TCR"*beta*" clonotype frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")), na.value='black')

FigS5f <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped),],
                 aes(x=orig.ident, fill=Clone_Frequency_TRB_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("") + ylab(expression("% of CD4"^+''*" T cells/tissue")) + 
  labs(fill = expression("TCR"*beta*" clonotype size"),
       title=expression(bold("TCR"*beta*" clonotype frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(limits=c(0,0.31), expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))

FigS5g <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped),],
                 aes(x=tissues_res.0.7_clusters, fill=Clone_Frequency_TRB_grouped)) + theme_classic() +
  geom_bar(position = "fill") + xlab("Cluster ID") + ylab("Fraction") + 
  ylab(expression("% of CD4"^+''*" T cells/tissue cluster")) + 
  labs(fill = expression("TCR"*beta*" clonotype size"),
       title=expression(bold("TCR"*beta*" clonotype frequency in fetal CD4"^+''*" T cells")))+
  theme(axis.text.x = element_text(size=22,  color=cluster_colors),
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        axis.title.x = element_text(size=22),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=8), strip.text = element_text(size=28))+
  scale_y_continuous(expand = c(0,0), limits=c(0,0.5),labels=scales::percent_format())+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))

heatmap_data_res.0.7_tissue_TRB <- overlapping_clones_res.0.7_tissue_TRB %>%
  group_by(tissues_res.0.7_manual_group1, orig.ident_group1, tissues_res.0.7_manual_group2, orig.ident_group2) %>%
  summarise(Mean_Jaccard = mean(JaccardIndex, na.rm = TRUE), summedClones = sum(SharedClones),
            meanClones=mean(SharedClones),.groups = "drop")%>%
  mutate(Mean_Jaccard_color=ifelse(Mean_Jaccard<0.004,T,F), 
         summedClones_color=ifelse(summedClones<8,T,F),
         meanClones_color=ifelse(meanClones<2,T,F),
         group1=paste0(tissues_res.0.7_manual_group1, '_',orig.ident_group1),
         group2=paste0(tissues_res.0.7_manual_group2,'_',orig.ident_group2),
         comparison=paste0(pmin(group1,group2), '-', pmax(group1,group2))) %>%
  distinct(comparison, .keep_all = TRUE)
heatmap_data_res.0.7_tissue_TRB <- heatmap_data_res.0.7_tissue_TRB[heatmap_data_res.0.7_tissue_TRB$group1!=heatmap_data_res.0.7_tissue_TRB$group2,]

FigS5h <- ggplot(heatmap_data_res.0.7_tissue_TRB, aes(x = group1, y = group2, fill = Mean_Jaccard)) +
  geom_tile() +
  scale_color_manual(values=c('black','grey40'))+
  guides(color='none')+ scale_fill_gradient(low='#310062', high='gold') + 
  theme_classic()+ 
  labs(title = expression(bold("TCR"*beta*" clonotype overlap")),
       fill = "Mean Jaccard Index")+
  xlab('Cluster ID_Tissue')+ylab('Cluster ID_Tissue')+
  theme(axis.text = element_text(size=16), 
        axis.text.x = element_text(color=c(rep(cluster_colors[1:5],each=3), rep(cluster_colors[6],2), rep(cluster_colors[7:9],each=3)),
                                   angle=90), 
        axis.text.y = element_text(color=c(rep(cluster_colors[1],each=2), rep(cluster_colors[2:5],each=3),rep(cluster_colors[6],2), rep(cluster_colors[7:9],each=3))), 
        axis.title.y = element_text(size=22, vjust=3),   axis.title.x = element_text(size=22, vjust=-1), 
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),
        plot.margin = margin(40,40,40,40))

#### per TRB clone, calculate % no TRA, % non-matching TRA, % matching TRA (full clone)
TRA_nonmatchingTRB <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID_TRA)) %>%
  group_by(Donor, CloneID_TRA) %>%
  summarise(
    total_cells = n(),
    n_no_TRB = sum(is.na(CloneID_TRB)),
    n_matching_TRB = sum(Clone_Frequency>1),  # full clone match
    n_nonmatching_TRB = sum(!is.na(CloneID_TRB) & Clone_Frequency==1),
    pct_nonmatching_TRB = 100 * n_nonmatching_TRB / total_cells,
    pct_matching_TRB = 100 * n_matching_TRB / total_cells,
    pct_missing_TRB = 100 * n_no_TRB / total_cells
  ) %>%
  arrange(desc(pct_nonmatching_TRB))
TRA_nonmatchingTRB[is.na(TRA_nonmatchingTRB$pct_matching_TRB), 'pct_matching_TRB'] <- 0
TRA_nonmatchingTRB_melt <- melt(TRA_nonmatchingTRB, 
                                measure.vars=c('pct_nonmatching_TRB', 'pct_missing_TRB', 'pct_matching_TRB'))

TRA_nonmatchingTRB %>% group_by(total_cells) %>% summarise(mean(pct_nonmatching_TRB))
TRA_nonmatchingTRB %>% ungroup() %>% summarise(mean(pct_nonmatching_TRB))
TRA_nonmatchingTRB %>% ungroup() %>% summarise(mean(pct_missing_TRB))
TRA_nonmatchingTRB %>% ungroup() %>% summarise(100*(sum(pct_missing_TRB==100)/n()))

FigS5i <- ggplot(TRA_nonmatchingTRB_melt, aes(x=as.factor(total_cells),y=value,fill=variable)) + 
  geom_bar(stat='identity', position='fill')+
  xlab(expression("TCR"*alpha*" clonotype size")) + ylab(expression("Average %")) + 
  labs(fill = expression("barcode-matched TCR"*beta*" chain"),
       title=expression(bold("TCR"*beta*" information per TCR"*alpha*" clonotype")))+
  theme_classic() + theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=22,vjust=-1.5),
        axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
        legend.text = element_text(size=18),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=c('gold','black','pink'),
                    labels=c('Non-matching', 'Missing', 'Matching'))

#### same for TRB 
TRB_nonmatchingTRA <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID_TRB)) %>%
  group_by(Donor, CloneID_TRB) %>%
  summarise(
    total_cells = n(),
    n_no_TRA = sum(is.na(CloneID_TRA)),
    n_matching_TRA = sum(Clone_Frequency>1),  # full clone match
    n_nonmatching_TRA = sum(!is.na(CloneID_TRA) & Clone_Frequency==1),
    pct_nonmatching_TRA = 100 * n_nonmatching_TRA / total_cells,
    pct_matching_TRA = 100 * n_matching_TRA / total_cells,
    pct_missing_TRA = 100 * n_no_TRA / total_cells
  ) %>%
  arrange(desc(pct_nonmatching_TRA))
TRB_nonmatchingTRA[is.na(TRB_nonmatchingTRA$pct_matching_TRA), 'pct_matching_TRA'] <- 0
TRB_nonmatchingTRA_melt <- melt(TRB_nonmatchingTRA, 
                                measure.vars=c('pct_nonmatching_TRA', 'pct_missing_TRA', 'pct_matching_TRA'))

TRB_nonmatchingTRA %>% group_by(total_cells) %>% summarise(mean(pct_nonmatching_TRA))
TRB_nonmatchingTRA %>% ungroup() %>% summarise(mean(pct_nonmatching_TRA))
TRB_nonmatchingTRA %>% ungroup() %>% summarise(mean(pct_missing_TRA))
TRB_nonmatchingTRA %>% ungroup() %>% summarise(100*(sum(pct_missing_TRA==100)/n()))

FigS5j <- ggplot(TRB_nonmatchingTRA_melt, aes(x=total_cells,y=value,fill=variable)) + 
  geom_bar(stat='identity', position='fill')+
  xlab(expression("TCR"*beta*" clonotype size")) + ylab(expression("Average %")) + 
  labs(fill = expression("barcode-matched TCR"*alpha*" chain"),
       title=expression(bold("TCR"*alpha*" information per TCR"*beta*" clonotype")))+
  theme_classic() + theme(axis.text.x = element_text(size=28), axis.title.x = element_text(size=22,vjust=-1.5),
                          axis.text.y = element_text(size=20),axis.title.y = element_text(size=22,vjust=1.5),
                          legend.text = element_text(size=18),legend.title = element_text(size=26),
                          plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
                          plot.margin = margin(40,40,40,40),
                          text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_x_continuous(breaks=(seq(1,10,1)))+
  scale_fill_manual(values=c('gold','black','pink'),
                    labels=c('Non-matching', 'Missing', 'Matching'))

FigS5 <- plot_grid(FigS5a,FigS5b,FigS5c, FigS5d, FigS5e, FigS5f, FigS5g, FigS5h, FigS5i,
                   FigS5j)

ggsave(plot=FigS5,filename=('CD4_tissues_wTreg_TCR_figS5_NEW.pdf'), height=20, width=50, limitsize=F)

#### Fig S6 ####
CD4_tissues_wTreg_2_TCR$Publicity_TRA <- as.factor(CD4_tissues_wTreg_2_TCR$Publicity_TRA)
FigS6a <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity_TRA),], 
                 aes(x=orig.ident, fill=Publicity_TRA)) +
  geom_bar(position = "fill") + xlab("Tissue") + ylab(expression("% of CD4"^+''*" T cells/tissue")) +
  labs(fill = expression("# of donors sharing TCR"*alpha*" clonotype"), 
       title=expression(bold('TCR'*alpha*' Publicity within fetal tissues')), subtitle='')+
  theme_classic()+
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),axis.title.y = element_text(size=26),
        axis.title.x = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=22), strip.text = element_text(size=28))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=c('black',palette.colors(palette='Dark2')))

CD4_tissues_wTreg_2_TCR$Publicity_TRB <- as.factor(CD4_tissues_wTreg_2_TCR$Publicity_TRB)
FigS6b <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity_TRB),], 
                 aes(x=orig.ident, fill=Publicity_TRB)) +
  geom_bar(position = "fill") + xlab("Tissue") + ylab(expression("% of CD4"^+''*" T cells/tissue")) +
  labs(fill = expression("# of donors sharing TCR"*beta*" clonotype"), 
       title=expression(bold('TCR'*beta*' Publicity within fetal tissues')), subtitle='')+
  theme_classic()+
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),axis.title.y = element_text(size=26),
        axis.title.x = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=22), strip.text = element_text(size=28))+
  scale_y_continuous(expand = c(0,0), labels=scales::percent_format())+
  scale_fill_manual(values=c('black',palette.colors(palette='Dark2')))

CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data[sample(nrow(CD4_tissues_wTreg_2_TCR@meta.data)), ]
FigS6c <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity),], 
                 aes(x=Clone_Frequency, y=Publicity)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.5, size=4.5, aes(color=orig.ident, fill=orig.ident), shape=21, stroke=0.3) +  # Add jitter for visibility
  labs(title=expression(bold("TCR"*alpha*beta*" publicity vs TCR"*alpha*beta*" clonality")),
       x=expression("Clone size (TCR"*alpha*beta*")"),
       y=expression("# of donors sharing TCR"*alpha*beta*" clone"),
       color='Fetal tissue',
       fill='Fetal tissue') +
  scale_color_manual(values=tissue_colors)+
  scale_fill_manual(values=tissue_colors)+
  scale_x_continuous(breaks=seq(0,10,1))+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme_classic() +
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,120,40,80), 
        text=element_text(size=22), plot.caption.position = 'plot')+
  guides(color = guide_legend(override.aes = list(size = 5)))

CD4_tissues_wTreg_2_TCR$Publicity_TRB <- as.numeric(CD4_tissues_wTreg_2_TCR$Publicity_TRB)
FigS6d_data <- CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity_TRB),c('Publicity_TRB','Clone_Frequency_TRB', 'orig.ident')]
FigS6d_data <- FigS6d_data %>% group_by(Publicity_TRB,Clone_Frequency_TRB,orig.ident) %>% summarize(freq=n())
FigS6d_data <- pivot_wider(FigS6d_data,names_from = orig.ident, values_from = freq)
FigS6d_data$size <- rowSums(FigS6d_data[,3:5], na.rm=T)

FigS6d_data[4,3] <- 0
FigS6d_data[4,5] <- 0
FigS6d_data$Clone_Frequency_TRB <- as.factor(FigS6d_data$Clone_Frequency_TRB)
FigS6d <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Publicity_TRB),], 
                aes(x=as.factor(Clone_Frequency_TRB), y=Publicity_TRB)) +
  geom_boxplot(median.linewidth = 2)+
  labs(title=expression(bold("TCR"*beta*" publicity vs TCR"*beta*" clonotype size")),
       x=expression("TCR"*beta*" clonotype size"),
       y=expression("# of donors sharing TCR"*beta*" clonotype")) +
  theme_classic() +
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        axis.text.x.top  = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top  = element_blank(),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,120,40,80), 
        text=element_text(size=22), plot.caption.position = 'plot')

FigS6 <- plot_grid(FigS6a, FigS6b,FigS6c,FigS6d, nrow=2,
                   align='h', axis='l')

ggsave(plot=FigS6,filename=('CD4_tissues_wTreg_TCRpublicity_figS6_NEW.pdf'), height=18, width=30)

#### Fig S7 ####
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data[sample(nrow(CD4_tissues_wTreg_2_TCR@meta.data)), ]
FigS7a <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Nadditions_TRB)&CD4_tissues_wTreg_2_TCR$Nadditions_TRB<50,], 
                 aes(x=Clone_Frequency_TRB, y=Nadditions_TRB)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.5, size=3, aes(color=orig.ident)) +  # Add jitter for visibility
  labs(title=expression(bold("N-additions vs TCR"*beta*" clonotype frequency")),
       x=expression("Clonotype frequency (TCR"*beta*")"),
       y=expression("# of N-additions (TCR"*beta*")"),
       color='Fetal tissue') +
  scale_color_manual(values=tissue_colors)+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9))+
  theme_classic() +
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40), 
        text=element_text(size=22), plot.caption.position = 'plot')+
  guides(color = guide_legend(override.aes = list(size = 5)))

FigS7b <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Nadditions_TRB),], 
                 aes(x=tissues_res.0.7_manual, y=Nadditions_TRB)) +
  geom_boxplot() +  
  theme_classic() +
  scale_y_sqrt(limits=c(0,35), expand=c(0.01,0), breaks=c(0,1,3,5,10,20,30))+
  labs(title=expression(bold("N-additions in TCR"*beta*" chain of fetal CD4"^+''*" T cells")),
       x="Cluster ID",
       y=expression("# of N-additions (TCR"*beta*")"))+
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        axis.text.x = element_text(color=cluster_colors),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=22), strip.text = element_text(size=28))+
  scale_fill_manual(values=rev(palette.colors(palette = "Okabe-Ito")))+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))
FigS7b

F7c_df <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(Nadditions_TRA)) %>%
  # keep ONLY TRUE cases
  mutate(Nadditions_TRA_zero = Nadditions_TRA == 0)

# summarise proportions + CI
plot_data_F7c <- F7c_df %>%
  group_by(orig.ident, tissues_res.0.7_clusters) %>%
  summarise(
    n = n(),
    true_count = sum(Nadditions_TRA_zero),
    prop_true = true_count / n,
    se = sqrt(prop_true * (1 - prop_true) / n),    # standard error
    lower = prop_true - 1.96 * se,                 # 95% CI
    upper = prop_true + 1.96 * se)%>%
  filter(n>10)

FigS7c <-ggplot(plot_data_F7c,
                aes(x = tissues_res.0.7_clusters, y = prop_true)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Cluster ID",
       y = expression("% of CD4"^+''*" T cells/tissue cluster"),
       title = expression(bold("Fraction of fetal TCR"*alpha*" chains with 0 N-additions"))) +
  theme_classic()+
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        axis.text.x = element_text(color=cluster_colors),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=8), strip.text = element_text(size=28))+
  scale_y_continuous(expand = c(0,0), limits=c(0,1),labels=scales::percent_format(), breaks=c(0,0.25,0.5,0.75,1))+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))
FigS7c

F7d_df <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(Nadditions_TRB)) %>%
  # keep ONLY TRUE cases
  mutate(Nadditions_TRB_zero = Nadditions_TRB == 0)

# summarise proportions + CI
plot_data_F7d <- F7d_df %>%
  group_by(orig.ident, tissues_res.0.7_clusters) %>%
  summarise(
    n = n(),
    true_count = sum(Nadditions_TRB_zero),
    prop_true = true_count / n,
    se = sqrt(prop_true * (1 - prop_true) / n),    # standard error
    lower = prop_true - 1.96 * se,                 # 95% CI
    upper = prop_true + 1.96 * se)%>%
  filter(n>10)

FigS7d <-ggplot(plot_data_F7d,
                aes(x = tissues_res.0.7_clusters, y = prop_true)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Cluster ID",
       y = expression("% of CD4"^+''*" T cells/tissue cluster"),
       title = expression(bold("Fraction of fetal TCR"*beta*" chains with 0 N-additions"))) +
  theme_classic()+
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        axis.text.x = element_text(color=cluster_colors),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=8), strip.text = element_text(size=28))+
  scale_y_continuous(expand = c(0,0), limits=c(0,1),labels=scales::percent_format(), breaks=c(0,0.25,0.5,0.75,1))+
  facet_wrap(facets=vars(orig.ident))+
  scale_x_discrete(labels=c('0','1','2','3','4','5','6','7','8'))
FigS7d

FigS7e <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Nadditions_TRA),], 
                 aes(x=orig.ident, y=Nadditions_TRA)) +
  geom_boxplot() +  
  theme_classic() +
  scale_y_sqrt(limits=c(0,35), expand=c(0.01,0), breaks=c(0,1,3,5,10,20,30))+
  labs(title=expression(bold("N-additions (TCR"*alpha*") within fetal tissues")),
       x="Tissue",
       y=expression("# of N-additions (TCR"*alpha*")"))+
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=22), strip.text = element_text(size=28))
FigS7e

FigS7f <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data[!is.na(CD4_tissues_wTreg_2_TCR$Nadditions_TRB),], 
                 aes(x=orig.ident, y=Nadditions_TRB)) +
  geom_boxplot() +  
  theme_classic() +
  scale_y_sqrt(limits=c(0,35), expand=c(0.01,0), breaks=c(0,1,3,5,10,20,30))+
  labs(title=expression(bold("N-additions (TCR"*beta*") within fetal tissues")),
       x="Tissue",
       y=expression("# of N-additions (TCR"*beta*")"))+
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=22), strip.text = element_text(size=28))
FigS7f



F7g_df <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(Nadditions_TRB)) %>%
  group_by(orig.ident, Nadditions_TRB) %>%  summarise(n = n(), Dfreq = mean(is.na(`D-gene TRB`)))

FigS7g <- ggplot(CD4_tissues_wTreg_2_TCR@meta.data%>%
                   filter(!is.na(Nadditions_TRB)), aes(x=is.na(`D-gene TRB`), y=Nadditions_TRB)) +
  geom_boxplot() +
  theme_classic() +
  scale_x_discrete(labels=c('Yes','No'))+
  scale_y_continuous(expand=c(0,0), breaks=c(0,1,3,5,10,20,30), limits=c(0,30))+
  labs(title=expression(bold("N-additions in TCR"*beta*" with/without D-gene")),
       x=expression("D-gene detected"),
       y=expression("# of N-additions (TCR"*beta*")"))+
  theme(axis.text = element_text(size=24), axis.title = element_text(size=26),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=28, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(40,40,40,40),panel.spacing = unit(2, "lines"),
        text=element_text(size=22), strip.text = element_text(size=28))
FigS7g




FigS7 <- plot_grid(FigS7a,FigS7b,FigS7c, FigS7d, FigS7e, FigS7f)
ggsave(plot=FigS7,filename=('CD4_tissues_wTreg_TCR_FigS7_NEW.pdf'), height=18, width=40, limitsize=F)

 
#### Fig S8 ####
FigS8a <- ggplot(subset(CD4_tissues_wTreg_2_TCR@meta.data[CD4_tissues_wTreg_2_TCR$`J-gene TRA`=='TRAJ18'&
                                                              CD4_tissues_wTreg_2_TCR$`V-gene TRA`=='TRAV10',],
                          !is.na(orig.ident)),
                   aes(x=tissues_res.0.7_manual, fill=orig.ident))+
  geom_bar(position = "stack") + xlab("") + ylab(expression("# of cells with TRAJ18-TRAV10 TCR"*alpha)) +
  labs(fill = "Fetal tissue", title=expression(bold('TRAJ18-TRAV10 "NKT cells"')), subtitle='')+
  theme_classic()+ theme(axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=24),axis.title.y = element_text(size=28,vjust=1.5),
        legend.text = element_text(size=26),legend.title = element_text(size=28),
        plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
        plot.margin = margin(20,20,20,40),
        text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), limits=c(0,30))+
  scale_fill_manual(values=tissue_colors)

FigS8b <- ggplot(subset(CD4_tissues_wTreg_2_TCR@meta.data[CD4_tissues_wTreg_2_TCR$`J-gene TRA`=='TRAJ33'&
                                                            CD4_tissues_wTreg_2_TCR$`V-gene TRA`=='TRAV1-2',],
                        !is.na(orig.ident)),
                 aes(x=tissues_res.0.7_manual, fill=orig.ident))+
  geom_bar(position = "stack") + xlab("") + ylab(expression("# of cells with TRAJ33-TRAV1-2 TCR"*alpha)) +
  labs(fill = "Fetal Tissue", title=expression(bold('TRAJ33-TRAV1-2 "MAIT cells"')), subtitle='')+
  theme_classic()+ theme(axis.text.x = element_text(size=28),
                         axis.text.y = element_text(size=24),axis.title.y = element_text(size=28,vjust=1.5),
                         legend.text = element_text(size=26),legend.title = element_text(size=28),
                         plot.title = element_text(size=30, face='bold', hjust=0.5),plot.subtitle = element_text(size=28, hjust=0.5),
                         plot.margin = margin(20,20,20,40),
                         text=element_text(size=8))+
  scale_y_continuous(expand = c(0,0), limits=c(0,5), breaks=c(0,1,2,3,4,5))+
  scale_fill_manual(values=tissue_colors)

FigS8 <- FigS8a + FigS8b + plot_layout(ncol=2)

ggsave(plot=FigS8,filename=('CD4_tissues_wTreg_NKTMAIT_figS8_NEW.pdf'), height=10, width=24, limitsize = F)

#### Fig S9 ####

DEG_MLN_C1vsC0 <- read_xlsx('~/PhD/Fetal 10X/DEG/RNA_CD4_wTreg_Harmony_RNA_MAST_0.7_manual_noD010_MLNonly_Naiveclusters_1vfs0.xlsx')
Fig_Volcano_C1vsC0_MLN <- ggplot(DEG_MLN_C1vsC0, 
                             aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(colour = abs(avg_log2FC)), size=4) +
  ggtitle(expression(bold('Differential gene expression between CD45RA'^+''*' naive-like clusters (MLN only)')), 
          subtitle = 'Cluster 0                                                                                          Cluster 1') +
  geom_text_repel(aes(label=gene,x = avg_log2FC, y = -log10(p_val_adj)), 
                  size=7, direction='both', nudge_y = 0.25,
                  max.overlaps = 15)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_gradient(low = "gold", high = "blue") +
  scale_y_continuous(limits=c(0,110),expand=c(0,0))+
  scale_x_continuous(limits=c(-4,4))+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        axis.text = element_text(size=32),
        axis.title.x = element_text(size=36, vjust=-2),axis.title.y = element_text(size=36, vjust=3),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=46, face='bold', hjust=0.5, vjust=2),
        plot.subtitle = element_text(size=48, hjust=0.5, vjust=1),
        plot.margin = margin(40,40,40,40))
Fig_Volcano_C1vsC0_MLN
ggsave(Fig_Volcano_C1vsC0_MLN,filename=('C1vsC0_naiveclusters_MLNonly_res.0.7_manual_CD45RAhighsubsets_Volcano_figS9a.pdf'), height=14, width=26)

DEG_FLNFSI_C1 <- read_xlsx('~/PhD/Fetal 10X/DEG/DESeq2results_pseudobulk_naive_C1_res.0.7_manual_FSIvsFLN_paired_nod006_CD45RAcutoff.xlsx')
Fig_Volcano_FLNvsFSI <- ggplot(DEG_FLNFSI_C1, aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point(aes(colour = abs(log2FoldChange)), size=4) +
  ggtitle(expression(bold('Pseudobulk differential gene expression within Naive cluster 1')), 
          subtitle = 'MLN                                                                                        Small intestine') +
  geom_text_repel(aes(label=Gene,x = log2FoldChange, y = -log10(padj)), 
                  size=7, direction='both', nudge_y = 0.25,
                  max.overlaps = 15)+
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_gradient(low = "gold", high = "blue") +
  scale_y_continuous(limits=c(0,25),expand=c(0,0))+
  scale_x_continuous(limits=c(-3,3))+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5*.5),
        axis.text = element_text(size=32),
        axis.title.x = element_text(size=36, vjust=-2),axis.title.y = element_text(size=36, vjust=3),
        legend.text = element_text(size=24),legend.title = element_text(size=26),
        plot.title = element_text(size=46, face='bold', hjust=0.5, vjust=2),
        plot.subtitle = element_text(size=48, hjust=0.5, vjust=1),
        plot.margin = margin(40,40,40,40))
Fig_Volcano_FLNvsFSI
ggsave(Fig_Volcano_FLNvsFSI,filename=('C1_FLNvsFSI_naiveclusters_res.0.7_manual_CD45RAcutoff_Volcano_figS9b.pdf'), height=14, width=26)
 

#### 6. DEG ####
#### ROC/MAST for cluster defining markers 
DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2 <- SetIdent(CD4_tissues_wTreg_2, value = "tissues_res.0.7_manual")

##Find Markers that are specific for each cluster
Batchedmarkers.mast_i_RNA_data_0.7=FindAllMarkers(CD4_tissues_wTreg_2, test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                                  min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_0.7$pct.fold <- Batchedmarkers.mast_i_RNA_data_0.7$pct.1/Batchedmarkers.mast_i_RNA_data_0.7$pct.2

## Create list
listDEgenes_i_RNA_MAST_0.7<- split(Batchedmarkers.mast_i_RNA_data_0.7, f=Batchedmarkers.mast_i_RNA_data_0.7$cluster)

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
listDEgenes_i_RNA_MAST_0.7 <-lapply(listDEgenes_i_RNA_MAST_0.7, function(x){dplyr::filter(x, p_val_adj<0.05)})
## Sort on logFC
listDEgenes_i_RNA_MAST_0.7<-lapply(listDEgenes_i_RNA_MAST_0.7,function(x){x<-x[order(x$avg_log2FC, decreasing=T),]})


##save as Robj
setwd("")
save(listDEgenes_i_RNA_MAST_0.7, file='SupplData5_listDEG_fetal_CD4wTreg_scRNA_0.7.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(listDEgenes_i_RNA_MAST_0.7, file='SupplData5_listDEG_fetal_CD4wTreg_scRNA_0.7.xlsx')
detach("package:openxlsx", unload=TRUE)

##Find Markers between the two naive clusters
DefaultAssay(object = CD4_tissues_wTreg_2) <- "RNA"
CD4_tissues_wTreg_2 <- JoinLayers(CD4_tissues_wTreg_2)
CD4_tissues_wTreg_2_Naive <- subset(CD4_tissues_wTreg_2, subset=adt_CD45RA>1)
CD4_tissues_wTreg_2_Naive <- SetIdent(CD4_tissues_wTreg_2_Naive, value = "tissues_res.0.7_manual")

Batchedmarkers.mast_i_RNA_data_Naive =FindMarkers(CD4_tissues_wTreg_2_Naive, ident.1='1', ident.2='0',
                                               test.use = "MAST", slot='data',logfc.threshold = 0.1,
                                               min.cells.feature = 5, only.pos = FALSE, min.diff.pct = 0.10)
Batchedmarkers.mast_i_RNA_data_Naive$gene <- rownames(Batchedmarkers.mast_i_RNA_data_Naive)
##optional: add pct.fold = how large is the absolute difference in percentage?
Batchedmarkers.mast_i_RNA_data_Naive$pct.fold <- Batchedmarkers.mast_i_RNA_data_Naive$pct.1/Batchedmarkers.mast_i_RNA_data_Naive$pct.2

## Filter on adj.P-value
##change name according to test used (MAST, roc, negbinom, et.c)
Batchedmarkers.mast_i_RNA_data_Naive <-dplyr::filter(Batchedmarkers.mast_i_RNA_data_Naive, p_val_adj<0.05)
## Sort on logFC
Batchedmarkers.mast_i_RNA_data_Naive<-Batchedmarkers.mast_i_RNA_data_Naive[order(Batchedmarkers.mast_i_RNA_data_Naive$avg_log2FC, decreasing=T),]


##save as Robj
setwd("")
save(Batchedmarkers.mast_i_RNA_data_Naive, file='SupplData6_listDEG_fetal_CD4wTreg_scRNA_0.7_NaiveC1vsC0.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(Batchedmarkers.mast_i_RNA_data_Naive, file='SupplData6_listDEG_fetal_CD4wTreg_scRNA_0.7_NaiveC1vsC0.xlsx')
detach("package:openxlsx", unload=TRUE)

#### 7. TCR -------------------------------------------------------------------------------------
setwd("C:/Users/elise/Documents/PhD/Fetal 10X/TCR")
FSI_TCR <- read.csv('filtered_contig_annotations_FSIab.csv')  
FSP_TCR <- read.csv('filtered_contig_annotations_FSPab.csv') 
FLN_TCR <- read.csv('filtered_contig_annotations_FLNab.csv') 


FSI_TCR_combined_all <- combineTCR(FSI_TCR, samples='FSI', removeNA = F, removeMulti = T)
FSP_TCR_combined_all <- combineTCR(FSP_TCR, samples='FSP', removeNA = F, removeMulti = T)
FLN_TCR_combined_all <- combineTCR(FLN_TCR, samples='FLN', removeNA = F, removeMulti = T)

rownames(CD4_tissues_wTreg_2@meta.data)

FSP_TCR_combined_all$FSP$barcode <- paste0(str_remove(FSP_TCR_combined_all$FSP$barcode, 'FSP_'),'_1')
FLN_TCR_combined_all$FLN$barcode <- paste0(str_remove(FLN_TCR_combined_all$FLN$barcode, 'FLN_'),'_2')
FSI_TCR_combined_all$FSI$barcode <- paste0(str_remove(FSI_TCR_combined_all$FSI$barcode, 'FSI_'),'_3')

tissues_TCR_combined_all <- rbind(FSI_TCR_combined_all$FSI, FSP_TCR_combined_all$FSP, 
                                  FLN_TCR_combined_all$FLN)

tissues_TCR_combined_all <- cbind(tissues_TCR_combined_all, 
                                  CD4_tissues_wTreg_2@meta.data[match(tissues_TCR_combined_all$barcode,
                                                                      rownames(CD4_tissues_wTreg_2@meta.data)),
                                                                'Donor'])
colnames(tissues_TCR_combined_all)[13] <- 'Donor'

colnames(CD4_tissues_wTreg_2@meta.data)
##
CD4_tissues_wTreg_2_TCR <- CD4_tissues_wTreg_2
CD4_tissues_wTreg_2_TCR@meta.data <-  cbind(CD4_tissues_wTreg_2@meta.data, 
                                            tissues_TCR_combined_all[match(rownames(CD4_tissues_wTreg_2@meta.data),
                                                                           tissues_TCR_combined_all$barcode),1:12])

#### calculate CloneID for TCRab (clonal expansion only recognized if both are available, otherwise unique cloneID)
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  mutate(CloneID = ifelse(!(grepl('NA_', CTaa)|grepl('_NA', CTaa))&!is.na(Donor),as.integer(factor(CTaa)), NA)) 


# Compute the frequency of each CloneID within each Donor
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  group_by(Donor, CloneID) %>%
  mutate(Clone_Frequency = sum(!is.na(CloneID))) %>%
  ungroup()

# make sure that if there is no TCR info there is also no clone frequency
CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$CloneID), 'Clone_Frequency'] <- NA
CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped <- CD4_tissues_wTreg_2_TCR$Clone_Frequency
clonefreq_groups <- setNames(c("Single", "Rare (2)", rep("Small (3-5)",3), 
                               rep("Medium (6-20)",15), rep("Large (21-100)",20)),
                             as.character(1:40))  # Ensure correct names as character
CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped <- as.character(CD4_tissues_wTreg_2_TCR$Clone_Frequency)
CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped <- unname(clonefreq_groups[CD4_tissues_wTreg_2_TCR$Clone_Frequency])
CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped <- factor(CD4_tissues_wTreg_2_TCR$Clone_Frequency_grouped, 
                                                          levels = unique(clonefreq_groups))


#### same but only for TCR-beta
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  mutate(CloneID_TRB = ifelse(!is.na(Donor),as.integer(factor(cdr3_aa2)), NA))

CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  group_by(Donor, CloneID_TRB) %>%
  mutate(Clone_Frequency_TRB = sum(!is.na(CloneID_TRB))) %>%
  ungroup()

CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$CloneID_TRB), 'Clone_Frequency_TRB'] <- NA
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped <- CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB
clonefreq_groups <- setNames(c("Single", "Rare (2)", rep("Small (3-5)",3), 
                               rep("Medium (6-20)",15), rep("Large (21-100)",20)),
                             as.character(1:40))  # Ensure correct names as character
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped <- as.character(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB)
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped <- unname(clonefreq_groups[CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB])
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped <- factor(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRB_grouped, 
                                                              levels = unique(clonefreq_groups))

#### same but only for TCR-alpha
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  mutate(CloneID_TRA = ifelse(!is.na(Donor),as.integer(factor(cdr3_aa1)), NA))

CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  group_by(Donor, CloneID_TRA) %>%
  mutate(Clone_Frequency_TRA = sum(!is.na(CloneID_TRA))) %>%
  ungroup()

CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$CloneID_TRA), 'Clone_Frequency_TRA'] <- NA
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped <- CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA
clonefreq_groups <- setNames(c("Single", "Rare (2)", rep("Small (3-5)",3), 
                               rep("Medium (6-20)",15), rep("Large (21-100)",20)),
                             as.character(1:40))  # Ensure correct names as character
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped <- as.character(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA)
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped <- unname(clonefreq_groups[CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA])
CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped <- factor(CD4_tissues_wTreg_2_TCR$Clone_Frequency_TRA_grouped, 
                                                              levels = unique(clonefreq_groups))

#### Calculate how many unique donors each CloneID appears in
clone_publicity <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(Donor)) %>%
  group_by(CloneID) %>%
  summarise(Publicity = n_distinct(Donor), .groups = "drop")

# Merge back into meta.data
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  left_join(clone_publicity, by = "CloneID")

CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$CloneID), 'Publicity'] <- NA
CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$Donor), 'Publicity'] <- NA
#### same but for TCR-beta
clone_publicity_TRB <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(Donor)) %>%
  group_by(CloneID_TRB) %>%
  summarise(Publicity_TRB = n_distinct(Donor), .groups = "drop")

# Merge back into meta.data
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  left_join(clone_publicity_TRB, by = "CloneID_TRB")

CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$CloneID_TRB), 'Publicity_TRB'] <- NA
CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$Donor), 'Publicity_TRB'] <- NA

#### same but for TCR-alpha
clone_publicity_TRA <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(Donor)) %>%
  group_by(CloneID_TRA) %>%
  summarise(Publicity_TRA = n_distinct(Donor), .groups = "drop")

# Merge back into meta.data
CD4_tissues_wTreg_2_TCR@meta.data <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  left_join(clone_publicity_TRA, by = "CloneID_TRA")

CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$CloneID_TRA), 'Publicity_TRA'] <- NA
CD4_tissues_wTreg_2_TCR@meta.data[is.na(CD4_tissues_wTreg_2_TCR$Donor), 'Publicity_TRA'] <- NA


##### add TCR VDJ gene family information in separate columns
newcols_TRA_VDJ <- as.data.frame(do.call(rbind,str_split(CD4_tissues_wTreg_2_TCR$TCR1, '\\.')))
colnames(newcols_TRA_VDJ) <- c('V-gene TRA', 'J-gene TRA','C-gene TRA')
CD4_tissues_wTreg_2_TCR@meta.data <- cbind(CD4_tissues_wTreg_2_TCR@meta.data, newcols_TRA_VDJ)

newcols_TRB_VDJ <- as.data.frame(do.call(rbind,str_split(CD4_tissues_wTreg_2_TCR$TCR2, '\\.')))
colnames(newcols_TRB_VDJ) <- c('V-gene TRB','D-gene TRB', 'J-gene TRB','C-gene TRB')
CD4_tissues_wTreg_2_TCR@meta.data <- cbind(CD4_tissues_wTreg_2_TCR@meta.data, newcols_TRB_VDJ)

CD4_tissues_wTreg_2_TCR@meta.data[grep('*NA',CD4_tissues_wTreg_2_TCR$`D-gene TRB`),'D-gene TRB'] <- NA

##### add number of N-additions (used IMGT tool, allele information unknown, so they will use *01)
# IMGT_input_B <- as.data.frame(cbind(paste0(CD4_tissues_wTreg_2_TCR$barcode,', ',
#                            paste0(CD4_tissues_wTreg_2_TCR$`V-gene TRB`,'*?'),', ',
#                     paste0(CD4_tissues_wTreg_2_TCR$`J-gene TRB`,'*?')),
#                     CD4_tissues_wTreg_2_TCR$cdr3_nt2))
# colnames(IMGT_input_B)
# setwd("~/")
# 
# df_to_faa(IMGT_input_B[1:5000,],file="TRB_IMGTinput_part1.fasta")
# df_to_faa(IMGT_input_B[5001:10000,],file="TRB_IMGTinput_part2.fasta")
# df_to_faa(IMGT_input_B[10000:10470,],file="TRB_IMGTinput_part3.fasta")

setwd("~/")
IMGT_output <- read.csv('IMGT_output.csv', sep = ';')
IMGT_output$Ngc <- str_replace(IMGT_output$Ngc,'\\-','/')

CD4_tissues_wTreg_2_TCR$Nadditions_TRB <- IMGT_output[match(CD4_tissues_wTreg_2_TCR$barcode,
                                                            IMGT_output$Input),'Ngc']
CD4_tissues_wTreg_2_TCR$Nadditions_TRB <- as.numeric(str_remove(CD4_tissues_wTreg_2_TCR$Nadditions_TRB,'.*/'))

##### add number of N-additions (used IMGT tool, allele information unknown, so they will use *01)
# IMGT_input_A <- as.data.frame(cbind(paste0(CD4_tissues_wTreg_2_TCR$barcode,', ',
#                                          paste0(CD4_tissues_wTreg_2_TCR$`V-gene TRA`,'*?'),', ',
#                                          paste0(CD4_tissues_wTreg_2_TCR$`J-gene TRA`,'*?')),
#                                   CD4_tissues_wTreg_2_TCR$cdr3_nt1))
# colnames(IMGT_input_A)
# setwd("~/PhD/Fetal 10X/TCR")
# 
# df_to_faa(IMGT_input_A[1:5000,],file="TRA_IMGTinput_part1.fasta")
# df_to_faa(IMGT_input_A[5001:10000,],file="TRA_IMGTinput_part2.fasta")
# df_to_faa(IMGT_input_A[10000:13045,],file="TRA_IMGTinput_part3.fasta")

setwd("~/")
IMGT_output_TRA <- read.csv('IMGT_output_TRA.csv', sep = ';')
IMGT_output_TRA$Ngc <- str_replace(IMGT_output_TRA$Ngc,'\\-','/')

CD4_tissues_wTreg_2_TCR$Nadditions_TRA <- IMGT_output_TRA[match(CD4_tissues_wTreg_2_TCR$barcode,
                                                                IMGT_output_TRA$Input),'Ngc']
CD4_tissues_wTreg_2_TCR$Nadditions_TRA <- as.numeric(str_remove(CD4_tissues_wTreg_2_TCR$Nadditions_TRA,'.*/'))


mean(CD4_tissues_wTreg_2_TCR@meta.data[CD4_tissues_wTreg_2_TCR$tissues_res.0.7_manual=='4','Nadditions_TRA'], na.rm=T)

##### calculate CDR3 length
CD4_tissues_wTreg_2_TCR$TRA_length <- nchar(CD4_tissues_wTreg_2_TCR$cdr3_aa1)
CD4_tissues_wTreg_2_TCR$TRB_length <- nchar(CD4_tissues_wTreg_2_TCR$cdr3_aa2)

#### calculate clonal overlap between clusters (res.0.7) (per donor)
CD4_tissues_wTreg_2_TCR$tissues_res.0.7_manual <- as.character(CD4_tissues_wTreg_2_TCR$tissues_res.0.7_manual)

clone_per_group_res.0.7 <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID)&!is.na(Donor)) %>%  # Ignore NA clone IDs
  group_by(Donor, tissues_res.0.7_manual) %>%
  summarise(CloneIDs = list(unique(CloneID)),  # list of all the unique CloneIDs per group
            CloneIDs_total = length(unique(CloneID)), # number of unique clones per group
            .groups = "drop")

overlapping_clones_res.0.7 <- clone_per_group_res.0.7 %>%
  inner_join(clone_per_group_res.0.7, by = "Donor", suffix = c("_group1", "_group2")) %>% # create pairwise comparisons
  filter(tissues_res.0.7_manual_group1 < tissues_res.0.7_manual_group2) %>%  # Remove duplicate comparisons
  mutate(SharedClones = map2_int(CloneIDs_group1, CloneIDs_group2, ~ length(intersect(.x, .y)))) %>% # find overlapping cloneIDs between both lists of each comparison
  mutate(Union_Clones = CloneIDs_total_group1 + CloneIDs_total_group2 - SharedClones) %>% # number of non-shared unique clones per group
  mutate(JaccardIndex = SharedClones / Union_Clones)  # Compute Jaccard Index before summarizing 

overlapping_clones_summary_res.0.7 <- overlapping_clones_res.0.7 %>%
  mutate(comparison = paste0(tissues_res.0.7_manual_group1, '_', tissues_res.0.7_manual_group2)) %>% # to enable grouping by comparison
  group_by(comparison) %>%
  summarise(
    summedClones = sum(SharedClones),
    meanClones = mean(SharedClones),
    medianClones = median(SharedClones),
    meanJaccard = mean(JaccardIndex), 
    .groups = "drop"
  ) 

#### calculate clonal overlap between clusters (res.0.7) (per donor) - TRB
clone_per_group_res.0.7_TRB <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID_TRB)&!is.na(Donor)) %>%  # Ignore NA clone IDs
  group_by(Donor, tissues_res.0.7_manual) %>%
  summarise(CloneIDs = list(unique(CloneID_TRB)),  # list of all the unique CloneIDs per group
            CloneIDs_total = length(unique(CloneID_TRB)), # number of unique clones per group
            .groups = "drop")

overlapping_clones_res.0.7_TRB <- clone_per_group_res.0.7_TRB %>%
  inner_join(clone_per_group_res.0.7_TRB, by = "Donor", suffix = c("_group1", "_group2")) %>% # create pairwise comparisons
  filter(tissues_res.0.7_manual_group1 < tissues_res.0.7_manual_group2) %>%  # Remove duplicate comparisons
  mutate(SharedClones = map2_int(CloneIDs_group1, CloneIDs_group2, ~ length(intersect(.x, .y)))) %>% # find overlapping cloneIDs between both lists of each comparison
  mutate(Union_Clones = CloneIDs_total_group1 + CloneIDs_total_group2 - SharedClones) %>% # number of non-shared unique clones per group
  mutate(JaccardIndex = SharedClones / Union_Clones)  # Compute Jaccard Index before summarizing 

overlapping_clones_summary_res.0.7_TRB <- overlapping_clones_res.0.7_TRB %>%
  mutate(comparison = paste0(tissues_res.0.7_manual_group1, '_', tissues_res.0.7_manual_group2)) %>% # to enable grouping by comparison
  group_by(comparison) %>%
  summarise(
    summedClones = sum(SharedClones),
    meanClones = mean(SharedClones),
    medianClones = median(SharedClones),
    meanJaccard = mean(JaccardIndex), 
    .groups = "drop"
  ) 
#### calculate clonal overlap between clusters (res.0.7) (per donor) - TRA
clone_per_group_res.0.7_TRA <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID_TRA)&!is.na(Donor)) %>%  # Ignore NA clone IDs
  group_by(Donor, tissues_res.0.7_manual) %>%
  summarise(CloneIDs = list(unique(CloneID_TRA)),  # list of all the unique CloneIDs per group
            CloneIDs_total = length(unique(CloneID_TRA)), # number of unique clones per group
            .groups = "drop")

overlapping_clones_res.0.7_TRA <- clone_per_group_res.0.7_TRA %>%
  inner_join(clone_per_group_res.0.7_TRA, by = "Donor", suffix = c("_group1", "_group2")) %>% # create pairwise comparisons
  filter(tissues_res.0.7_manual_group1 < tissues_res.0.7_manual_group2) %>%  # Remove duplicate comparisons
  mutate(SharedClones = map2_int(CloneIDs_group1, CloneIDs_group2, ~ length(intersect(.x, .y)))) %>% # find overlapping cloneIDs between both lists of each comparison
  mutate(Union_Clones = CloneIDs_total_group1 + CloneIDs_total_group2 - SharedClones) %>% # number of non-shared unique clones per group
  mutate(JaccardIndex = SharedClones / Union_Clones)  # Compute Jaccard Index before summarizing 

overlapping_clones_summary_res.0.7_TRA <- overlapping_clones_res.0.7_TRA %>%
  mutate(comparison = paste0(tissues_res.0.7_manual_group1, '_', tissues_res.0.7_manual_group2)) %>% # to enable grouping by comparison
  group_by(comparison) %>%
  summarise(
    summedClones = sum(SharedClones),
    meanClones = mean(SharedClones),
    medianClones = median(SharedClones),
    meanJaccard = mean(JaccardIndex), 
    .groups = "drop"
  ) 

#### calculate clonal overlap between clusters and tissues (per donor) - full
clone_per_group_res.0.7_tissue <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID)&!is.na(Donor)) %>%  # Ignore NA clone IDs
  group_by(Donor, tissues_res.0.7_manual, orig.ident) %>%
  summarise(CloneIDs = list(unique(CloneID)),  # list of all the unique CloneIDs per group
            CloneIDs_total = length(unique(CloneID)), # number of unique clones per group
            .groups = "drop")%>%
  mutate(tissues_res.0.7_manual = as.character(tissues_res.0.7_manual))

overlapping_clones_res.0.7_tissue <- clone_per_group_res.0.7_tissue %>%
  inner_join(clone_per_group_res.0.7_tissue, by = "Donor", suffix = c("_group1", "_group2")) %>% # create pairwise comparisons
  filter(tissues_res.0.7_manual_group1 <= tissues_res.0.7_manual_group2) %>%  # Remove duplicate comparisons
  mutate(SharedClones = map2_int(CloneIDs_group1, CloneIDs_group2, ~ length(intersect(.x, .y)))) %>% # find overlapping cloneIDs between both lists of each comparison
  mutate(Union_Clones = CloneIDs_total_group1 + CloneIDs_total_group2 - SharedClones) %>% # number of non-shared unique clones per group
  mutate(JaccardIndex = SharedClones / Union_Clones)  # Compute Jaccard Index before summarizing 

overlapping_clones_summary_res.0.7_tissue <- overlapping_clones_res.0.7_tissue %>%
  mutate(comparison = paste0(tissues_res.0.7_manual_group1, '_',orig.ident_group1,'-', 
                             tissues_res.0.7_manual_group2,'_',orig.ident_group2)) %>% # to enable grouping by comparison
  group_by(comparison) %>%
  summarise(
    summedClones = sum(SharedClones),
    meanClones = mean(SharedClones),
    medianClones = median(SharedClones),
    meanJaccard = mean(JaccardIndex), 
    .groups = "drop"
  ) 


#### calculate clonal overlap between clusters and tissues (per donor) - TRB
clone_per_group_res.0.7_tissue_TRB <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID_TRB)&!is.na(Donor)) %>%  # Ignore NA clone IDs
  group_by(Donor, tissues_res.0.7_manual, orig.ident) %>%
  summarise(CloneIDs = list(unique(CloneID_TRB)),  # list of all the unique CloneIDs per group
            CloneIDs_total = length(unique(CloneID_TRB)), # number of unique clones per group
            .groups = "drop") %>%
  mutate(tissues_res.0.7_manual = as.character(tissues_res.0.7_manual))

overlapping_clones_res.0.7_tissue_TRB <- clone_per_group_res.0.7_tissue_TRB %>%
  inner_join(clone_per_group_res.0.7_tissue_TRB, by = "Donor", suffix = c("_group1", "_group2")) %>% # create pairwise comparisons
  filter(tissues_res.0.7_manual_group1 <= tissues_res.0.7_manual_group2) %>%  # Remove duplicate comparisons
  mutate(SharedClones = map2_int(CloneIDs_group1, CloneIDs_group2, ~ length(intersect(.x, .y)))) %>% # find overlapping cloneIDs between both lists of each comparison
  mutate(Union_Clones = CloneIDs_total_group1 + CloneIDs_total_group2 - SharedClones) %>% # number of non-shared unique clones per group
  mutate(JaccardIndex = SharedClones / Union_Clones)  # Compute Jaccard Index before summarizing 

overlapping_clones_summary_res.0.7_tissue_TRB <- overlapping_clones_res.0.7_tissue_TRB %>%
  mutate(comparison = paste0(tissues_res.0.7_manual_group1, '_',orig.ident_group1,'-', 
                             tissues_res.0.7_manual_group2,'_',orig.ident_group2)) %>% # to enable grouping by comparison
  group_by(comparison) %>%
  summarise(
    summedClones = sum(SharedClones),
    meanClones = mean(SharedClones),
    medianClones = median(SharedClones),
    meanJaccard = mean(JaccardIndex), 
    .groups = "drop"
  ) 

#### calculate clonal overlap between clusters & tissues (per donor) - TRA
clone_per_group_res.0.7_tissue_TRA <- CD4_tissues_wTreg_2_TCR@meta.data %>%
  filter(!is.na(CloneID_TRA)&!is.na(Donor)) %>%  # Ignore NA clone IDs
  group_by(Donor, tissues_res.0.7_manual, orig.ident) %>%
  summarise(CloneIDs = list(unique(CloneID_TRA)),  # list of all the unique CloneIDs per group
            CloneIDs_total = length(unique(CloneID_TRA)), # number of unique clones per group
            .groups = "drop") %>%
  mutate(tissues_res.0.7_manual = as.character(tissues_res.0.7_manual))

overlapping_clones_res.0.7_tissue_TRA <- clone_per_group_res.0.7_tissue_TRA %>%
  inner_join(clone_per_group_res.0.7_tissue_TRA, by = "Donor", suffix = c("_group1", "_group2")) %>% # create pairwise comparisons
  filter(tissues_res.0.7_manual_group1 <= tissues_res.0.7_manual_group2) %>%  # Remove duplicate comparisons
  mutate(SharedClones = map2_int(CloneIDs_group1, CloneIDs_group2, ~ length(intersect(.x, .y)))) %>% # find overlapping cloneIDs between both lists of each comparison
  mutate(Union_Clones = CloneIDs_total_group1 + CloneIDs_total_group2 - SharedClones) %>% # number of non-shared unique clones per group
  mutate(JaccardIndex = SharedClones / Union_Clones)  # Compute Jaccard Index before summarizing 

overlapping_clones_summary_res.0.7_tissue_TRA <- overlapping_clones_res.0.7_tissue_TRA %>%
  mutate(comparison = paste0(tissues_res.0.7_manual_group1, '_',orig.ident_group1,'-', 
                             tissues_res.0.7_manual_group2,'_',orig.ident_group2)) %>% # to enable grouping by comparison
  group_by(comparison) %>%
  summarise(
    summedClones = sum(SharedClones),
    meanClones = mean(SharedClones),
    medianClones = median(SharedClones),
    meanJaccard = mean(JaccardIndex), 
    .groups = "drop"
  ) 

