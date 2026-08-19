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

load('EA012SCT_pctmtribo_nocellcylce.Robj')
load('EA014SCT_pctmtribo_nocellcylce.Robj')
load('EA016SCT_pctmtribo_nocellcylce.Robj')

load('EA021SCT_pctmtribo_nocellcylce.Robj')

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

#### Paper Figures ####
setwd("")

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

setwd("")
C2genes <- read_excel("SupplData7_listDEG_childrenadult_CD4Naive_scRNA_0.6_C2vsC10.xlsx")
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
setwd("")
save(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata, file='SupplData7_listDEG_childrenadult_CD4Naive_scRNA_0.6_C2vsC10.xlsx.Robj')

## Write to Excel
library('openxlsx')
write.xlsx(Naive.gutblood.merged.all_res0.6_manual_markers_MASTRNAdata, file='SupplData7_listDEG_childrenadult_CD4Naive_scRNA_0.6_C2vsC10.xlsx.xlsx')

