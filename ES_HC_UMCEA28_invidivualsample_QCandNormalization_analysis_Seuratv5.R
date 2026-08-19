#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
  library(Seurat)
  
library(dplyr)
  
  library(sctransform)
  library(ggplot2)
  library(flowCore)
  library(Matrix)
library(glmGamPoi)

### 1. Setting working directory and loading data----------------------------------------------------------------
 
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/2025-2412-UMC-EA-s023-s030/raw_count_tables/non_poisson_corrected/UMC-EA-s028-raw")
  
  # Read in matrix.mtx
  counts <- readMM("matrix.mtx")
  
  # Read in genes.tsv
  genes <- read.table("features.tsv")
  gene_ids <- genes$V2
  
  # Read in barcodes.tsv
  cells <- read.table("barcodes.tsv")
  cell_ids <- cells$V1
  
  # Make the column names as the cell IDs and the row names as the gene IDs
  rownames(counts) <- gene_ids
  colnames(counts) <- cell_ids
  count_matrix <- as.data.frame(counts)
  
  #couple the rigth cell ID to the barcodes
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq")
  
  barcodes <- read.table("celseq2_barcodes_position_wells.tsv", skip=1, col.names = c('Well position', 'Name', 'Barcode'))
  
  colnames(count_matrix) <- match(colnames(count_matrix),barcodes$Barcode)
  count_matrix <- count_matrix[,order(as.numeric(colnames(count_matrix)))]
  
  #Check count matrix
  dim(count_matrix) # 36693  384
  

### 2. Cleaning dataframe-----------------------------------------------------------------------------------
  
   # Remove columns of spike ins using ERCCs (=External RNA Control Consortium) controls
   # O21-O24 and P21-P24 -> see barcodes file for the colnumbers
  ##from E5 till L17 nonTreg
  count_matrix <- count_matrix[,-c(357:360, 381:384)]
    dim(count_matrix) 
  
  # Check the names and number of spike ins in the gene list
    ##Should be 92 or lower (because not all spikeins have to be transcribed)
  rownames(count_matrix)[grep('ERCC\\.', rownames(count_matrix))] #gives the list of all ERCCs
  length(rownames(count_matrix)[grep('ERCC\\.', rownames(count_matrix))]) # gives nr of rows with ERCC somewhere in the name
  length(rownames(count_matrix)[grep('^ERCC\\.', rownames(count_matrix))]) # gives nr of rows with ERRC at the start of the name
  
   # Remove all the rows of spike ins from the gene list
  rawData <- count_matrix[-grep('^ERCC\\.', rownames(count_matrix)),]
  
   # Check if and how many rows are removed
  length(grep('ERCC\\.', rownames(rawData))) #Should be 0, because all rownames with ERCC. should be removed by now
  length(rownames(rawData))

  
### 3. Broad overview of the data-------------------------------------------------------------------------------------

  # Look at first five rows and first five columns
  rawData[1:5,1:5]
  
  # Check dimensions of your data
  dim(rawData) #36601  376
  
  # Number (sum) of zero's in the dataset
  sum(rawData==0) #13627887
  
  # Calculate gene counts -> so in how many cells a gene is counted
  geneCounts <- apply(rawData,1,function (x){sum(x>0)}) # 1 indicates rows
  min(geneCounts) # 0
  max(geneCounts) # 365
  head(geneCounts)
  
  mean(geneCounts) #3.663534
  
  length(geneCounts[geneCounts<3]) #28245
  length(geneCounts[geneCounts==0]) #23212
  
  # Calculate cell counts -> so how many genes are expressed per cell
  cellCounts<-apply(rawData,2,function (x){sum(x>0)}) #2 indicates columns
  
  min(cellCounts) # 26
  max(cellCounts) # 2657
  
  mean(cellCounts) #356.6197
  
  length(cellCounts[cellCounts<300]) #211
  length(cellCounts[cellCounts == 0]) # 0
  
  # Remove the values created to get an overview
  rm(geneCounts)
  rm(cellCounts)
  
### 4. Read Flow data -------------------------------------------------------------------------------------------------------------
  setwd("T:/cff-data/WKZ/Group-vanWijk/Elise/Sort experiments/HINT/241129_HINT148_gut_scRNAseq_sort/")
  ## 
  EA028_flow_x_1 <-read.FCS('HINT148_gut_INX_Tube_001_013_compensated.fcs', truncate_max_range = FALSE)
  EA028_flow_x_2 <-read.FCS('HINT148_gut_INX_Tube_002_014_compensated.fcs', truncate_max_range = FALSE)
  View(EA028_flow_x@exprs)
  
  EA028_flow_1 <- getIndexSort(EA028_flow_x_1)
  EA028_flow_2 <- getIndexSort(EA028_flow_x_2)
  EA028_flow <- rbind(EA028_flow_1, EA028_flow_2)
  ## see barcode file: scRNAseq data are ordered from A1->A24 and then B1->B24
  ## but cells were sorted (+index sort data) from A1->A24 and then B24->B1
  ## xloc = rows yloc = columns
  EA028_flow <- EA028_flow[order(EA028_flow$XLoc,EA028_flow$YLoc),]
  hist(EA028_flow$YLoc, breaks=seq(0,23,1))
  hist(EA028_flow$XLoc, breaks=seq(0,15,1))
 ## many rowss missing
  dummy_df <- EA028_flow
  dummy_df[1:255,1:26] <- NA
  EA028_flow <- rbind(EA028_flow[1:168,],dummy_df[1:4,],EA028_flow[169:188,])
  
  rownames(EA028_flow)  <- colnames(count_matrix)[1:192]
  colnames(EA028_flow)[8:23] <- EA028_flow_x@parameters@data$desc[7:22]
  colnames(EA028_flow)[is.na(colnames(EA028_flow))] <- 'empty'
  ##check compensation
  plot(EA028_flow$CD4~EA028_flow$CD31)
  plot(EA028_flow$CD95~EA028_flow$CD8)
  plot(EA028_flow$CCR7~EA028_flow$CD127)
  
  View(EA028_flow)
  dim(EA028_flow)
  
  
### 5. Seurat object with flow data -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data. remove genes that are expressed in less than 3 cells
  EA028 <- CreateSeuratObject(counts = rawData, project = "EA028", min.cells=3)
  EA028
  #add flowdata as metadata
  EA028@meta.data$flowCD31 <- c(EA028_flow$CD31,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CD31)))
  EA028@meta.data$flowCD25 <- c(EA028_flow$CD25,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CD25)))
  EA028@meta.data$flowCD45RA <- c(EA028_flow$CD45RA,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CD45RA)))
  EA028@meta.data$flowCD27 <- c(EA028_flow$CD27,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CD27)))
  EA028@meta.data$flowCCR7 <- c(EA028_flow$CCR7,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CCR7)))
  EA028@meta.data$flowCD127 <- c(EA028_flow$CD127,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CD127)))
  EA028@meta.data$flowCD69 <- c(EA028_flow$CD69,rep(NA, length(EA028@meta.data$orig.ident)-length(EA028_flow$CD69)))
  EA028 

### 6. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/QC")

  ## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  EA028 <- PercentageFeatureSet(EA028, pattern = "^MT\\.", col.name = 'percent.mt')
  EA028 <- PercentageFeatureSet(EA028, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(count_matrix)[grep('MT', rownames(count_matrix))]) #somewhere in name
  length(rownames(count_matrix)[grep('^MT\\.', rownames(count_matrix))]) #start of name
  
  # Show QC metrics for the first 5 cells
  head(EA028@meta.data, 5)
  
  #save(EA028, file = "EA028.Robj")
  
  ## Plot histograms of QC data pre-filtering
  toPlot <- EA028@meta.data
  
  # percent.mito
  png(file="EA028_pre_percMito.png", width=850)
  par(mfrow=c(1,2))
  tmp <- toPlot[order(toPlot$percent.mt),]
  hist(tmp$percent.mt, breaks=30)
  barplot(tmp$percent.mt)
  dev.off()
  
  # nFeature_RNA
  png(file="EA028_pre_nFeature_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nFeature_RNA),]
  hist(tmp$nFeature_RNA, breaks=30)
  barplot(tmp$nFeature_RNA)
  dev.off()
  
  # nCount_RNA
  png(file="EA028_pre_nCount_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nCount_RNA),]
  hist(tmp$nCount_RNA, breaks=30)
  barplot(tmp$nCount_RNA)
  dev.off()
  
  # Visualize QC metrics as a violin plot
  pdf(file="EA028_beforeQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA028, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="EA028_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(EA028, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(EA028, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(EA028, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  plot1 + plot2 + plot3
  dev.off()
  
  #With filtering cut-offs 
  pdf(file="EA028_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
  plot1 <- ggplot(EA028@meta.data, aes(EA028$nCount_RNA, EA028$nFeature_RNA))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 2000)+ annotate("text",x=2300,y=1,label=c("2000"),hjust=0, size=2.8)
    plot2 <- ggplot(EA028@meta.data, aes(EA028$nCount_RNA, EA028$percent.mt))+geom_point(size=.9)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 350)+ annotate("text",x=380,y=1,label=c("350"),hjust=0, size=2.8)+
    geom_hline(yintercept = 13)+ annotate("text",x=700,y=13.5,label=c("13"),vjust=0, size=2.8)
  plot3 <- ggplot(EA028@meta.data, aes(EA028$nFeature_RNA, EA028$percent.mt))+geom_point(size=.9)+theme_bw()+ 
    geom_vline(xintercept = 300)+annotate("text",x=320,y=1,label=c("300"),hjust=0, size=2.8)+
    geom_hline(yintercept = 13)+ annotate("text",x=500,y=13.5,label=c("13"),vjust=0, size=2.8)
  plot4 <- ggplot(EA028@meta.data, aes(EA028$nCount_RNA, EA028$percent.ribo))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 350)+ annotate("text",x=380,y=1,label=c("350"),hjust=0, size=2.8)+
    geom_hline(yintercept = 10)+ annotate("text",x=4000,y=11,label=c("10"),vjust=0, size=2.8)
  plot5 <- ggplot(EA028@meta.data, aes(EA028$percent.mt, EA028$percent.ribo))+geom_point(size=.9)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(EA028, expression = nFeature_RNA > 300 & nFeature_RNA < 2000 &nCount_RNA > 350 & percent.mt < 13 & percent.ribo>10)
  length(selected) #How many cells are left
  
  # Filter cells
  EA028 <- subset(EA028, subset = nFeature_RNA > 300 & nFeature_RNA < 2000 &nCount_RNA > 350 & percent.mt < 13 & percent.ribo>1)
  
  # Visualize QC metrics post-filtering
  pdf(file="EA028_afterQC_Vlns.pdf", width = 13, height = 8)
  VlnPlot(EA028, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
### 7. SCTransform and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  EA028 <- NormalizeData(EA028, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  EA028 <- CellCycleScoring(object = EA028, s.features = s.genes, g2m.features = g2m.genes, 
                            set.ident = FALSE)
  head(x = EA028@meta.data)
  
  
  ## SCTransform is a code for UMI corrected matrices
  #optional (new version, don't know if it's default already, haven't tried yet): vst.flavor='v2'
  ## I think I don't want to regress out cell-cycle, as you don't except many cells to be dividing, and if so, that's actually interesting info
  EA028_mt <- SCTransform(object = EA028, verbose = TRUE, 
                       vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
  
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
  save(EA028_mt, file = 'EA028SCT_pctmt_nocellcylce.Robj')
  
  #with percent.ribo
  EA028_mtrb <- SCTransform(object = EA028, verbose = TRUE, 
                       vars.to.regress = c("percent.mt", "percent.ribo"), return.only.var.genes = FALSE)
  save(EA028_mtrb, file = 'EA028SCT_pctmtribo_nocellcylce.Robj')
  
  #with cell cycle
  EA028_mtcc <- SCTransform(object = EA028, verbose = TRUE, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), return.only.var.genes = FALSE)
  save(EA028_mtcc, file = 'EA028SCT_pctmtcc_cellcylce.Robj')
  
### 8. Optional: First check of clustering --------------------------------------------------------------------------------------------
 ## weird that the KRT genes are not there now (specific to that sequencing run?)
  EA028_mt <- RunPCA(EA028_mt, verbose = TRUE)
  EA028_mt <- RunUMAP(EA028_mt, dims = 1:30)
  
  
  ## 
  EA028_mt <- FindNeighbors(object = EA028_mt, reduction = "pca", dims = 1:30)
  ## not more clusters than with naive plate
  EA028_mt <- FindClusters(object = EA028_mt, resolution = .9)
  
  DimPlot(object = EA028_mt, reduction = "umap", pt.size=1.2)

  FeaturePlot(EA028_mt, 'flowCD31', pt.size=2)
  FeaturePlot(EA028_mt, 'flowCD127', pt.size=2)
  VlnPlot(EA028_mt, 'flowCD31')
  
  ## better: less MT and RP genes
  EA028_mtrb <- RunPCA(EA028_mtrb, verbose = TRUE)
  EA028_mtrb <- RunUMAP(EA028_mtrb, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA028_mtrb <- FindNeighbors(object = EA028_mtrb, reduction = "pca", dims = 1:30)
  ##
  EA028_mtrb <- FindClusters(object = EA028_mtrb, resolution = .9)
  
  DimPlot(object = EA028_mtrb, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA028_mtrb, 'flowCD31', pt.size=2)
  FeaturePlot(EA028_mtrb, 'flowCD127', pt.size=2)
  VlnPlot(EA028_mtrb, 'flowCD31')
  
  ## similar to mt
  EA028_mtcc <- RunPCA(EA028_mtcc, veccose = TRUE)
  EA028_mtcc <- RunUMAP(EA028_mtcc, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA028_mtcc <- FindNeighbors(object = EA028_mtcc, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA028_mtcc <- FindClusters(object = EA028_mtcc, resolution = .95)
  
  DimPlot(object = EA028_mtcc, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA028_mtcc, 'flowCD31', pt.size=2)
  FeaturePlot(EA028_mtcc, 'flowCD127', pt.size=2)
  VlnPlot(EA028_mtcc, 'flowCD31')