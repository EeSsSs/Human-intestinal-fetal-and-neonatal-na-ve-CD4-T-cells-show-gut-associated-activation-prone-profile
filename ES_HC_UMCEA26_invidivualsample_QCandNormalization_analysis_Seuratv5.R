#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
  library(Seurat)
  
library(dplyr)
  
  library(sctransform)
  library(ggplot2)
  library(flowCore)
  library(Matrix)
library(glmGamPoi)

### 1. Setting working directory and loading data----------------------------------------------------------------
 
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/2025-2412-UMC-EA-s023-s030/raw_count_tables/non_poisson_corrected/UMC-EA-s026-raw")
  
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
  sum(rawData==0) #13563656
  
  # Calculate gene counts -> so in how many cells a gene is counted
  geneCounts <- apply(rawData,1,function (x){sum(x>0)}) # 1 indicates rows
  min(geneCounts) # 0
  max(geneCounts) # 376
  head(geneCounts)
  
  mean(geneCounts) #5.418431
  
  length(geneCounts[geneCounts<3]) #27484
  length(geneCounts[geneCounts==0]) #23114
  
  # Calculate cell counts -> so how many genes are expressed per cell
  cellCounts<-apply(rawData,2,function (x){sum(x>0)}) #2 indicates columns
  
  min(cellCounts) # 31
  max(cellCounts) # 1270
  
  mean(cellCounts) #527.4468
  
  length(cellCounts[cellCounts<300]) #76
  length(cellCounts[cellCounts == 0]) # 0
  
  # Remove the values created to get an overview
  rm(geneCounts)
  rm(cellCounts)
  
### 4. Read Flow data -------------------------------------------------------------------------------------------------------------
  setwd("T:/cff-data/WKZ/Group-vanWijk/Elise/Sort experiments/HINT/240926_HINT140_blood_scRNAseq_sort/")
  ## 
  EA026_flow_x <-read.FCS('HINT140_blood_INX_Tube_001_013_compensated.fcs', truncate_max_range = FALSE)
  
  View(EA026_flow_x@exprs)
  
  EA026_flow <- getIndexSort(EA026_flow_x)
  ## see barcode file: scRNAseq data are ordered from A1->A24 and then B1->B24
  ## but cells were sorted (+index sort data) from A1->A24 and then B24->B1
  ## xloc = rows yloc = columns
  EA026_flow <- EA026_flow[order(EA026_flow$XLoc,EA026_flow$YLoc),]
  hist(EA026_flow$YLoc)
 ## many rowss missing
  dummy_df <- EA026_flow
  dummy_df[1:255,1:26] <- NA
  EA026_flow <- rbind(dummy_df[1:7,],EA026_flow[1,],dummy_df[1:18,],EA026_flow[2,],dummy_df[1:3,],
                      EA026_flow[3,],dummy_df[1:10,],EA026_flow[4:5,],dummy_df[1:13,],
                      EA026_flow[6,], dummy_df[1,],EA026_flow[7,],dummy_df[1,],
                      EA026_flow[8:9,],dummy_df[1:60,],EA026_flow[10,],dummy_df[1:10,],
                      EA026_flow[11,],dummy_df[1:4,],EA026_flow[12,],dummy_df[1,],
                      EA026_flow[13:15,],dummy_df[1:6,],EA026_flow[16,],dummy_df[1:16,],
                      EA026_flow[17:19,],dummy_df[1:4,],EA026_flow[20,],dummy_df[1:2,],
                      EA026_flow[21,],dummy_df[1:5,],EA026_flow[22,],dummy_df[1,],
                      EA026_flow[23,],dummy_df[1:4,],EA026_flow[24,],dummy_df[1:7,],
                      EA026_flow[25,],dummy_df[1:2,],EA026_flow[26,],dummy_df[1:7,],
                      EA026_flow[27,],dummy_df[1:12,],EA026_flow[28,],dummy_df[1:8,],
                      EA026_flow[29:32,],dummy_df[1:7,],EA026_flow[33,],dummy_df[1:8,],
                      EA026_flow[34,],dummy_df[1:5,],EA026_flow[35,],dummy_df[1:16,],
                      EA026_flow[36,],dummy_df[1:6,],EA026_flow[37:39,],dummy_df[1:4,],
                      EA026_flow[40,],dummy_df[1:3,],EA026_flow[41:42,],dummy_df[1:2,],
                      EA026_flow[43,],dummy_df[1:3,],EA026_flow[44,],dummy_df[1:12,],
                      EA026_flow[45,],dummy_df[1:3,],EA026_flow[46,],dummy_df[1:13,],
                      EA026_flow[47,],dummy_df[1:2,],EA026_flow[48,],dummy_df[1:10,],
                      EA026_flow[49,],dummy_df[1:2,],EA026_flow[50,],dummy_df[1:20,],
                      EA026_flow[51:52,],dummy_df[1:6,])
  
  rownames(EA026_flow)  <- colnames(count_matrix)
  colnames(EA026_flow)[8:23] <- EA026_flow_x@parameters@data$desc[7:22]
  colnames(EA026_flow)[is.na(colnames(EA026_flow))] <- 'empty'
  ##check compensation
  plot(EA026_flow$CD4~EA026_flow$CD31)
  plot(EA026_flow$CD95~EA026_flow$CD8)
  plot(EA026_flow$CCR7~EA026_flow$CD127)
  
  View(EA026_flow)
  dim(EA026_flow)
  
  
### 5. Seurat object with flow data -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data. remove genes that are expressed in less than 3 cells
  EA026 <- CreateSeuratObject(counts = rawData, project = "EA026", min.cells=3)
  EA026
  #add flowdata as metadata
  EA026@meta.data$flowCD31 <- EA026_flow$CD31
  EA026@meta.data$flowCD25 <- EA026_flow$CD25
  EA026@meta.data$flowCD45RA <- EA026_flow$CD45RA
  EA026@meta.data$flowCD27 <- EA026_flow$CD27
  EA026@meta.data$flowCCR7 <- EA026_flow$CCR7
  EA026@meta.data$flowCD127 <- EA026_flow$CD127
  EA026@meta.data$flowCD69 <- EA026_flow$CD69
  EA026 

### 6. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/QC")

  ## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  EA026 <- PercentageFeatureSet(EA026, pattern = "^MT\\.", col.name = 'percent.mt')
  EA026 <- PercentageFeatureSet(EA026, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(count_matrix)[grep('MT', rownames(count_matrix))]) #somewhere in name
  length(rownames(count_matrix)[grep('^MT\\.', rownames(count_matrix))]) #start of name
  
  # Show QC metrics for the first 5 cells
  head(EA026@meta.data, 5)
  
  #save(EA026, file = "EA026.Robj")
  
  ## Plot histograms of QC data pre-filtering
  toPlot <- EA026@meta.data
  
  # percent.mito
  png(file="EA026_pre_percMito.png", width=850)
  par(mfrow=c(1,2))
  tmp <- toPlot[order(toPlot$percent.mt),]
  hist(tmp$percent.mt, breaks=30)
  barplot(tmp$percent.mt)
  dev.off()
  
  # nFeature_RNA
  png(file="EA026_pre_nFeature_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nFeature_RNA),]
  hist(tmp$nFeature_RNA, breaks=30)
  barplot(tmp$nFeature_RNA)
  dev.off()
  
  # nCount_RNA
  png(file="EA026_pre_nCount_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nCount_RNA),]
  hist(tmp$nCount_RNA, breaks=30)
  barplot(tmp$nCount_RNA)
  dev.off()
  
  # Visualize QC metrics as a violin plot
  pdf(file="EA026_beforeQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA026, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="EA026_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(EA026, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(EA026, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(EA026, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  plot1 + plot2 + plot3
  dev.off()
  
  #With filtering cut-offs 
  pdf(file="EA026_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
  plot1 <- ggplot(EA026@meta.data, aes(EA026$nCount_RNA, EA026$nFeature_RNA))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 2000)+ annotate("text",x=2250,y=1,label=c("2000"),hjust=0, size=2.8)
    plot2 <- ggplot(EA026@meta.data, aes(EA026$nCount_RNA, EA026$percent.mt))+geom_point(size=.9)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 300)+ annotate("text",x=380,y=1,label=c("300"),hjust=0, size=2.8)+
    geom_hline(yintercept = 16)+ annotate("text",x=700,y=16.5,label=c("16"),vjust=0, size=2.8)
  plot3 <- ggplot(EA026@meta.data, aes(EA026$nFeature_RNA, EA026$percent.mt))+geom_point(size=.9)+theme_bw()+ 
    geom_vline(xintercept = 250)+annotate("text",x=320,y=1,label=c("250"),hjust=0, size=2.8)+
    geom_hline(yintercept = 16)+ annotate("text",x=500,y=16.5,label=c("16"),vjust=0, size=2.8)
  plot4 <- ggplot(EA026@meta.data, aes(EA026$nCount_RNA, EA026$percent.ribo))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 300)+ annotate("text",x=380,y=1,label=c("300"),hjust=0, size=2.8)+
    geom_hline(yintercept = 10)+ annotate("text",x=4000,y=11,label=c("10"),vjust=0, size=2.8)
  plot5 <- ggplot(EA026@meta.data, aes(EA026$percent.mt, EA026$percent.ribo))+geom_point(size=.9)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(EA026, expression = nFeature_RNA > 250 & nFeature_RNA < 2000 &nCount_RNA > 300 & percent.mt < 16 & percent.ribo > 10)
  length(selected) #How many cells are left
  
  # Filter cells
  EA026 <- subset(EA026, subset = nFeature_RNA > 250 & nFeature_RNA < 2000 &nCount_RNA > 300 & percent.mt < 16 & percent.ribo > 10)
  
  # Visualize QC metrics post-filtering
  pdf(file="EA026_afterQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA026, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
### 7. SCTransform and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  EA026 <- NormalizeData(EA026, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  EA026 <- CellCycleScoring(object = EA026, s.features = s.genes, g2m.features = g2m.genes, 
                            set.ident = FALSE)
  head(x = EA026@meta.data)
  
  
  ## SCTransform is a code for UMI corrected matrices
  #optional (new version, don't know if it's default already, haven't tried yet): vst.flavor='v2'
  ## I think I don't want to regress out cell-cycle, as you don't except many cells to be dividing, and if so, that's actually interesting info
  EA026_mt <- SCTransform(object = EA026, verbose = TRUE, 
                       vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
  
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq")
  save(EA026_mt, file = 'EA026SCT_pctmt_nocellcylce.Robj')
  
  #with percent.ribo
  EA026_mtrb <- SCTransform(object = EA026, verbose = TRUE, 
                       vars.to.regress = c("percent.mt", "percent.ribo"), return.only.var.genes = FALSE)
  save(EA026_mtrb, file = 'EA026SCT_pctmtribo_nocellcylce.Robj')
  
  #with cell cycle
  EA026_mtcc <- SCTransform(object = EA026, verbose = TRUE, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), return.only.var.genes = FALSE)
  save(EA026_mtcc, file = 'EA026SCT_pctmtcc_cellcylce.Robj')
  
### 8. Optional: First check of clustering --------------------------------------------------------------------------------------------
 ## weird that the KRT genes are not there now (specific to that sequencing run?)
  EA026_mt <- RunPCA(EA026_mt, verbose = TRUE)
  EA026_mt <- RunUMAP(EA026_mt, dims = 1:30)
  
  
  ## 
  EA026_mt <- FindNeighbors(object = EA026_mt, reduction = "pca", dims = 1:30)
  ## not more clusters than with naive plate
  EA026_mt <- FindClusters(object = EA026_mt, resolution = .9)
  
  DimPlot(object = EA026_mt, reduction = "umap", pt.size=1.2)

  FeaturePlot(EA026_mt, 'flowCD31', pt.size=2)
  FeaturePlot(EA026_mt, 'flowCD127', pt.size=2)
  VlnPlot(EA026_mt, 'flowCD31')
  
  ## better: less MT and RP genes
  EA026_mtrb <- RunPCA(EA026_mtrb, verbose = TRUE)
  EA026_mtrb <- RunUMAP(EA026_mtrb, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA026_mtrb <- FindNeighbors(object = EA026_mtrb, reduction = "pca", dims = 1:30)
  ##
  EA026_mtrb <- FindClusters(object = EA026_mtrb, resolution = .9)
  
  DimPlot(object = EA026_mtrb, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA026_mtrb, 'flowCD31', pt.size=2)
  FeaturePlot(EA026_mtrb, 'flowCD127', pt.size=2)
  VlnPlot(EA026_mtrb, 'flowCD31')
  
  ## similar to mt
  EA026_mtcc <- RunPCA(EA026_mtcc, veccose = TRUE)
  EA026_mtcc <- RunUMAP(EA026_mtcc, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA026_mtcc <- FindNeighbors(object = EA026_mtcc, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA026_mtcc <- FindClusters(object = EA026_mtcc, resolution = .95)
  
  DimPlot(object = EA026_mtcc, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA026_mtcc, 'flowCD31', pt.size=2)
  FeaturePlot(EA026_mtcc, 'flowCD127', pt.size=2)
  VlnPlot(EA026_mtcc, 'flowCD31')