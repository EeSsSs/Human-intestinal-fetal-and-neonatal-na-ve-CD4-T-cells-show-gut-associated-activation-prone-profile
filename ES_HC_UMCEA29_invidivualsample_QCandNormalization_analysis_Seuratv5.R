#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
  library(Seurat)
  
library(dplyr)
  
  library(sctransform)
  library(ggplot2)
  library(flowCore)
  library(Matrix)
library(glmGamPoi)

### 1. Setting working directory and loading data----------------------------------------------------------------
 
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/2025-2412-UMC-EA-s023-s030/raw_count_tables/non_poisson_corrected/UMC-EA-s029-raw")
  
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
  sum(rawData==0) #13497649
  
  # Calculate gene counts -> so in how many cells a gene is counted
  geneCounts <- apply(rawData,1,function (x){sum(x>0)}) # 1 indicates rows
  min(geneCounts) # 0
  max(geneCounts) # 373
  head(geneCounts)
  
  mean(geneCounts) #7.221852
  
  length(geneCounts[geneCounts<3]) #26465
  length(geneCounts[geneCounts==0]) #22398
  
  # Calculate cell counts -> so how many genes are expressed per cell
  cellCounts<-apply(rawData,2,function (x){sum(x>0)}) #2 indicates columns
  
  min(cellCounts) # 23
  max(cellCounts) # 1293
  
  mean(cellCounts) #702.9973
  
  length(cellCounts[cellCounts<300]) #44
  length(cellCounts[cellCounts == 0]) # 0
  
  # Remove the values created to get an overview
  rm(geneCounts)
  rm(cellCounts)
  
### 4. Read Flow data -------------------------------------------------------------------------------------------------------------
  setwd("T:/cff-data/WKZ/Group-vanWijk/Elise/Sort experiments/HINT/241129_HINT148_blood_scRNAseq_sort/")
  ## 
  EA029_flow_x <-read.FCS('HINT148_blood_INX_Tube_001_013_compensated.fcs', truncate_max_range = FALSE)
  View(EA029_flow_x@exprs)
  
  EA029_flow <- getIndexSort(EA029_flow_x)
  ## see barcode file: scRNAseq data are ordered from A1->A24 and then B1->B24
  ## but cells were sorted (+index sort data) from A1->A24 and then B24->B1
  ## xloc = rows yloc = columns
  EA029_flow <- EA029_flow[order(EA029_flow$XLoc,EA029_flow$YLoc),]
  hist(EA029_flow$YLoc)
  hist(EA029_flow$XLoc, breaks=seq(0,15,1))
 ## many rowss missing
  dummy_df <- EA029_flow
  dummy_df[1:255,1:26] <- NA
  EA029_flow <- rbind(EA029_flow[1:54,],dummy_df[1,],EA029_flow[55:56,],dummy_df[1,],EA029_flow[57:61,],
                      dummy_df[1,],EA029_flow[62:99,], dummy_df[1,],EA029_flow[100:146,],dummy_df[1,],
                      EA029_flow[147:190,],dummy_df[1,],  EA029_flow[191:194,],dummy_df[1,],
                      EA029_flow[195:214,],dummy_df[1,],EA029_flow[215:272,],dummy_df[1,],
                      EA029_flow[273:362,],dummy_df[1,],EA029_flow[363:366,])
  
  rownames(EA029_flow)  <- colnames(count_matrix)
  colnames(EA029_flow)[8:23] <- EA029_flow_x@parameters@data$desc[7:22]
  colnames(EA029_flow)[is.na(colnames(EA029_flow))] <- 'empty'
  ##check compensation
  plot(EA029_flow$CD4~EA029_flow$CD31)
  plot(EA029_flow$CD95~EA029_flow$CD8)
  plot(EA029_flow$CCR7~EA029_flow$CD127)
  
  View(EA029_flow)
  dim(EA029_flow)
  
  
### 5. Seurat object with flow data -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data. remove genes that are expressed in less than 3 cells
  EA029 <- CreateSeuratObject(counts = rawData, project = "EA029", min.cells=3)
  EA029
  #add flowdata as metadata
  EA029@meta.data$flowCD31 <- EA029_flow$CD31
  EA029@meta.data$flowCD25 <- EA029_flow$CD25
  EA029@meta.data$flowCD45RA <- EA029_flow$CD45RA
  EA029@meta.data$flowCD27 <- EA029_flow$CD27
  EA029@meta.data$flowCCR7 <- EA029_flow$CCR7
  EA029@meta.data$flowCD127 <- EA029_flow$CD127
  EA029@meta.data$flowCD69 <- EA029_flow$CD69
  EA029 

### 6. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/QC")

  ## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  EA029 <- PercentageFeatureSet(EA029, pattern = "^MT\\.", col.name = 'percent.mt')
  EA029 <- PercentageFeatureSet(EA029, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(count_matrix)[grep('MT', rownames(count_matrix))]) #somewhere in name
  length(rownames(count_matrix)[grep('^MT\\.', rownames(count_matrix))]) #start of name
  
  # Show QC metrics for the first 5 cells
  head(EA029@meta.data, 5)
  
  #save(EA029, file = "EA029.Robj")
  
  ## Plot histograms of QC data pre-filtering
  toPlot <- EA029@meta.data
  
  # percent.mito
  png(file="EA029_pre_percMito.png", width=850)
  par(mfrow=c(1,2))
  tmp <- toPlot[order(toPlot$percent.mt),]
  hist(tmp$percent.mt, breaks=30)
  barplot(tmp$percent.mt)
  dev.off()
  
  # nFeature_RNA
  png(file="EA029_pre_nFeature_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nFeature_RNA),]
  hist(tmp$nFeature_RNA, breaks=30)
  barplot(tmp$nFeature_RNA)
  dev.off()
  
  # nCount_RNA
  png(file="EA029_pre_nCount_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nCount_RNA),]
  hist(tmp$nCount_RNA, breaks=30)
  barplot(tmp$nCount_RNA)
  dev.off()
  
  # Visualize QC metrics as a violin plot
  pdf(file="EA029_beforeQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA029, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="EA029_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(EA029, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(EA029, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(EA029, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  plot1 + plot2 + plot3
  dev.off()
  
  #With filtering cut-offs 
  pdf(file="EA029_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
  plot1 <- ggplot(EA029@meta.data, aes(EA029$nCount_RNA, EA029$nFeature_RNA))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 3350)+ annotate("text",x=2300,y=1,label=c("3350"),hjust=0, size=2.8)
    plot2 <- ggplot(EA029@meta.data, aes(EA029$nCount_RNA, EA029$percent.mt))+geom_point(size=.9)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 350)+ annotate("text",x=380,y=1,label=c("350"),hjust=0, size=2.8)+
    geom_hline(yintercept = 16)+ annotate("text",x=700,y=16.5,label=c("16"),vjust=0, size=2.8)
  plot3 <- ggplot(EA029@meta.data, aes(EA029$nFeature_RNA, EA029$percent.mt))+geom_point(size=.9)+theme_bw()+ 
    geom_vline(xintercept = 300)+annotate("text",x=320,y=1,label=c("300"),hjust=0, size=2.8)+
    geom_hline(yintercept = 16)+ annotate("text",x=500,y=16.5,label=c("16"),vjust=0, size=2.8)
  plot4 <- ggplot(EA029@meta.data, aes(EA029$nCount_RNA, EA029$percent.ribo))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 350)+ annotate("text",x=380,y=1,label=c("350"),hjust=0, size=2.8)
  plot5 <- ggplot(EA029@meta.data, aes(EA029$percent.mt, EA029$percent.ribo))+geom_point(size=.9)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(EA029, expression = nFeature_RNA > 300 & nFeature_RNA < 3350 &nCount_RNA > 350 & percent.mt < 16)
  length(selected) #How many cells are left
  
  # Filter cells
  EA029 <- subset(EA029, subset = nFeature_RNA > 300 & nFeature_RNA < 3350 &nCount_RNA > 350 & percent.mt < 16)
  
  # Visualize QC metrics post-filtering
  pdf(file="EA029_afterQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA029, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
### 7. SCTransform and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  EA029 <- NormalizeData(EA029, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  EA029 <- CellCycleScoring(object = EA029, s.features = s.genes, g2m.features = g2m.genes, 
                            set.ident = FALSE)
  head(x = EA029@meta.data)
  
  
  ## SCTransform is a code for UMI corrected matrices
  #optional (new version, don't know if it's default already, haven't tried yet): vst.flavor='v2'
  ## I think I don't want to regress out cell-cycle, as you don't except many cells to be dividing, and if so, that's actually interesting info
  EA029_mt <- SCTransform(object = EA029, verbose = TRUE, 
                       vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
  
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq")
  save(EA029_mt, file = 'EA029SCT_pctmt_nocellcylce.Robj')
  
  #with percent.ribo
  EA029_mtrb <- SCTransform(object = EA029, verbose = TRUE, 
                       vars.to.regress = c("percent.mt", "percent.ribo"), return.only.var.genes = FALSE)
  save(EA029_mtrb, file = 'EA029SCT_pctmtribo_nocellcylce.Robj')
  
  #with cell cycle
  EA029_mtcc <- SCTransform(object = EA029, verbose = TRUE, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), return.only.var.genes = FALSE)
  save(EA029_mtcc, file = 'EA029SCT_pctmtcc_cellcylce.Robj')
  
### 8. Optional: First check of clustering --------------------------------------------------------------------------------------------
 ## weird that the KRT genes are not there now (specific to that sequencing run?)
  EA029_mt <- RunPCA(EA029_mt, verbose = TRUE)
  EA029_mt <- RunUMAP(EA029_mt, dims = 1:30)
  
  
  ## 
  EA029_mt <- FindNeighbors(object = EA029_mt, reduction = "pca", dims = 1:30)
  ## not more clusters than with naive plate
  EA029_mt <- FindClusters(object = EA029_mt, resolution = .9)
  
  DimPlot(object = EA029_mt, reduction = "umap", pt.size=1.2)

  FeaturePlot(EA029_mt, 'flowCD31', pt.size=2)
  FeaturePlot(EA029_mt, 'flowCD127', pt.size=2)
  VlnPlot(EA029_mt, 'flowCD31')
  
  ## better: less MT and RP genes
  EA029_mtrb <- RunPCA(EA029_mtrb, verbose = TRUE)
  EA029_mtrb <- RunUMAP(EA029_mtrb, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA029_mtrb <- FindNeighbors(object = EA029_mtrb, reduction = "pca", dims = 1:30)
  ##
  EA029_mtrb <- FindClusters(object = EA029_mtrb, resolution = .9)
  
  DimPlot(object = EA029_mtrb, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA029_mtrb, 'flowCD31', pt.size=2)
  FeaturePlot(EA029_mtrb, 'flowCD127', pt.size=2)
  VlnPlot(EA029_mtrb, 'flowCD31')
  
  ## similar to mt
  EA029_mtcc <- RunPCA(EA029_mtcc, veccose = TRUE)
  EA029_mtcc <- RunUMAP(EA029_mtcc, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA029_mtcc <- FindNeighbors(object = EA029_mtcc, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA029_mtcc <- FindClusters(object = EA029_mtcc, resolution = .95)
  
  DimPlot(object = EA029_mtcc, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA029_mtcc, 'flowCD31', pt.size=2)
  FeaturePlot(EA029_mtcc, 'flowCD127', pt.size=2)
  VlnPlot(EA029_mtcc, 'flowCD31')