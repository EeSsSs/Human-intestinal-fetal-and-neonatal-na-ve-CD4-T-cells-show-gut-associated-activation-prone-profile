#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
  library(Seurat)
  library(dplyr)
  
  library(sctransform)
  library(ggplot2)
  library(flowCore)
  library(Matrix)

### 1. Setting working directory and loading data----------------------------------------------------------------
 
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/2022-1754-umc-ea-s012-s022/raw_count_tables/non_poisson_corrected/UMC-EA-s014-raw")

  
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
  ##sorted until H1
  count_matrix <- count_matrix[,-c(193:384)]
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
  dim(rawData) #36601  192
  
  # Number (sum) of zero's in the dataset
  sum(rawData==0) #6924681
  
  # Calculate gene counts -> so in how many cells a gene is counted
  geneCounts <- apply(rawData,1,function (x){sum(x>0)}) # 1 indicates rows
  min(geneCounts) # 0
  max(geneCounts) # 192
  head(geneCounts)
  
  mean(geneCounts) #2.806235
  
  length(geneCounts[geneCounts<3]) #28954
  length(geneCounts[geneCounts==0]) #23996
  
  # Calculate cell counts -> so how many genes are expressed per cell
  cellCounts<-apply(rawData,2,function (x){sum(x>0)}) #2 indicates columns
  min(cellCounts) # 56
  max(cellCounts) # 1295
  
  mean(cellCounts) #534.9531
  
  length(cellCounts[cellCounts<300]) #60
  length(cellCounts[cellCounts == 0]) # 0
  
  # Remove the values created to get an overview
  rm(geneCounts)
  rm(cellCounts)
  
### 4. Read Flow data -------------------------------------------------------------------------------------------------------------
  setwd("T:/cff-data/WKZ/Group-vanWijk/Elise/Sort experiments/HINT/230704_HINT130_gut_scRNAseq_sort")
  ## apparently there are two index files, one with 241 cells and one with 29 cells
  EA014_flow_x <-read.FCS('HINT130_gut_INX_Tube_009_021_compensated.fcs')
  
  View(EA014_flow_x@exprs)
  
  EA014_flow <- getIndexSort(EA014_flow_x)
  
  ## see barcode file: scRNAseq data are ordered from A1->A24 and then B1->B24
  ## but cells were sorted (+index sort data) from A1->A24 and then B24->B1
  ## xloc = rows yloc = columns
  EA014_flow <- EA014_flow[order(EA014_flow$XLoc,EA014_flow$YLoc),]
  
  rownames(EA014_flow)  <- colnames(count_matrix)
  colnames(EA014_flow)[8:23] <- EA014_flow_x@parameters@data$desc[7:22]
  colnames(EA014_flow)[is.na(colnames(EA014_flow))] <- 'empty'
  ##check compensation
  plot(EA014_flow$CD4~EA014_flow$CD31)
  plot(EA014_flow$CD95~EA014_flow$CD8)
  
  View(EA014_flow)
  dim(EA014_flow)
  
  
### 5. Seurat object with flow data -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data. remove genes that are expressed in less than 3 cells
  EA014 <- CreateSeuratObject(counts = rawData, project = "EA014", min.cells=3)
  EA014
  #add flowdata as metadata
  EA014@meta.data$flowCD31 <- EA014_flow$CD31
  EA014@meta.data$flowCD25 <- EA014_flow$CD25
  EA014@meta.data$flowCD45RA <- EA014_flow$CD45RA
  EA014@meta.data$flowCD27 <- EA014_flow$CD27
  EA014@meta.data$flowCCR7 <- EA014_flow$CCR7
  EA014@meta.data$flowCD127 <- EA014_flow$CD127
  EA014 

### 6. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/QC")

  ## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  EA014 <- PercentageFeatureSet(EA014, pattern = "^MT\\.", col.name = 'percent.mt')
  EA014 <- PercentageFeatureSet(EA014, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(count_matrix)[grep('MT', rownames(count_matrix))]) #somewhere in name
  length(rownames(count_matrix)[grep('^MT\\.', rownames(count_matrix))]) #start of name
  
  # Show QC metrics for the first 5 cells
  head(EA014@meta.data, 5)
  
  #save(EA014, file = "EA014.Robj")
  
  ## Plot histograms of QC data pre-filtering
  toPlot <- EA014@meta.data
  
  # percent.mito
  png(file="EA014_pre_percMito.png", width=850)
  par(mfrow=c(1,2))
  tmp <- toPlot[order(toPlot$percent.mt),]
  hist(tmp$percent.mt, breaks=30)
  barplot(tmp$percent.mt)
  dev.off()
  
  # nFeature_RNA
  png(file="EA014_pre_nFeature_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nFeature_RNA),]
  hist(tmp$nFeature_RNA, breaks=30)
  barplot(tmp$nFeature_RNA)
  dev.off()
  
  # nCount_RNA
  png(file="EA014_pre_nCount_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nCount_RNA),]
  hist(tmp$nCount_RNA, breaks=30)
  barplot(tmp$nCount_RNA)
  dev.off()
  
  # Visualize QC metrics as a violin plot
  pdf(file="EA014_beforeQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA014, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="EA014_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(EA014, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(EA014, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(EA014, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  plot1 + plot2 + plot3
  dev.off()
  
  #With filtering cut-offs
  pdf(file="EA014_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
  plot1 <- ggplot(EA014@meta.data, aes(EA014$nCount_RNA, EA014$nFeature_RNA))+geom_point(size=.9)+theme_bw()
  plot2 <- ggplot(EA014@meta.data, aes(EA014$nCount_RNA, EA014$percent.mt))+geom_point(size=.9)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 400)+ annotate("text",x=430,y=1,label=c("400"),hjust=0, size=2.8)+
    geom_hline(yintercept = 19)+ annotate("text",x=700,y=19.5,label=c("19"),vjust=0, size=2.8)
  plot3 <- ggplot(EA014@meta.data, aes(EA014$nFeature_RNA, EA014$percent.mt))+geom_point(size=.9)+theme_bw()+ 
    geom_vline(xintercept = 350)+annotate("text",x=370,y=1,label=c("350"),hjust=0, size=2.8)+
    geom_hline(yintercept = 19)+ annotate("text",x=500,y=19.5,label=c("19"),vjust=0, size=2.8)
  plot4 <- ggplot(EA014@meta.data, aes(EA014$nCount_RNA, EA014$percent.ribo))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 400)+ annotate("text",x=430,y=1,label=c("400"),hjust=0, size=2.8)
  plot5 <- ggplot(EA014@meta.data, aes(EA014$percent.mt, EA014$percent.ribo))+geom_point(size=.9)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(EA014, expression = nFeature_RNA > 350 & nCount_RNA < 5000 & nCount_RNA > 400 & percent.mt < 19)
  length(selected) #How many cells are left
  
  # Filter cells
  EA014 <- subset(EA014, subset = nFeature_RNA > 350 & nCount_RNA < 5000 & nCount_RNA > 400 & percent.mt < 19)
  
  # Visualize QC metrics post-filtering
  pdf(file="EA014_afterQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA014, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
### 7. SCTransform and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  EA014 <- NormalizeData(EA014, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  EA014 <- CellCycleScoring(object = EA014, s.features = s.genes, g2m.features = g2m.genes, 
                            set.ident = FALSE)
  head(x = EA014@meta.data)
  
  
  ## SCTransform is a code for UMI corrected matrices
  #optional (new version, don't know if it's default already, haven't tried yet): vst.flavor='v2'
  ## I think I don't want to regress out cell-cycle, as you don't except many cells to be dividing, and if so, that's actually interesting info
  EA014_mt <- SCTransform(object = EA014, verbose = TRUE, 
                       vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
  
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
  save(EA014_mt, file = 'EA014_SCT_pctmt_nocellcylce.Robj')
  
  #with percent.ribo
  EA014_mtrb <- SCTransform(object = EA014, verbose = TRUE, 
                       vars.to.regress = c("percent.mt", "percent.ribo"), return.only.var.genes = FALSE)
  save(EA014_mtrb, file = 'EA014SCT_pctmtribo_nocellcylce.Robj')
  
  #with cell cycle
  EA014_mtcc <- SCTransform(object = EA014, verbose = TRUE, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), return.only.var.genes = FALSE)
  save(EA014_mtcc, file = 'EA014SCT_pctmtcc_cellcylce.Robj')
  
### 8. Optional: First check of clustering --------------------------------------------------------------------------------------------
  EA014_mt <- RunPCA(EA014_mt, verbose = TRUE)
  EA014_mt <- RunUMAP(EA014_mt, dims = 1:30)
  
  
  ## set back defaultassay to integrated if reclustering
  EA014_mt <- FindNeighbors(object = EA014_mt, reduction = "pca", dims = 1:30)
  EA014_mt <- FindClusters(object = EA014_mt, resolution = 0.9)
  
  DimPlot(object = EA014_mt, reduction = "umap", pt.size=1.2)

  FeaturePlot(EA014_mt, 'flowCD31', pt.size=2)
  FeaturePlot(EA014_mt, 'flowCD127', pt.size=2)
  VlnPlot(EA014_mt, 'flowCD31')
  
  ## better: less MT and RP genes
  EA014_mtrb <- RunPCA(EA014_mtrb, verbose = TRUE)
  EA014_mtrb <- RunUMAP(EA014_mtrb, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA014_mtrb <- FindNeighbors(object = EA014_mtrb, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA014_mtrb <- FindClusters(object = EA014_mtrb, resolution = .95)
  
  DimPlot(object = EA014_mtrb, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA014_mtrb, 'flowCD31', pt.size=2)
  FeaturePlot(EA014_mtrb, 'flowCD127', pt.size=2)
  VlnPlot(EA014_mtrb, 'flowCD31')
  
  ## similar to mt
  EA014_mtcc <- RunPCA(EA014_mtcc, veccose = TRUE)
  EA014_mtcc <- RunUMAP(EA014_mtcc, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA014_mtcc <- FindNeighbors(object = EA014_mtcc, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA014_mtcc <- FindClusters(object = EA014_mtcc, resolution = .95)
  
  DimPlot(object = EA014_mtcc, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA014_mtcc, 'flowCD31', pt.size=2)
  FeaturePlot(EA014_mtcc, 'flowCD127', pt.size=2)
  VlnPlot(EA014_mtcc, 'flowCD31')