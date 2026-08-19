#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
  library(Seurat)
  library(dplyr)
  
  library(sctransform)
  library(ggplot2)
  library(flowCore)
  library(Matrix)

### 1. Setting working directory and loading data----------------------------------------------------------------
 
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/umc-ea-s013-s017-20231025/raw_count_tables/non_poisson_corrected/UMC-EA-s013-raw/")
  
  # Read in matrix.mtx
  counts <- readMM("matrix.mtx")
  counts
  
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
  count_matrix <- count_matrix[,-c(357:360, 381:384)]
  dim(count_matrix) #Should leave 376 columns (384-8 for the spike ins)
  
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
  sum(rawData==0) #13522071
  
  # Calculate gene counts -> so in how many cells a gene is counted
  geneCounts <- apply(rawData,1,function (x){sum(x>0)}) # 1 indicates rows
  min(geneCounts) # 0
  max(geneCounts) # 376
  head(geneCounts)
  
  mean(geneCounts) #6.554602
  
  length(geneCounts[geneCounts<3]) #26493
  length(geneCounts[geneCounts==0]) #21741
  
  # Calculate cell counts -> so how many genes are expressed per cell
  cellCounts<-apply(rawData,2,function (x){sum(x>0)}) #2 indicates columns
  min(cellCounts) # 26
  max(cellCounts) # 1709
  
  mean(cellCounts) #638.0452
  
  length(cellCounts[cellCounts<300]) #61
  length(cellCounts[cellCounts == 0]) # 0
  
  # Remove the values created to get an overview
  rm(geneCounts)
  rm(cellCounts)
  
### 4. Read Flow data -------------------------------------------------------------------------------------------------------------
  setwd("T:/cff-data/WKZ/Group-vanWijk/Elise/Sort experiments/HINT/230510_HINT129_blood_scRNAseq_sort")
  EA013_flow_x <-read.FCS('HINT129_blood_INX_Tube_007_019_compensated.fcs')
  View(EA013_flow_x@exprs)
  
  EA013_flow <- getIndexSort(EA013_flow_x)
  
  ## see barcode file: scRNAseq data are ordered from A1->A24 and then B1->B24
  ## but cells were sorted (+index sort data) from A1->A24 and then B24->B1
  ## xloc = rows yloc = columns
  EA013_flow <- EA013_flow[order(EA013_flow$XLoc,EA013_flow$YLoc),]
  
  ## only 375 rows: second well of third row is missing -> create empty row
  EA013_flow <- add_row(EA013_flow, XLoc=2, YLoc=1, .after = 49)
  
  rownames(EA013_flow)  <- colnames(count_matrix)
  colnames(EA013_flow)[8:23] <- EA013_flow_x@parameters@data$desc[7:22]
  colnames(EA013_flow)[is.na(colnames(EA013_flow))] <- 'empty'
  ##check compensation
  plot(EA013_flow$CD4~EA013_flow$CD31)
  plot(EA013_flow$CD95~EA013_flow$CD8)
  
  View(EA013_flow)
  dim(EA013_flow)
  
  
### 5. Seurat object with flow data -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data. remove genes that are expressed in less than 3 cells
  EA013 <- CreateSeuratObject(counts = rawData, project = "EA013", min.cells=3)
  EA013
  #add flowdata as metadata
  EA013@meta.data$flowCD31 <- EA013_flow$CD31
  EA013@meta.data$flowCD25 <- EA013_flow$CD25
  EA013@meta.data$flowCD45RA <- EA013_flow$CD45RA
  EA013@meta.data$flowCD27 <- EA013_flow$CD27
  EA013@meta.data$flowCCR7 <- EA013_flow$CCR7
  EA013@meta.data$flowCD127 <- EA013_flow$CD127
  EA013 

### 6. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq/QC")

  ## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  EA013 <- PercentageFeatureSet(EA013, pattern = "^MT\\.", col.name = 'percent.mt')
  EA013 <- PercentageFeatureSet(EA013, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(count_matrix)[grep('MT', rownames(count_matrix))]) #somewhere in name
  length(rownames(count_matrix)[grep('^MT\\.', rownames(count_matrix))]) #start of name
  
  # Show QC metrics for the first 5 cells
  head(EA013@meta.data, 5)
  
  save(EA013, file = "EA013.Robj")
  
  ## Plot histograms of QC data pre-filtering
  toPlot <- EA013@meta.data
  
  # percent.mito
  png(file="EA013_pre_percMito.png", width=850)
  par(mfrow=c(1,2))
  tmp <- toPlot[order(toPlot$percent.mt),]
  hist(tmp$percent.mt, breaks=30)
  barplot(tmp$percent.mt)
  dev.off()
  
  # nFeature_RNA
  png(file="EA013_pre_nFeature_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nFeature_RNA),]
  hist(tmp$nFeature_RNA, breaks=30)
  barplot(tmp$nFeature_RNA)
  dev.off()
  
  # nCount_RNA
  png(file="EA013_pre_nCount_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nCount_RNA),]
  hist(tmp$nCount_RNA, breaks=30)
  barplot(tmp$nCount_RNA)
  dev.off()
  
  # Visualize QC metrics as a violin plot
  pdf(file="EA013_beforeQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA013, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="EA013_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(EA013, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(EA013, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(EA013, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  plot1 + plot2 + plot3
  dev.off()
  
  #With filtering cut-offs
  pdf(file="EA013_correlations QC_ggplot_withfilterlingablines.pdf", width = 24, height = 12)
  plot1 <- ggplot(EA013@meta.data, aes(EA013$nCount_RNA, EA013$nFeature_RNA))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 3000)+ annotate("text",x=3100,y=1,label=c("3000"),hjust=0, size=2.8)
  plot2 <- ggplot(EA013@meta.data, aes(EA013$nCount_RNA, EA013$percent.mt))+geom_point(size=.9)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 250)+ annotate("text",x=280,y=1,label=c("250"),hjust=0, size=2.8)+
    geom_hline(yintercept = 18)+ annotate("text",x=700,y=18.5,label=c("18"),vjust=0, size=2.8)
  plot3 <- ggplot(EA013@meta.data, aes(EA013$nFeature_RNA, EA013$percent.mt))+geom_point(size=.9)+theme_bw()+ 
    geom_vline(xintercept = 200)+annotate("text",x=220,y=1,label=c("200"),hjust=0, size=2.8)+
    geom_hline(yintercept = 18)+ annotate("text",x=500,y=18.5,label=c("18"),vjust=0, size=2.8)
  plot4 <- ggplot(EA013@meta.data, aes(EA013$nCount_RNA, EA013$percent.ribo))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 250)+ annotate("text",x=280,y=1,label=c("250"),hjust=0, size=2.8)
  plot5 <- ggplot(EA013@meta.data, aes(EA013$percent.mt, EA013$percent.ribo))+geom_point(size=.9)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(EA013, expression = nFeature_RNA > 200 & nCount_RNA < 3000 & nCount_RNA > 250 & percent.mt < 18)
  length(selected) #How many cells are left
  
  # Filter cells
  EA013 <- subset(EA013, subset = nFeature_RNA > 200 & nCount_RNA < 3000 & nCount_RNA > 250 & percent.mt < 18)
  
  # Visualize QC metrics post-filtering
  pdf(file="EA013_afterQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA013, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
### 7. SCTransform and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  EA013 <- NormalizeData(EA013, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  EA013 <- CellCycleScoring(object = EA013, s.features = s.genes, g2m.features = g2m.genes, 
                            set.ident = FALSE)
  head(x = EA013@meta.data)
  
  
  ## SCTransform is a code for UMI corrected matrices
  #optional (new version, don't know if it's default already, haven't tried yet): vst.flavor='v2'
  ## I think I don't want to regress out cell-cycle, as you don't except many cells to be dividing, and if so, that's actually interesting info
  EA013_mt <- SCTransform(object = EA013, verbose = TRUE, 
                       vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
  
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/HC SORT-seq")
  save(EA013_mt, file = 'EA013_SCT_pctmt_nocellcylce.Robj')
  
  #with percent.ribo
  EA013_mtrb <- SCTransform(object = EA013, verbose = TRUE, 
                       vars.to.regress = c("percent.mt", "percent.ribo"), return.only.var.genes = FALSE)
  save(EA013_mtrb, file = 'EA013SCT_ML_pctmtribo_nocellcylce.Robj')
  
### 8. Optional: First check of clustering --------------------------------------------------------------------------------------------
  ## KRT genes in PCA?????
  EA013_mt <- RunPCA(EA013_mt, verbose = TRUE)
  EA013_mt <- RunUMAP(EA013_mt, dims = 1:30)
  
  
  ## set back defaultassay to integrated if reclustering
  EA013_mt <- FindNeighbors(object = EA013_mt, reduction = "pca", dims = 1:30)
  EA013_mt <- FindClusters(object = EA013_mt, resolution = 0.7)
  
  DimPlot(object = EA013_mt, reduction = "umap", pt.size=1.2)

  FeaturePlot(EA013_mt, 'flowCD31', pt.size=2)
  FeaturePlot(EA013_mt, 'flowCD127', pt.size=2)
  VlnPlot(EA013_mt, 'flowCD31')
  
  ## KRT genes in PCA?????
  EA013_mtrb <- RunPCA(EA013_mtrb, verbose = TRUE)
  EA013_mtrb <- RunUMAP(EA013_mtrb, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA013_mtrb <- FindNeighbors(object = EA013_mtrb, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters
  EA013_mtrb <- FindClusters(object = EA013_mtrb, resolution = 0.9)
  
  DimPlot(object = EA013_mtrb, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA013_mtrb, 'flowCD31', pt.size=2)
  FeaturePlot(EA013_mtrb, 'flowCD127', pt.size=2)
  VlnPlot(EA013_mtrb, 'flowCD31')
  