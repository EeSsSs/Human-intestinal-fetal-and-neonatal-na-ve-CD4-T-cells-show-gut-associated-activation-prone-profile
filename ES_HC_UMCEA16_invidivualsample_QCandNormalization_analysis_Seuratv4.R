#### 0. Loading required packages to run this R script for Seurat (single cell RNA-seq) ####
  library(Seurat)
  library(dplyr)
  
  library(sctransform)
  library(ggplot2)
  library(flowCore)
  library(Matrix)

### 1. Setting working directory and loading data----------------------------------------------------------------
 
setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/2022-1754-umc-ea-s012-s022/raw_count_tables/non_poisson_corrected/UMC-EA-s016-raw")

  
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
  ##sorted until H24 written on lab journal, but based on analysis report of seurat I would say J24? (>500 UMIs/well after H24)
  ##see below index sort data -> also 217 indicating that I stopped later than H24
  count_matrix <- count_matrix[,-c(217:239, 241:384)]
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
  dim(rawData) #36601  217
  
  # Number (sum) of zero's in the dataset
  sum(rawData==0) #7822985
  
  # Calculate gene counts -> so in how many cells a gene is counted
  geneCounts <- apply(rawData,1,function (x){sum(x>0)}) # 1 indicates rows
  min(geneCounts) # 0
  max(geneCounts) # 217
  head(geneCounts)
  
  mean(geneCounts) #3.26308
  
  length(geneCounts[geneCounts<3]) #28219
  length(geneCounts[geneCounts==0]) #22627
  
  # Calculate cell counts -> so how many genes are expressed per cell
  cellCounts<-apply(rawData,2,function (x){sum(x>0)}) #2 indicates columns
  min(cellCounts) # 54
  max(cellCounts) # 2376
  
  mean(cellCounts) #550.3779
  
  length(cellCounts[cellCounts<300]) #65
  length(cellCounts[cellCounts == 0]) # 0
  
  # Remove the values created to get an overview
  rm(geneCounts)
  rm(cellCounts)
  
### 4. Read Flow data -------------------------------------------------------------------------------------------------------------
  setwd("T:/cff-data/WKZ/Group-vanWijk/Elise/Sort experiments/HINT/230711_HINT131_gut_scRNAseq_sort")
  ## apparently there are two index files, one with 241 cells and one with 29 cells
  EA016_flow_x <-read.FCS('HINT131_gut_INX_Tube_010_022_compensated.fcs')
  
  View(EA016_flow_x@exprs)
  
  EA016_flow <- getIndexSort(EA016_flow_x)
  
  ## see barcode file: scRNAseq data are ordered from A1->A24 and then B1->B24
  ## but cells were sorted (+index sort data) from A1->A24 and then B24->B1
  ## xloc = rows yloc = columns
  EA016_flow <- EA016_flow[order(EA016_flow$XLoc,EA016_flow$YLoc),]
  
  rownames(EA016_flow)  <- colnames(count_matrix)
  colnames(EA016_flow)[8:23] <- EA016_flow_x@parameters@data$desc[7:22]
  colnames(EA016_flow)[is.na(colnames(EA016_flow))] <- 'empty'
  ##check compensation
  plot(EA016_flow$CD4~EA016_flow$CD31)
  plot(EA016_flow$CD95~EA016_flow$CD8)
  
  View(EA016_flow)
  dim(EA016_flow)
  
  
### 5. Seurat object with flow data -------------------------------------------------------------------------------------
  ## Initialize the Seurat object with the raw (non-normalized) data. remove genes that are expressed in less than 3 cells
  EA016 <- CreateSeuratObject(counts = rawData, project = "EA016", min.cells=3)
  EA016
  #add flowdata as metadata
  EA016@meta.data$flowCD31 <- EA016_flow$CD31
  EA016@meta.data$flowCD25 <- EA016_flow$CD25
  EA016@meta.data$flowCD45RA <- EA016_flow$CD45RA
  EA016@meta.data$flowCD27 <- EA016_flow$CD27
  EA016@meta.data$flowCCR7 <- EA016_flow$CCR7
  EA016@meta.data$flowCD127 <- EA016_flow$CD127
  EA016 

### 6. Quality control in Seurat--------------------------------------------------------------------------------------
  ## Change working directory where you want to save the QC images per plate
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq/QC")

  ## Pre-processing workflow
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  EA016 <- PercentageFeatureSet(EA016, pattern = "^MT\\.", col.name = 'percent.mt')
  EA016 <- PercentageFeatureSet(EA016, pattern = '^RP', col.name = 'percent.ribo')
  
  length(rownames(count_matrix)[grep('MT', rownames(count_matrix))]) #somewhere in name
  length(rownames(count_matrix)[grep('^MT\\.', rownames(count_matrix))]) #start of name
  
  # Show QC metrics for the first 5 cells
  head(EA016@meta.data, 5)
  
  #save(EA016, file = "EA016.Robj")
  
  ## Plot histograms of QC data pre-filtering
  toPlot <- EA016@meta.data
  
  # percent.mito
  png(file="EA016_pre_percMito.png", width=850)
  par(mfrow=c(1,2))
  tmp <- toPlot[order(toPlot$percent.mt),]
  hist(tmp$percent.mt, breaks=30)
  barplot(tmp$percent.mt)
  dev.off()
  
  # nFeature_RNA
  png(file="EA016_pre_nFeature_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nFeature_RNA),]
  hist(tmp$nFeature_RNA, breaks=30)
  barplot(tmp$nFeature_RNA)
  dev.off()
  
  # nCount_RNA
  png(file="EA016_pre_nCount_RNA.png", width=850)
  par(mfrow=c(1,2))
  tmp<-toPlot[order(toPlot$nCount_RNA),]
  hist(tmp$nCount_RNA, breaks=30)
  barplot(tmp$nCount_RNA)
  dev.off()
  
  # Visualize QC metrics as a violin plot
  pdf(file="EA016_beforeQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA016, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf(file="EA016_correlationsQC.pdf", width = 24, height = 8)
  plot1 <- FeatureScatter(EA016, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(EA016, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(EA016, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  plot1 + plot2 + plot3
  dev.off()
  
  #With filtering cut-offs --> difficult!
  pdf(file="EA016_correlations QC_ggplot_withfilterlingablines_revisit.pdf", width = 24, height = 12)
  plot1 <- ggplot(EA016@meta.data, aes(EA016$nCount_RNA, EA016$nFeature_RNA))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 3000)+ annotate("text",x=3100,y=1,label=c("3000"),hjust=0, size=2.8)
  plot2 <- ggplot(EA016@meta.data, aes(EA016$nCount_RNA, EA016$percent.mt))+geom_point(size=.9)+
    theme(axis.text.x = element_text(size=0.3))+theme_bw()+
    geom_vline(xintercept = 300)+ annotate("text",x=330,y=1,label=c("300"),hjust=0, size=2.8)+
    geom_hline(yintercept = 19)+ annotate("text",x=700,y=19.5,label=c("19"),vjust=0, size=2.8)
  plot3 <- ggplot(EA016@meta.data, aes(EA016$nFeature_RNA, EA016$percent.mt))+geom_point(size=.9)+theme_bw()+ 
    geom_vline(xintercept = 250)+annotate("text",x=270,y=1,label=c("250"),hjust=0, size=2.8)+
    geom_hline(yintercept = 19)+ annotate("text",x=500,y=19.5,label=c("19"),vjust=0, size=2.8)
  plot4 <- ggplot(EA016@meta.data, aes(EA016$nCount_RNA, EA016$percent.ribo))+geom_point(size=.9)+theme_bw()+
    geom_vline(xintercept = 300)+ annotate("text",x=330,y=1,label=c("300"),hjust=0, size=2.8)
  plot5 <- ggplot(EA016@meta.data, aes(EA016$percent.mt, EA016$percent.ribo))+geom_point(size=.9)+theme_bw()
  plot1 + plot2 + plot3+plot4+plot5
  dev.off()
  
  ## Filtering based on QC metrics
  # Check filtering to be set
  selected <- WhichCells(EA016, expression = nFeature_RNA > 250 & nCount_RNA < 3000 & nCount_RNA > 300 & percent.mt < 19)
  length(selected) #How many cells are left
  
  # Filter cells
  EA016 <- subset(EA016, subset = nFeature_RNA > 250 & nCount_RNA < 3000 & nCount_RNA > 300 & percent.mt < 19)
  
  # Visualize QC metrics post-filtering
  pdf(file="EA016_afterQC_Vlns.pdf", width = 16, height = 8)
  VlnPlot(EA016, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  dev.off()
  
### 7. SCTransform and cell cycle scoring---------------------------------------------------------------------------------------------------------
  
  ## Full regression for cell cycle
  # Normalize to enable cell cycle scoring (new normalization (SCT) will follow on count data)
  EA016 <- NormalizeData(EA016, normalization.method = "LogNormalize", scale.factor = 10000)
  # Segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Assign cell cycle scores to the cells 
  EA016 <- CellCycleScoring(object = EA016, s.features = s.genes, g2m.features = g2m.genes, 
                            set.ident = FALSE)
  head(x = EA016@meta.data)
  
  
  ## SCTransform is a code for UMI corrected matrices
  #optional (new version, don't know if it's default already, haven't tried yet): vst.flavor='v2'
  ## I think I don't want to regress out cell-cycle, as you don't except many cells to be dividing, and if so, that's actually interesting info
  EA016_mt <- SCTransform(object = EA016, verbose = TRUE, 
                       vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
  
  setwd("T:/lti/_Clinical/Group-van-Wijk/Elise/scRNAseq/Gut SORT-seq")
  save(EA016_mt, file = 'EA016_SCT_pctmt_nocellcylce.Robj')
  
  #with percent.ribo
  EA016_mtrb <- SCTransform(object = EA016, verbose = TRUE, 
                       vars.to.regress = c("percent.mt", "percent.ribo"), return.only.var.genes = FALSE)
  save(EA016_mtrb, file = 'EA016SCT_pctmtribo_nocellcylce.Robj')
  
  #with cell cycle
  EA016_mtcc <- SCTransform(object = EA016, verbose = TRUE, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), return.only.var.genes = FALSE)
  save(EA016_mtcc, file = 'EA016SCT_pctmtcc_cellcylce.Robj')
  
### 8. Optional: First check of clustering --------------------------------------------------------------------------------------------
 ## lot of ribosomal genes in first PC
   EA016_mt <- RunPCA(EA016_mt, verbose = TRUE)
  EA016_mt <- RunUMAP(EA016_mt, dims = 1:30)
  
  
  ## set back defaultassay to integrated if reclustering
  EA016_mt <- FindNeighbors(object = EA016_mt, reduction = "pca", dims = 1:30)
  EA016_mt <- FindClusters(object = EA016_mt, resolution = 0.9)
  
  DimPlot(object = EA016_mt, reduction = "umap", pt.size=1.2)

  FeaturePlot(EA016_mt, 'flowCD31', pt.size=2)
  FeaturePlot(EA016_mt, 'flowCD127', pt.size=2)
  VlnPlot(EA016_mt, 'flowCD31')
  
  ## better: less MT and RP genes
  EA016_mtrb <- RunPCA(EA016_mtrb, verbose = TRUE)
  EA016_mtrb <- RunUMAP(EA016_mtrb, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA016_mtrb <- FindNeighbors(object = EA016_mtrb, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA016_mtrb <- FindClusters(object = EA016_mtrb, resolution = .95)
  
  DimPlot(object = EA016_mtrb, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA016_mtrb, 'flowCD31', pt.size=2)
  FeaturePlot(EA016_mtrb, 'flowCD127', pt.size=2)
  VlnPlot(EA016_mtrb, 'flowCD31')
  
  ## similar to mt
  EA016_mtcc <- RunPCA(EA016_mtcc, veccose = TRUE)
  EA016_mtcc <- RunUMAP(EA016_mtcc, dims = 1:30)
  
  ## set back defaultassay to integrated if reclustering
  EA016_mtcc <- FindNeighbors(object = EA016_mtcc, reduction = "pca", dims = 1:30)
  ## bit higher resolution needed to get 2 clusters -> jumps from 1 to 4..
  EA016_mtcc <- FindClusters(object = EA016_mtcc, resolution = .95)
  
  DimPlot(object = EA016_mtcc, reduction = "umap", pt.size=1.2)
  
  FeaturePlot(EA016_mtcc, 'flowCD31', pt.size=2)
  FeaturePlot(EA016_mtcc, 'flowCD127', pt.size=2)
  VlnPlot(EA016_mtcc, 'flowCD31')