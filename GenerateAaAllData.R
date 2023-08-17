setwd("C:/Users/Alison/bioinformatics/R_analyses/Aurelia") ## Where you want to work
#NOTE: this stretches the limits of my local capacity; perhaps better run on the server.
#load libraries: 
library(easypackages)
libraries("Seurat", "RColorBrewer",'patchwork','viridis','ggplot2',
          'pals','Rmagic', "readxl",'SeuratWrappers')
setup=F
# setup genes ----
if (setup)
{
  #load and update gene names...
  Aa_genes <- read_excel("DataS3.xls",sheet='features_CellRangerOutput')
  Aa_genes <- as.data.frame(Aa_genes)
  Aa_annotations <- read_excel("DataS3.xls",sheet='GeneAnnotations')
  #merge will work if column names are the same.
  Aa_genes=merge(Aa_genes, Aa_annotations, by="Cellranger.Aaur2", all.x=T, sort = F)
  Aa_genes = Aa_genes[c(1,3,4,6,10,11,17)]
  
  mito.genes <- grep(pattern = "mitochondrial", Aa_genes$best_OG_desc)
  mito_genes = Aa_genes[mito.genes,2]
  rm(mito.genes)
  
  ribo.genes <- grep(pattern = "ribosome", Aa_genes$best_OG_desc)
  ribo_genes = Aa_genes[ribo.genes,2]
  rm(ribo.genes)
  
  AaTF_list <- which(!(Aa_genes$TF_Fam=='-'))
  AaTF_list = Aa_genes$name[AaTF_list]
  
  save(OMA,Aa_annotations,Aa_genes,AaTF_list, file="Aaur2.newnames.RData")
  
  #optional: save workspace here and then just load '.RData' each time you start and skip the begining
}else{load (file="Aaur2.newnames.RData")}

load.rawdata = F
run.integration = T

if (load.rawdata)
{  #direct to the matrix files of interest here.
  
  raw.data1 <- Read10X(data.dir = "/TR2 mapping/polyp")
  raw.data2 <- Read10X(data.dir = "/TR2 mapping/polyp2")
  raw.data3 <- Read10X(data.dir = "/TR2 mapping/clover.polyp")
  raw.data4 <- Read10X(data.dir = "/TR2 mapping/strobila")
  raw.data5 <- Read10X(data.dir = "/TR2 mapping/ephyra")
  raw.data6 <- Read10X(data.dir = "/TR2 mapping/medusa")
  raw.data7 <- Read10X(data.dir = "/TR2 mapping/mT_gastrodermis")
  raw.data8 <- Read10X(data.dir = "/TR2 mapping/mT_manubrium")
  raw.data9 <- Read10X(data.dir = "/TR2 mapping/mT_margin")
  raw.data10 <- Read10X(data.dir = "/TR2 mapping/mT_umbrella")
  
  # # set the gene names to the Aa_annotations
  rownames(raw.data1) <- Aa_genes$gene_short_name
  rownames(raw.data2) <- Aa_genes$gene_short_name
  rownames(raw.data3) <- Aa_genes$gene_short_name
  rownames(raw.data4) <- Aa_genes$gene_short_name
  rownames(raw.data5) <- Aa_genes$gene_short_name
  rownames(raw.data6) <- Aa_genes$gene_short_name
  rownames(raw.data7) <- Aa_genes$gene_short_name
  rownames(raw.data8) <- Aa_genes$gene_short_name
  rownames(raw.data9) <- Aa_genes$gene_short_name
  rownames(raw.data10) <- Aa_genes$gene_short_name
  
  #generate Seurat object    
  polyp.1  <- CreateSeuratObject(counts = raw.data1, project = "polyp")
  polyp.2  <- CreateSeuratObject(counts = raw.data2, project = "polyp2")
  polyp.clover  <- CreateSeuratObject(counts = raw.data3, project = "polyp.clover")
  strobila  <- CreateSeuratObject(counts = raw.data4, project = "strobila")
  ephyra  <- CreateSeuratObject(counts = raw.data5, project = "ephyra")
  medusa  <- CreateSeuratObject(counts = raw.data6, project = "medusa")
  gastric  <- CreateSeuratObject(counts = raw.data7, project = "mT.gastrodermis")
  manubrium  <- CreateSeuratObject(counts = raw.data8, project = "mT.manubrium")
  margin  <- CreateSeuratObject(counts = raw.data9, project = "mT.margin")
  umbrella  <- CreateSeuratObject(counts = raw.data10, project = "mT.umbrella")
  
  #add library info to names for later identification
  
  polyp.1  <- RenameCells(polyp.1, add.cell.id = "AaPolyp1")
  polyp.2  <- RenameCells(polyp.2, add.cell.id = "AaPolyp2")
  polyp.clover  <- RenameCells(polyp.clover, add.cell.id ="polyp.clover")
  strobila <- RenameCells(strobila, add.cell.id = "strobila")
  ephyra  <- RenameCells(ephyra, add.cell.id = "ephyra")
  medusa  <- RenameCells(medusa, add.cell.id = "medusa")
  gastric  <- RenameCells(gastric, add.cell.id = "mT.gastrodermis")
  manubrium  <- RenameCells(manubrium, add.cell.id ="mT.manubrium")
  margin  <- RenameCells(margin, add.cell.id ="mT.margin")
  umbrella  <- RenameCells(umbrella, add.cell.id ="mT.umbrella")
  
  #clean up the workspace  
  rm (raw.data1, raw.data2, raw.data3,raw.data4,raw.data5,raw.data6,raw.data7,raw.data8,raw.data9,raw.data10)
  
  Aa.Alldata.unprocessed=merge(x=polyp.1,y= c(polyp.2,polyp.clover,strobila,ephyra,medusa,gastric,
                                  manubrium,margin,umbrella), 
                   merge.data = F)
  length(Aa.Alldata.unprocessed@assays$RNA@counts@Dimnames[[1]])
  
  VlnPlot(object = Aa.Alldata.unprocessed, features = c("nFeature_RNA", "nCount_RNA"),# "percent.mt"),
          group.by = 'orig.ident',cols = LibCP)#
  Aa.Alldata.unprocessed <- subset(x = Aa.Alldata.unprocessed, subset = nFeature_RNA > 250 & nCount_RNA > 750)#I think min 500 gene filter removes a LOT of the polyp data
  save(Aa.Alldata.unprocessed,file='AaAlldata.raw.updatenames.Robj')
  # filter each separately:
 
    single.AllData=SplitObject(Aa.Alldata.unprocessed,split.by = 'orig.ident')
    for (i in 1:length(single.AllData)) 
    {
      print(VlnPlot(single.AllData[[i]],c('nFeature_RNA','nCount_RNA'),group.by = 'orig.ident')+
              labs(subtitle=c(i,names(single.AllData)[i])))
      print(single.AllData[[i]]@meta.data %>%
              ggplot(aes(x = nFeature_RNA)) +
              geom_density() +
              scale_x_log10() +
              theme_classic() +
              ylab("Cell density") +
              geom_vline(xintercept = c(300, 400, 500))+
        labs(subtitle=names(single.AllData)[[i]]))
    }
    cut=c(10000,4000,10000,20000,30000,10000,50000,15000,40000,40000) #polyp2 seems low??
    cut.genes=c(1500,1200,2000,4000,4000,1500,4000,3000,4000,4000)
    cut.low=c(350,450,400,400,400,400,500,450,500,425)
    for (i in 1:length(single.AllData))
      #   
      single.AllData[[i]]<-subset(x = single.AllData[[i]],
                                  subset = nFeature_RNA > cut.low[i] & nFeature_RNA < cut.genes[i] & nCount_RNA < cut[i])# &  & 
    Aa.Alldata.unprocessed=merge(single.AllData[[1]],c(single.AllData[[2]],single.AllData[[3]],
                                                       single.AllData[[4]],single.AllData[[5]],single.AllData[[6]],single.AllData[[7]],single.AllData[[8]],
                                                       single.AllData[[9]],single.AllData[[10]]))
    
    save(Aa.Alldata.unprocessed,file='AaAlldata.raw.updatenames.Robj')
    save(single.AllData,file='single.Alldata.updatenames.ROBJ')
    }

if (run.integration)
{
  load(file='single.Alldata.updatenames.ROBJ')
    # normalize and identify variable features for each dataset independently
    single.AllData <- lapply(X = single.AllData, FUN = function(x) {
      x <- NormalizeData(x, scale.factor = 5000)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    
    #make a list of ALL variable features"
    vargenelist= NULL
    for (i in 1:length(single.AllData)) {
      x <- single.AllData[[i]]@assays$RNA@var.features
      vargenelist=c(vargenelist,x)}
    vargenelist=unique(vargenelist)
    length(vargenelist)
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = single.AllData,nfeatures=2000)
    # remove cell.cycle genes from the anchor list
    features.red= setdiff(features,cell.cycle.goi)
    
    #scale data and calculate PCA individually
    for (i in 1:length(single.AllData)) {
      single.AllData[[i]]@assays$RNA@var.features = features.red}
    
    single.AllData <- lapply(X = single.AllData, FUN = function(x) {
      x <- ScaleData(x,features=features.red)
      x <- RunPCA(x)
    })
    save (single.AllData, file= 'single.Alldata.updatenames.ROBJ')
    #* Perform integration
    anchors.rpca <- FindIntegrationAnchors(object.list = single.AllData,
                                           scale=F,
                                           normalization.method = 'LogNormalize',
                                           reduction='rpca',
                                           anchor.features = features.red)
    save(anchors.rpca,file='anchors.single.updatenames.RObj')
}
if (1)
{
    load(file='anchors.single.updatenames.RObj')
    #if run in parts, this would be needed to reload anchors and integrate
      vargenelist= NULL
      for (i in 1:10) {
        x <- anchors.rpca@object.list[[i]]@assays$RNA@var.features
        vargenelist=c(vargenelist,x)}
      vargenelist=unique(vargenelist)
      length(vargenelist)
  # this command creates an 'integrated' data assay
    rm(single.AllData)
  goi.integrate = unique(c(vargenelist,#all variable genes from the separate datasets
                           AaTF_list,cell.cycle.goi))#,#all expressed TFs
  
  Aa.Alldata <- IntegrateData(anchorset = anchors.rpca, 
                              features.to.integrate = goi.integrate)
  DefaultAssay(Aa.Alldata) <- "integrated"
  Aa.Alldata@assays$integrated@var.features=(Aa.Alldata@assays$integrated@var.features)
  
  save(Aa.Alldata,file='AaAlldata.integrated.updatenames.Robj')
}
    