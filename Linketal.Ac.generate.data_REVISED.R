#set library paths and load libraries:
.libPaths(c("/lisc/data/scratch/molevo/agcole/R/libs/seurat4/","/lisc/data/scratch/molevo/agcole/R/libs/course24/","/lisc/opt/sw/software/R/4.5.0/lib64/R/library"))

setwd("/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/Ac_manuscript_revision_ACOE")

library(Seurat,quietly=T)
packageVersion('Seurat') #check that it switched!
library(RColorBrewer,quietly=T)
library(patchwork,quietly=T)
library(ggplot2,quietly=T)
library(pals,quietly=T) #had to add this

library(readxl,quietly=T)
library(tidyr,quietly=T)

LibCP =brewer.paired(16)
#pal.bands(LibCP)
LibCP.stages=LibCP[c(2,7,10,12)]

#* load rawdata
  raw.data1 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/AcPolyp1/outs/filtered_feature_bc_matrix")
  raw.data2 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/A9_AcPolyp2/outs/filtered_feature_bc_matrix")
  raw.data3 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/A10_AcClover/outs/filtered_feature_bc_matrix")
  raw.data4 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/A7_AcStrobila/outs/filtered_feature_bc_matrix")
  raw.data5 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/AcEphyra/outs/filtered_feature_bc_matrix")
  raw.data6 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/A8_AcMedusa/outs/filtered_feature_bc_matrix")
  raw.data7 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/D1_AcM_Gast/outs/filtered_feature_bc_matrix")
  raw.data8 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/B2_AcM_Manub/outs/filtered_feature_bc_matrix")
  raw.data9 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/A1_AcM_Marg/outs/filtered_feature_bc_matrix")
  raw.data10 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/C2_AcM_Umb/outs/filtered_feature_bc_matrix")
  raw.data11 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/F4_AcStrobila2//outs/filtered_feature_bc_matrix")
  raw.data12 <- Read10X(data.dir = "/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/AcEphyra2/outs/filtered_feature_bc_matrix")
  raw.data13 <- Read10X(data.dir = as.character('/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/F3_AcPolypXte/outs/filtered_feature_bc_matrix/'))
  raw.data14 <- Read10X(data.dir = as.character('/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/E10_AcPolypTe/outs/filtered_feature_bc_matrix/'))
  raw.data15 <- Read10X(data.dir = as.character('/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/F7_AcPolypXteRep/outs/filtered_feature_bc_matrix/'))
  raw.data16 <- Read10X(data.dir = as.character('/lisc/data/scratch/molevo/agcole/data/cellranger/results/aurelia/ac.kostyaPlus/F6_AcPolypTeRep/outs/filtered_feature_bc_matrix/'))
  
  #* update the gene names to the annotations
  load(file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE.genes.RData') 
  if(identical(rownames(raw.data1),as.character(genes.ac$cellranger)))
  {
  rownames(raw.data1) <- genes.ac$gene_short_name
  rownames(raw.data2) <- genes.ac$gene_short_name
  rownames(raw.data3) <- genes.ac$gene_short_name
  rownames(raw.data4) <- genes.ac$gene_short_name
  rownames(raw.data5) <- genes.ac$gene_short_name
  rownames(raw.data6) <- genes.ac$gene_short_name
  rownames(raw.data7) <- genes.ac$gene_short_name
  rownames(raw.data8) <- genes.ac$gene_short_name
  rownames(raw.data9) <- genes.ac$gene_short_name
  rownames(raw.data10) <- genes.ac$gene_short_name
  rownames(raw.data11) <- genes.ac$gene_short_name
  rownames(raw.data12) <- genes.ac$gene_short_name
  rownames(raw.data13) <- genes.ac$gene_short_name
  rownames(raw.data14) <- genes.ac$gene_short_name
  rownames(raw.data15) <- genes.ac$gene_short_name
  rownames(raw.data16) <- genes.ac$gene_short_name
  }
  
  #generate Seurat object    
  polyp.1  <- CreateSeuratObject(counts = raw.data1, project = "polyp.1")
  polyp.2  <- CreateSeuratObject(counts = raw.data2, project = "polyp.2")
  polyp.clover  <- CreateSeuratObject(counts = raw.data3, project = "polyp.clover")
  polyp.4b  <- CreateSeuratObject(counts = raw.data13, project = "polyp.clover.body1")
  polyp.5b  <- CreateSeuratObject(counts = raw.data15, project = "polyp.clover.body2")#this is ten1
  polyp.4t  <- CreateSeuratObject(counts = raw.data14, project = "polyp.clover.tent1")#this is body2
  polyp.5t  <- CreateSeuratObject(counts = raw.data16, project = "polyp.clover.tent2")
  
  strobila  <- CreateSeuratObject(counts = raw.data4, project = "strobila.early")
  strobila.2  <- CreateSeuratObject(counts = raw.data11, project = "strobila.late")
  ephyra  <- CreateSeuratObject(counts = raw.data5, project = "ephyra.late")
  ephyra.2  <- CreateSeuratObject(counts = raw.data12, project = "ephyra.early")
  medusa  <- CreateSeuratObject(counts = raw.data6, project = "medusa.early")
  gastric  <- CreateSeuratObject(counts = raw.data7, project = "medusa.gastrodermis")
  manubrium  <- CreateSeuratObject(counts = raw.data8, project = "medusa.manubrium")
  margin  <- CreateSeuratObject(counts = raw.data9, project = "medusa.margin")
  umbrella  <- CreateSeuratObject(counts = raw.data10, project = "medusa.umbrella")
  
  #add library info to names for later identification
  
  polyp.1  <- RenameCells(polyp.1, add.cell.id = "polyp1")
  polyp.2  <- RenameCells(polyp.2, add.cell.id = "polyp2")
  polyp.clover  <- RenameCells(polyp.clover, add.cell.id ="polyp3")
  polyp.4b  <- RenameCells(polyp.4b, add.cell.id ="polyp4b")
  polyp.5b  <- RenameCells(polyp.5b, add.cell.id ="polyp5b")
  polyp.4t  <- RenameCells(polyp.4t, add.cell.id ="polyp4t")
  polyp.5t  <- RenameCells(polyp.5t, add.cell.id ="polyp5t")
  strobila <- RenameCells(strobila, add.cell.id = "strobila1")
  strobila.2 <- RenameCells(strobila.2, add.cell.id = "strobila2")
  ephyra  <- RenameCells(ephyra, add.cell.id = "ephyra1")
  ephyra.2  <- RenameCells(ephyra.2, add.cell.id = "ephyra2")  
  medusa  <- RenameCells(medusa, add.cell.id = "medusa")
  gastric  <- RenameCells(gastric, add.cell.id = "gast")
  manubrium  <- RenameCells(manubrium, add.cell.id ="manub")
  margin  <- RenameCells(margin, add.cell.id ="margin")
  umbrella  <- RenameCells(umbrella, add.cell.id ="umbrel")
  
  #clean up the workspace  
  rm (raw.data1, raw.data2, raw.data3,raw.data4,raw.data5,raw.data6,raw.data7,raw.data8,raw.data9,raw.data10,raw.data11,raw.data12,raw.data13,raw.data14,raw.data15,raw.data16)
  
  Ac.Alldata=merge(x=polyp.1,y= c(polyp.2,polyp.4b,polyp.4t,polyp.5b,polyp.5t,polyp.clover,strobila,strobila.2,ephyra,ephyra.2,medusa,gastric,
                                  manubrium,margin,umbrella), 
                   merge.data = F)
  length(Ac.Alldata@assays$RNA@counts@Dimnames[[1]])
  rm (polyp.1,polyp.2,polyp.clover,polyp.4b,polyp.4t,polyp.5b,polyp.5t,strobila,strobila.2,ephyra,ephyra.2,medusa,gastric,manubrium,margin,umbrella)
  #minimum filter:
  Ac.Alldata <- subset(x = Ac.Alldata, subset = nFeature_RNA > 200 & nCount_RNA > 500)

  single.AllData=SplitObject(Ac.Alldata,split.by = 'orig.ident')
  for (i in 1:length(single.AllData)) 
  {
    print(VlnPlot(single.AllData[[i]],c('nFeature_RNA','nCount_RNA'),group.by = 'orig.ident')+labs(subtitle=paste(names(single.AllData)[i],i))& coord_trans(y = "log10"))#&
    # print(ggplot(single.AllData[[i]]@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_density_2d_filled()+scale_x_log10()+scale_y_log10()+labs(subtitle=names(single.AllData)[i]))
  }

  cut=c(6000,5000,
        5000,5000,4000,3000,
        10000,15000,
        15000,30000,30000,10000,
        30000,10000,30000,30000) #polyp2 seems low??
  cut.genes=c(1200,1200,
              1200,1000,1000,1500,
              2000,3000,
              4000,4000,4000,2000,
              3000,2500,3000,3000)
  cut.low=c(300,450,
            200,200,200,200,
            300,300,
            300,400,400,300,
            300,400,300,425)
  for (i in 1:length(single.AllData))
    #   
    single.AllData[[i]]<-subset(x = single.AllData[[i]],
                                subset = nFeature_RNA > cut.low[i] & nFeature_RNA < cut.genes[i] & nCount_RNA < cut[i])# &  & 
  Ac.Alldata.unprocessed=merge(single.AllData[[1]],c(single.AllData[[2]],single.AllData[[3]],
                                                     single.AllData[[4]],single.AllData[[5]],single.AllData[[6]],single.AllData[[7]],single.AllData[[8]],
                                                     single.AllData[[9]],single.AllData[[10]],single.AllData[[11]],single.AllData[[12]],single.AllData[[13]],single.AllData[[14]],single.AllData[[15]],single.AllData[[16]]))
  
  save(Ac.Alldata.unprocessed,file='AcAlldata.raw.16.Robj')
  save(single.AllData,file='single.Alldata.16.ROBJ')

# pull out the total genes detected and median gene for each
  x=matrix(0,16,3)
  for (i in 1:length(single.AllData))
    {x[i,2]=length(which(rowSums(single.AllData[[i]]@assays$RNA@counts)>=1))
    x[i,3]=median(single.AllData[[i]]$nCount_RNA)
    x[i,1]=names(single.AllData)[i]}
