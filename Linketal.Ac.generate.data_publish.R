
library(Matrix)
library(Seurat,quietly=T)
packageVersion('Seurat') #check that it switched!
library(RColorBrewer,quietly=T)
library(patchwork,quietly=T)
library(ggplot2,quietly=T)
library(pals,quietly=T) #had to add this

library(readxl,quietly=T)
library(tidyr,quietly=T)

load (file="Aaur2.newnames.RData")

LibCP =brewer.paired(12)
LibCP.stages=LibCP[c(2,4,5,8)]
{
  #* load rawdata
  raw.data1 <- Read10X(data.dir = "acpolyp1/outs/filtered_feature_bc_matrix")
  raw.data2 <- Read10X(data.dir = "aapolyp2_20000cells/outs/filtered_feature_bc_matrix")
  raw.data3 <- Read10X(data.dir = "cloverpolyp/outs/filtered_feature_bc_matrix")
  raw.data4 <- Read10X(data.dir = "strobila_2000/outs/filtered_feature_bc_matrix")
  raw.data5 <- Read10X(data.dir = "ephyra_1000/outs/filtered_feature_bc_matrix")
  raw.data6 <- Read10X(data.dir = "medusa_1000/outs/filtered_feature_bc_matrix")
  raw.data7 <- Read10X(data.dir = "medusaGast/outs/filtered_feature_bc_matrix")
  raw.data8 <- Read10X(data.dir = "medusaManubrium/outs/filtered_feature_bc_matrix")
  raw.data9 <- Read10X(data.dir = "medusaMargin/outs/filtered_feature_bc_matrix")
  raw.data10 <- Read10X(data.dir = "medusaUmbrella/outs/filtered_feature_bc_matrix")
  raw.data11 <- Read10X(data.dir = "strobila2/outs/filtered_feature_bc_matrix")
  raw.data12 <- Read10X(data.dir = "ephyra_ACME/outs/filtered_feature_bc_matrix")
  
  #* update the gene names to the Aa_annotations
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
  rownames(raw.data11) <- Aa_genes$gene_short_name
  rownames(raw.data12) <- Aa_genes$gene_short_name
  
  #generate Seurat object    
  polyp.1  <- CreateSeuratObject(counts = raw.data1, project = "Ac.polyp.1")
  polyp.2  <- CreateSeuratObject(counts = raw.data2, project = "Ac.polyp.2")
  polyp.clover  <- CreateSeuratObject(counts = raw.data3, project = "Ac.polyp.clover")
  strobila  <- CreateSeuratObject(counts = raw.data4, project = "Ac.strobila.early")
  strobila.2  <- CreateSeuratObject(counts = raw.data11, project = "Ac.strobila.late")
  ephyra  <- CreateSeuratObject(counts = raw.data5, project = "Ac.ephyra.late")
  ephyra.2  <- CreateSeuratObject(counts = raw.data12, project = "Ac.ephyra.early")
  medusa  <- CreateSeuratObject(counts = raw.data6, project = "Ac.medusa.early")
  gastric  <- CreateSeuratObject(counts = raw.data7, project = "Ac.medusa.gastrodermis")
  manubrium  <- CreateSeuratObject(counts = raw.data8, project = "Ac.medusa.manubrium")
  margin  <- CreateSeuratObject(counts = raw.data9, project = "Ac.medusa.margin")
  umbrella  <- CreateSeuratObject(counts = raw.data10, project = "Ac.medusa.umbrella")
  
  #add library info to names for later identification
  
  polyp.1  <- RenameCells(polyp.1, add.cell.id = "polyp1")
  polyp.2  <- RenameCells(polyp.2, add.cell.id = "polyp2")
  polyp.clover  <- RenameCells(polyp.clover, add.cell.id ="polyp3")
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
  rm (raw.data1, raw.data2, raw.data3,raw.data4,raw.data5,raw.data6,raw.data7,raw.data8,raw.data9,raw.data10,raw.data11)
  
  Aa.Alldata=merge(x=polyp.1,y= c(polyp.2,polyp.clover,strobila,strobila.2,ephyra,ephyra.2,medusa,gastric,
                                  manubrium,margin,umbrella), 
                   merge.data = F)
  length(Aa.Alldata@assays$RNA@counts@Dimnames[[1]])
  rm (polyp.1,polyp.2,polyp.clover,strobila,strobila.2,ephyra,medusa,gastric,manubrium,margin,umbrella)
  #minimum filter:
  Aa.Alldata <- subset(x = Aa.Alldata, subset = nFeature_RNA > 300 & nCount_RNA > 500)
  
  single.AllData=SplitObject(Aa.Alldata,split.by = 'orig.ident')
  for (i in 1:length(single.AllData)) 
  {
    print(VlnPlot(single.AllData[[i]],c('nFeature_RNA','nCount_RNA'),group.by = 'orig.ident')+labs(subtitle=c(i,names(single.AllData)[i]))& coord_trans(y = "log10"))#&
   }
  cut=c(10000,4000,10000,20000,25000,30000,40000,10000,40000,15000,40000,40000) #polyp2 seems low??
  cut.genes=c(1500,1200,2000,4000,4000,4000,4000,2000,4000,3000,4000,4000)
  cut.low=c(300,450,300,300,300,400,400,300,300,400,300,425)
  for (i in 1:length(single.AllData))
    #   
    single.AllData[[i]]<-subset(x = single.AllData[[i]],
                                subset = nFeature_RNA > cut.low[i] & nFeature_RNA < cut.genes[i] & nCount_RNA < cut[i])# &  & 
  Aa.Alldata.unprocessed=merge(single.AllData[[1]],c(single.AllData[[2]],single.AllData[[3]],
                                                     single.AllData[[4]],single.AllData[[5]],single.AllData[[6]],single.AllData[[7]],single.AllData[[8]],
                                                     single.AllData[[9]],single.AllData[[10]],single.AllData[[11]],single.AllData[[12]]))
  
  save(Aa.Alldata.unprocessed,file='Ac_manuscript_final/AaAlldata.raw.new.Robj')
  save(single.AllData,file='Ac_manuscript_final/single.Alldata.new.ROBJ')
}