#set library paths and load libraries:
.libPaths(c("/lisc/data/scratch/molevo/agcole/R/libs/seurat4/","/lisc/data/scratch/molevo/agcole/R/libs/course24/"))

setwd("/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/Ac_manuscript_revision_ACOE")

library(Seurat,quietly=T)
packageVersion('Seurat') #check that it switched!
library(RColorBrewer,quietly=T)
library(patchwork,quietly=T)
library(ggplot2,quietly=T)
library(pals,quietly=T) #had to add this
library(readxl,quietly=T)
library(tidyr,quietly=T)
library(dplyr,quietly=T)

load (file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE.genes.RData')

LibCP = c(stepped3(3),stepped(20)[c(13:16)],stepped(20)[c(11,9)],stepped2(20)[c(10,11)],stepped3(20)[6],stepped(20)[c(1,2)],stepped2(20)[c(17,18)])
LibCP.stages=LibCP[c(1,9,10,13)]

gene.cp = c('lightgrey', rev(brewer.pal(11 , "Spectral")))
clust.cp.separate = unique (c(cols25(25), glasbey(32), alphabet(26)))
clust.cp.graded = unique(c(stepped3(20), stepped2(20), stepped(20)))
get_marker_list <- function(all.markers, N=5)
{
  #generate a collated list of unique DE genes
  list = NULL
  for (i in 1:length(unique(all.markers$cluster)))
  {
    x = all.markers[as.numeric(all.markers$cluster) == i, ][1:min(N, length(which(as.numeric(all.markers$cluster) == i))), 7]
    
    list = c(list, x)
    list= list[!is.na(list)]
  }
  
  return(list)
}
update_marker_list <- function(all.markers)
{ 
  all.markers[,8:11] = 'NA'
  names(all.markers)[8:11]=names(genes.ac)[c(3:6)]
  # add GO terms associated with this list:
  for (i in 1:length(levels(data1@active.ident))) # 
  {
    x=all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),7]
    x2=genes.ac[match(x,genes.ac$gene_short_name),c(3:6)]
    all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),8:11]<-x2
  }  
  
  return(all.markers)
}
do.over=T
# subset analysis ----
if (do.over)
{ 
load (file='Ac.Alldata.Robj')
data1<-SetIdent(Ac.Alldata,value='IDs')
DimPlot(data1,cols=clust.cp.separate)+coord_flip()
data1.subsets <- SplitObject(Ac.Alldata,split.by = 'IDs')

order = match(levels(Ac.Alldata),names(data1.subsets))
data1.subsets=data1.subsets[order]
{
  for (i in 1:length(names (data1.subsets)))
  {
    data1=data1.subsets[[i]]
    data1@active.assay = 'RNA'
    data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
    # #use the full dataset scaling:
    Ac.Alldata@active.assay='RNA'
    coi = colnames(data1)
    t=ScaleData(Ac.Alldata,features = data1@assays$RNA@var.features,split.by = 'orig.ident')

    data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
    rm(t)
    # data1 <- ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50,verbose = F) #this is now failing...
    data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')

    pct <- data1[["harmony"]]@stdev / sum(data1[["harmony"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    d = 1:co2
    
    data1 <- RunUMAP(data1, reduction ="harmony", n.neighbors = 10L,spread =1,seed.use = 0, 
                     dims = d,reduction.name ='umap',
                     reduction.key ='umap',min.dist = 0.3,
                     local.connectivity = 10,verbose = F)#
    data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                           nn.method = 'annoy',  annoy.metric = 'cosine',verbose = F)
    data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0,verbose = F)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
    data1$IDs.fine = data1@active.ident
    data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0,verbose = F)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
    data1$IDs.coarse = data1@active.ident
    data1@reductions=data1@reductions[c(1,2,4,3)]
    data1.subsets[[i]] = data1
  }#this sets coarse clustering as default
}
for (i in 1:length(names (data1.subsets)))
{
  data1=data1.subsets[[i]]
print(DimPlot(data1,cols=clust.cp.separate,group.by = 'lifehistory',label=F,reduction = 'umap')+labs(title=names(data1.subsets[i]))&DimPlot(data1,cols=clust.cp.separate,#group.by = 'IDs',
label=F,reduction = 'umap')+labs(title=names(data1.subsets[i])))
}
##refine and label ----
####outer ----
{

  data1=data1.subsets$epidermis
  data1<-SetIdent(data1,value='IDs.coarse')
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'OuterClusters')
  
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
   data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='RNA'
  print(DimPlot(data1,group.by = c('IDs'),reduction='umap',label=F,cols = (clust.cp.separate))+NoAxes())
  data1.subsets$epidermis=data1
}
####st.muscle ----
{
  data1=data1.subsets$muscle.st
  #do-over clustering
  data1<-SetIdent(data1,value = 'IDs.coarse')
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'Muscle.ST') 
  
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1.subsets$muscle.st=data1
}
##### all .muscle ----
{
  Ac.Alldata<-SetIdent(Ac.Alldata,value ='IDs')
  coi1=WhichCells(Ac.Alldata,idents = levels(Ac.Alldata)[1])
  coi2 = WhichCells(data1.subsets$epidermis,idents=levels(data1.subsets$epidermis)[1])
  allmuscle=subset(Ac.Alldata,cells=c(coi1,coi2))
  data1=SetIdent(allmuscle,value = 'IDs')
  data1@active.assay='RNA'
  data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
  # #use the full dataset scaling:
  Ac.Alldata@active.assay='RNA'
  coi = colnames(data1)
  t=ScaleData(Ac.Alldata,features = data1@assays$RNA@var.features)
  data1@active.assay='RNA'
  data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
  rm(t)
  
  data1<-ScaleData(data1)
  data1 <- RunPCA(data1, pcs.compute = 50,verbose = F)
  #add harmony integration here!
  data1<- harmony::RunHarmony(data1,'orig.ident')
  d = 1:20
  data1 <- RunUMAP(data1, reduction ="harmony", n.neighbors = 10L,spread =1,seed.use = 42, 
                   dims = d,reduction.name ='umap',
                   reduction.key ='umap',min.dist = 0.4,
                   local.connectivity = 1,verbose = F)#
  DimPlot(data1,cols=clust.cp.separate,reduction = 'umap')&NoAxes()
  DimPlot(data1,group.by = 'lifehistory',reduction ='umap')
  data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                         nn.method = 'annoy',  annoy.metric = 'cosine',verbose = F)
  data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0,verbose = F)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
  data1$IDs.fine = data1@active.ident
  data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0,verbose = F)
  # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
  data1$IDs.coarse = data1@active.ident
  
  DimPlot(data1,cols=clust.cp.separate,reduction='umap')
  
  # manually re-annotate the clusters:
  cl.names= c('muscle.st','muscle.st','muscle.st','smooth.muscle','smooth.muscle','smooth.muscle','smooth.muscle')
  levels(data1@active.ident) = cl.names
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  DimPlot(data1,cols=clust.cp.separate,reduction='umap')
  
  allmuscle=data1
  save(allmuscle,file='New.allmuscle.RObj')
  
  # update into AllData:
  load(file='New.allmuscle.RObj')
  coi.st=WhichCells(allmuscle,ident='muscle.st')
  coi.sm = WhichCells(allmuscle,ident = 'smooth.muscle')
  # Ac.Alldata$ID.separate = as.character(Ac.Alldata$ID.separate)
  # Ac.Alldata$ID.separate[coi.sm] = as.character(allmuscle@active.ident[coi.sm])
  Ac.Alldata$IDs = as.character(Ac.Alldata$IDs)
  Ac.Alldata$IDs[coi.st] = as.character(allmuscle@active.ident[coi.st])
  Ac.Alldata$IDs[coi.sm] = rep("epidermis", length(coi.sm))
  Ac.Alldata$IDs <- as.factor(Ac.Alldata$IDs)
  Ac.Alldata$IDs<-factor(Ac.Alldata$IDs,levels(Ac.Alldata$IDs)[c(6,2,3,5,7,1,4,8)])
  }
#* re-run the st. muscle AND outer with updated subsets:
####outer update----
{
  Ac.Alldata<-SetIdent(Ac.Alldata,value='IDs')
  data1=subset(Ac.Alldata,idents='epidermis')
  {data1@active.assay = 'RNA'
    data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
    # #use the full dataset scaling:
    Ac.Alldata@active.assay='RNA'
    coi = colnames(data1)
    t=ScaleData(Ac.Alldata,features = data1@assays$RNA@var.features,split.by = 'orig.ident')
    
    data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
    rm(t)
    # data1 <- ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50,verbose = F) #this is now failing...
    data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')
    
    pct <- data1[["harmony"]]@stdev / sum(data1[["harmony"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    d = 1:co2
    
    data1 <- RunUMAP(data1, reduction ="harmony", n.neighbors = 10L,spread =1,seed.use = 0, 
                     dims = d,reduction.name ='umap',
                     reduction.key ='umap',min.dist = 0.3,
                     local.connectivity = 10,verbose = F)#
    data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                           nn.method = 'annoy',  annoy.metric = 'cosine',verbose = F)
    data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0,verbose = F)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
    data1$IDs.fine = data1@active.ident
    data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0,verbose = F)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
    data1$IDs.coarse = data1@active.ident
    data1@reductions=data1@reductions[c(1,2,4,3)]
    data1.subsets$epidermis=data1
    }             
  data1<-SetIdent(data1,value='IDs.coarse')
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'OuterClusters')
  
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='RNA'
  print(DimPlot(data1,group.by = c('IDs'),reduction='umap',label=F,cols = (clust.cp.separate))+NoAxes())
  data1.subsets$epidermis=data1
}
####st.muscle update----
{
  Ac.Alldata<-SetIdent(Ac.Alldata,value='IDs')
  data1=subset(Ac.Alldata,idents='muscle.st')
  {data1@active.assay = 'RNA'
    data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
    # #use the full dataset scaling:
    Ac.Alldata@active.assay='RNA'
    coi = colnames(data1)
    t=ScaleData(Ac.Alldata,features = data1@assays$RNA@var.features,split.by = 'orig.ident')
    
    data1@assays$RNA@scale.data = t@assays$RNA@scale.data[,coi]
    rm(t)
    # data1 <- ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50,verbose = F) #this is now failing...
    data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')
    
    pct <- data1[["harmony"]]@stdev / sum(data1[["harmony"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    d = 1:co2
    
    data1 <- RunUMAP(data1, reduction ="harmony", n.neighbors = 10L,spread =1,seed.use = 0, 
                     dims = d,reduction.name ='umap',
                     reduction.key ='umap',min.dist = 0.3,
                     local.connectivity = 10,verbose = F)#
    data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                           nn.method = 'annoy',  annoy.metric = 'cosine',verbose = F)
    data1 <- FindClusters(object = data1,resolution = 1,random.seed = 0,verbose = F)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
    data1$IDs.fine = data1@active.ident
    data1 <- FindClusters(object = data1,resolution = 0.2,random.seed = 0,verbose = F)
    # #look at relationship between clusters
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
    data1$IDs.coarse = data1@active.ident
    data1@reductions=data1@reductions[c(1,2,4,3)]
    data1.subsets$muscle.st=data1
  } 
  #do-over clustering
  data1<-SetIdent(data1,value = 'IDs.coarse')
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'Muscle.ST') 
  
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1.subsets$muscle.st=data1
}
####inner ----
{
  data1=data1.subsets$gastrodermis
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'InnerClusters') 
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  # data1@reductions=data1@reductions[c(1,5,4,2,3)]
  print(DimPlot(data1,group.by = c('IDs'),label=F,cols = (clust.cp.separate),reduction='umap')+NoAxes())
  data1.subsets$gastrodermis=data1
}
####cnido ----
{
  
  data1=data1.subsets$cnidocyte #
  data1<-SetIdent(data1,value = 'IDs.coarse')
  data1 <- RunUMAP(data1, reduction ="harmony", n.neighbors = 10L,spread =1,seed.use = 15, 
                   dims = 1:30,reduction.name ='umap',
                   reduction.key ='umap',min.dist = 0.3,
                   local.connectivity = 1,verbose = F)
  # recalculated just to have a similar topology as the nematostella for comparasons.
  
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'CnidoClusters')
  
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')

    #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident

  data1@active.assay='RNA'
  data1.subsets$cnidocyte=data1
  print(DimPlot(data1,cols=clust.cp.separate,reduction='umap',group.by = 'IDs')+NoAxes()+scale_x_reverse())#+
}
####mucin ----
{
  data1=data1.subsets$gland.muc
  data1<-SetIdent(data1,value = 'IDs.coarse')#fine picks up small clusters but don't know what to do with these
  # print(DimPlot(data1,cols = clust.cp.separate,label = T)&NoAxes())
  # data1@reductions <- data1@reductions[c(1,3,2)]
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'MucinClusters')
  
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  DimPlot(data1,label = T,cols=clust.cp.separate,reduction='umap')+NoAxes()
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='RNA'
  data1.subsets$gland.muc =data1
  
  print(DimPlot(data1,group.by = c('IDs'),label=F,cols = (clust.cp.separate),)+NoAxes())
  
  
}
####dig.gland ----
{
  data1=data1.subsets$gland.dig
  # data1<-SetIdent(data1,value = 'IDs.coarse')        
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'Dig.Gland')
  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='RNA'
  # data1@reductions=data1@reductions[c(1,3,2)]
print(DimPlot(data1,group.by = c('IDs'),label=F,reduction = 'umap', cols = (clust.cp.separate))+NoAxes())
  
  data1.subsets$gland.dig=data1
  
}
####neuron----
{
  data1=data1.subsets$neural
  data1<-SetIdent(data1,value='IDs.fine')#fine or coarse???
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx', sheet = 'Neuro')

  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0,split.by = 'lifehistory',cols=LibCP.stages)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  #make sure this is the expected size:
  if (length (rownames(m$RNA)) == length (goi))
  {  
    for (j in 1:length(levels(data1@active.ident))) #for each cluster set
    {
      clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
    }
    sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
    clusterNames$ID[clName]
    #Note: here we merge the cnidocytes, dig. gland, and neural
  } else
    print('error - check gene lists')
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='RNA'
  # data1@reductions=data1@reductions[c(1,3,2)]
  data1.subsets$neural=data1
  
  print(DimPlot(data1,group.by = c('IDs'),label=T,reduction = 'umap',
  cols = (clust.cp.separate))+NoAxes()+labs(title = 'neural population'))
}
#### OTHER ----
{
  data1=data1.subsets$unchar.1 #
 data1<-SetIdent(data1,value = 'IDs')
 data1.subsets$unchar.1 = data1

}

## add ids to alldata ----
{
  cl.order=NULL
  cl.ind = NULL
  # Ac.Alldata$IDs = Ac.Alldata@active.ident 
  Ac.Alldata$ID.separate = as.character(Ac.Alldata$IDs)
  
  for (i in 1:length(names (data1.subsets)))
    
  {
    coi=NULL
    coi=colnames(data1.subsets[[i]])
    Ac.Alldata$ID.separate[coi] = as.character(data1.subsets[[i]]@active.ident[coi])
    cl.order = c(cl.order,levels(data1.subsets[[i]]))
  }
  cl.ind = match(cl.order,levels(as.factor(Ac.Alldata$ID.separate))) #doesn't account for 'skipped' clusters...
  # cl.ind=c(cl.ind,setdiff(c(1:62),cl.ind))
  Ac.Alldata$ID.separate = as.factor(Ac.Alldata$ID.separate)
  Ac.Alldata$ID.separate = factor(Ac.Alldata$ID.separate,levels(Ac.Alldata$ID.separate)[cl.ind])
  
}
save(Ac.Alldata,file='Ac.Alldata.Robj')
save(data1.subsets,file='Ac.subsets.RObj')
}else{
  load (file='Ac.subsets.RObj')
  load (file='Ac.Alldata.Robj')
}
##calculate medus-only ----
{
  data1=SetIdent(Ac.Alldata,value='ID.separate')
  
  #  use the idents table:
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory))
  colnames(ids.cluster.library) = c('ID','Lifehistory','CellCount')
 
  medusa.cl=as.character(ids.cluster.library$ID[which(ids.cluster.library$CellCount[ids.cluster.library$Lifehistory=='polyp'] <= 3)])
  medusa.cl
}

## all.idents.DEGs ----
{
  data1=SetIdent(Ac.Alldata,value = 'ID.separate')
data1@active.assay='RNA'
all.markers.allidents <- FindAllMarkers(data1,
                                        logfc.threshold = 1,
                                        return.thresh = 0.0001,
                                        min.pct = 0.2,
                                        only.pos = TRUE, 
                                        verbose=T)

all.markers.allidents <- update_marker_list(all.markers.allidents)
#also for transcription factors only:
all.markers.TF.allidents <- FindAllMarkers(data1,
                                           logfc.threshold = 0.2,
                                           features = intersect(AcTF_list,rownames(data1)),
                                           min.pct = 0.05,
                                           only.pos = TRUE,
                                           return.thresh = 0.001,verbose=T)

all.markers.TF.allidents <- update_marker_list(all.markers.TF.allidents)
save (all.markers.TF.allidents,all.markers.allidents,file='all.idents.markers.RData')

#list = get_marker_list(all.markers.allidents,5)
#list2 = get_marker_list(all.markers.TF.allidents,5)

# print(DotPlot(Ac.Alldata,'RNA', features = unique(list), group.by = 'ID.separate',
#               scale.by='radius' , col.min = 0, col.max = 4,
#               cols = c('lightgrey','darkred'))+
#         RotatedAxis() +FontSize(6,8) +labs(title='DEG')+theme(legend.position = 'bottom',legend.title = element_text(size=8)))+coord_flip()+theme(panel.grid = element_line(color='black',linewidth = 0.2))
#print(
  # DotPlot(data1,'RNA', features = unique(list2),
  #         scale.by='radius' , col.min = 0, col.max = 4,
  #         cols = c('lightgrey','darkred'))+
  #   RotatedAxis() +FontSize(6,8) +labs(title ='DETFs'))+theme(panel.grid = element_line(color='black',linewidth = 0.2))
}
## Generate gene lists subsets ----
{
  markers = NULL
  markers.TF = NULL# 
  for (i in 1:7)#length(names (data1.subsets)))
  {
    data1=data1.subsets[[i]]
    data1@active.assay='RNA'
    all.markers <- FindAllMarkers(data1,
                                  logfc.threshold = 0.6,
                                  return.thresh = 0.001,
                                  min.pct = 0.2,
                                  only.pos = TRUE,verbose = F)

    all.markers <- update_marker_list(all.markers)
    markers[[i]] = all.markers

    #also for transcription factors only: this did not work...
    all.markers_TF <- FindAllMarkers(data1,
                                     features = intersect(AcTF_list,rownames(data1)),
                                     min.pct = 0.01,
                                     only.pos = TRUE,
                                     return.thresh = 0.001,verbose=F)
    
    all.markers_TF <- update_marker_list(all.markers_TF) 
    markers.TF[[i]] = all.markers_TF
  }
    

  # View(markers[[i]])
  names(markers) = names(data1.subsets)[1:7] #'unchar.1 is not subclustered
  names(markers.TF) = names(data1.subsets)[1:7]
}  
  save(markers,file = 'Ac.AllData.subsets.markers.RObj')
  save(markers.TF,file = 'Ac.AllData.subsets.TF.markers.RObj')
  for (i in 1:7)
  {
    write.csv(markers[[i]],file=paste('DS1.',paste(i+4,".csv",sep=''),sep = ''))
    write.csv(markers.TF[[i]],file=paste('DS1.',paste(i+4,"b.TF.csv",sep=''),sep = ''))
   }

  
## cell type tree neurons / lifehistory----
data1<-data1.subsets$neural#Ac.Alldata#
data1$split.ID = paste (data1$IDs,data1$lifehistory)
# data1$split.ID.separate = paste (data1$ID.separate,data1$lifehistory)
data1$split.ID.separate = data1$split.ID
data1$split.ID.separate=as.factor(data1$split.ID.separate)
order = sort(levels(data1$split.ID.separate),index.return=T)
data1$split.ID.separate=factor(data1$split.ID.separate,levels = levels(data1$split.ID.separate)[order$ix])
data1$split.ID.separate=droplevels(data1$split.ID.separate)
data1<-SetIdent(data1,value='split.ID.separate')
levels(data1)
embeddings <- Embeddings(object = data1, reduction = 'harmony')[,1:50]
data.dims=matrix('0',50L,length(levels(data1)))
for (i in 1:length(levels(data1)))  
{  cells <- WhichCells(object = data1, idents = levels(data1)[i])
if (length(cells) == 1)
  data.dims[,i] <- (embeddings[cells,]) else
data.dims[,i] <- colMeans(embeddings[cells,])
}
library(ape)
colnames(x = data.dims) <- levels(x = data1)
data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
nj.tree <- ape::nj(data.dist)
nj.boot.tree <- ape::boot.phylo(nj.tree, t(data.dims), FUN = function(xx) nj(dist(xx)), B = 100)

# add bootstraps:
nj.tree$node.label <- nj.boot.tree/1
nj.tree$node.label <- round(nj.tree$node.label)

node_col <- nj.tree$node.label
node_col[nj.tree$node.label < 80] <- "slateblue"
node_col[nj.tree$node.label < 50] <- "red3"
node_col[nj.tree$node.label >= 80] <- "black"
nj.tree$edge.length <- sqrt(nj.tree$edge.length)

# plot tree

plot.phylo(nj.tree, type = "phylogram", label.offset = 0.5,show.node.label = T,use.edge.length = F,lab4ut = 'axial', cex = 0.5,no.margin = T)
# for (i in 1:length(node_col)) {
#   nodelabels(node = length(nj.tree$tip.label)+i, pch=21, col="black", bg=node_col[i], cex=2)
# }

ape::plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = T,use.edge.length = F,lab4ut = 'axial', cex = 1,no.margin = T, rotate.tree = 120)
for (i in 1:length(node_col)) {
  nodelabels(node = length(nj.tree$tip.label)+i, pch=21, col="black", bg=node_col[i], cex=2)
}

## neuron top go ----
library(topGO)

Unfold <- annot.red.Plus %>%
dplyr::mutate(`go` = strsplit(as.character(`go`), ",")) %>%
tidyr::unnest(`go`)
geneID2GO <- Unfold %>% split(x = .$`go`, f = .$gene_short_name.x)
geneNames <- names(geneID2GO)

all.markers <- FindAllMarkers(data1.subsets$neural,assay='RNA',  return.thresh = 0.01,  min.pct = 0.2,  only.pos = TRUE)

GOdata = NULL
myInterestingGenes = NULL
geneList = NULL
# this takes a really long time:
for(i in 1:length(unique(all.markers$cluster)))
{
#identify the GOI from the DE list
myInterestingGenes[[i]] <- all.markers$gene[all.markers$cluster==unique(all.markers$cluster)[i]] #list of genes.ac you want to perform GO enrichment for
geneList[[i]] <- factor(as.integer(geneNames %in% myInterestingGenes[[i]]))
names(geneList[[i]]) <- geneNames
GOdata[[i]] <- new("topGOdata",ontology = "BP", allGenes = geneList[[i]],annot = annFUN.gene2GO, gene2GO = geneID2GO)
}
# filter and generate figures.
allRes_intgenes = NULL
resultFis=NULL
topGO.IDs=NULL
for(i in 1:length(unique(all.markers$cluster)))
{
#run the test...
resultFis[[i]] <- runTest(GOdata[[i]], algorithm = "classic", statistic = "fisher")
allRes_intgenes[[i]]<- GenTable(GOdata[[i]], pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=50)
allRes_intgenes[[i]]$pvalues<-as.numeric(allRes_intgenes[[i]]$pvalues)
#convert NA values to zero; only if p value is so small it cannot be displayed in r
allRes_intgenes[[i]][is.na(allRes_intgenes[[i]])]<-0.00000001
}
#plot GOenrichment
for(i in 1:length(unique(all.markers$cluster)))
{
topGO.IDs[[i]]=
ggplot2::ggplot(allRes_intgenes[[i]][1:10,], #pulling out only the top 10 terms for plotting
                ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
col=neur.cp[i],fill=neur.cp[i]) +
ggplot2::coord_flip()+
ggplot2::xlab("Biological Process") +
ggplot2::ylab("Enrichment -log10 p-value") +
ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
}

#save for later:
save(all.markers,topGO.IDs,allRes_intgenes,resultFis,file='topGO.neurons.RData')

## generate muscle gene lists (DataS1.13) ----
Ac.Alldata<-SetIdent(Ac.Alldata,value='IDs')
coi.st=WhichCells(Ac.Alldata,idents ='muscle.st')

Ac.Alldata<-SetIdent(Ac.Alldata,value='ID.separate')
coi.sm=WhichCells(Ac.Alldata,idents ='smooth.muscle')
data1=Ac.Alldata
Idents(data1) <- data1$IDs
data1$reduced=as.character(data1@active.ident)

data1$reduced[coi.sm]=as.character(Ac.Alldata$ID.separate[coi.sm])
data1<-SetIdent(data1,value='reduced')
Ac.Alldata$reduced=data1@active.ident
#' **update** the script... save into AcAlldata...
clust.cp.reduced=c('goldenrod2','#FB6496','seagreen3','purple3','grey','black',"#ff7f00",'red2','#1F78C8')
names(clust.cp.reduced)=levels(SetIdent(Ac.Alldata,value='reduced'))

levels(data1)
data1@active.assay='RNA'

data1<-FindVariableFeatures(data1,nfeatures = 6000,verbose=F)
goi=data1@assays$RNA@var.features
#use dotplot to generate pct exp data:
expression=DotPlot(data1,'RNA',goi,group.by = 'reduced')
expression.data=expression$data
expression.st=expression.data[expression.data$id=='muscle.st',]
expression.sm=expression.data[expression.data$id=='smooth.muscle',]
expression.other=expression.data[!expression.data$id=='muscle.st' & !expression.data$id=='smooth.muscle',]

goi.st=goi[which(rowSums(data1@assays$RNA@counts[goi,WhichCells(data1,idents = 'muscle.st')])>=3)] #list of all genes detectable in the st.muscle

goi.sm=goi[which(rowSums(data1@assays$RNA@counts[goi,WhichCells(data1,idents = 'smooth.muscle')])>=3)] #list of all genes detectable in the sm.muscle

# use the dot plot data to filter the genes:

# filter for at least 10% of the cells in either cluster
muscle.st=expression.data[(expression.data$pct.exp >= 10 & expression.data$avg.exp.scaled >0 & expression.data$id == 'muscle.st'),]
muscle.sm=expression.data[(expression.data$pct.exp >= 10 & expression.data$avg.exp.scaled >0 & expression.data$id == 'smooth.muscle'),]

#generate a specificity profile: not expressed in other clusters:

specificity.index.st=  expression.other[expression.other$features.plot %in% muscle.st$features.plot,] %>% group_by(features.plot) %>% count(pct.exp <=10 & avg.exp.scaled <=0) %>% arrange(desc(n))  
muscle.st.specific=as.character(specificity.index.st$features.plot[specificity.index.st$n>=7&specificity.index.st[,2]==TRUE])

specificity.index.sm=  expression.other[expression.other$features.plot %in% muscle.sm$features.plot,] %>% group_by(features.plot) %>% count(pct.exp <=10 & avg.exp.scaled <=0) %>% arrange(desc(n))  
muscle.sm.specific=as.character(specificity.index.sm$features.plot[specificity.index.sm$n>=7&specificity.index.sm[,2]==TRUE])

muscle.all=intersect(muscle.sm.specific,muscle.st.specific)

goi=unique(c(muscle.all,muscle.sm.specific,muscle.st.specific))
#check these lists:
DotPlot(data1,'RNA',features=goi,cols=c('lightgrey','red'),col.min = 0)&RotatedAxis()&theme(panel.grid =element_line(colour = 'grey90',linewidth = 0.25))

save(data1,muscle.all,muscle.sm.specific,muscle.st.specific,file='muscle.profiles.RData')

# Muscle Genes ----
DataS2.sh=genes.ac[match(muscle.all,genes.ac$gene_short_name),] 
DataS2.sm=genes.ac[match(setdiff(muscle.sm.specific,muscle.all),genes.ac$gene_short_name),] 
DataS2.st=genes.ac[match(setdiff(muscle.st.specific,muscle.all),genes.ac$gene_short_name),] 

DataS2.sh$type = 'shared'
DataS2.sm$type = 'smooth'
DataS2.st$type = 'striated'

save (DataS2.sh,DataS2.sm,DataS2.st,file='DataS2.RData')
xlsx::write.xlsx(DataS2.st,'DataS1.13.xlsx',sheetName='striated.specific',col.names = T,row.names = F)
xlsx::write.xlsx(DataS2.sm,'DataS1.13.xlsx',sheetName='smooth.specific',col.names = T,row.names = F,append=T)
xlsx::write.xlsx(DataS2.sh,'DataS1.13.xlsx',sheetName='shared',col.names = T,row.names = F,append=T)

## DEGs per cell state per lifehistory: ----

{data1<-Ac.Alldata
data1$split.ID = paste (data1$IDs,data1$lifehistory)
data1$split.ID.separate = paste (data1$ID.separate,data1$lifehistory)
data1$split.ID.separate=as.factor(data1$split.ID.separate)
order = sort(levels(data1$split.ID.separate),index.return=T)
data1$split.ID.separate=factor(data1$split.ID.separate,levels = levels(data1$split.ID.separate)[order$ix])
data1$split.ID.separate=droplevels(data1$split.ID.separate)
data1<-SetIdent(data1,value='split.ID.separate')

data1@active.assay = 'RNA' # set RNA assay
all.markers <- FindAllMarkers(
  data1,
  logfc.threshold = 1,
  return.thresh = 0.001,
  min.pct = 0.5,
  only.pos = TRUE,
)

all.markers<-update_marker_list(all.markers)
# View(all.markers)
goi=get_marker_list(all.markers,2)

DotPlot(data1,assay='RNA', features = unique(rev(goi)),
        scale.by='radius' , col.min = 0, col.max = 4,
        cols = c('lightgrey','darkred'))+  RotatedAxis() +FontSize(6,8) +labs(title ='DETFs')+theme(panel.grid = element_line(color='black',linewidth = 0.2))
}
save(all.markers,file='DEG_allIDsbystage.RObj')
