
library(Matrix)
library(Seurat,quietly=T)
packageVersion('Seurat') #check that it switched!
library(RColorBrewer,quietly=T)
library(patchwork,quietly=T)
library(ggplot2,quietly=T)
library(pals,quietly=T) #had to add this
library(readxl,quietly=T)
library(tidyr,quietly=T)
library(dplyr,quietly=T)

load (file="Aaur2.newnames.RData")

LibCP =brewer.paired(12)
LibCP.stages=LibCP[c(2,4,5,8)]
gene.cp = c('lightgrey', rev(brewer.pal(11 , "Spectral")))
# maximally unique colour palette for clustering:
clust.cp.separate = unique (c(cols25(25), glasbey(32), alphabet(26)))
# gradient colour palette for clustering:
clust.cp.graded = unique(c(stepped3(20), stepped2(20), stepped(20)))
do.over=F
# subset analysis ----
if (do.over)
{ 
load (file='Ac_manuscript_final/AaAlldata.Robj')
data1<-SetIdent(Aa.Alldata,value='IDs')
DimPlot(data1,cols=clust.cp.separate)
data1.subsets <- SplitObject(Aa.Alldata,split.by = 'IDs')

order = match(levels(Aa.Alldata),names(data1.subsets))
data1.subsets=data1.subsets[order]
{
  for (i in 1:length(names (data1.subsets)))
  {
    data1=data1.subsets[[i]]
    data1@active.assay = 'RNA'
    
    data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
    data1@active.assay = 'integrated'
    
    # #use the full dataset scaling:
    Aa.Alldata@active.assay='integrated'
    coi = colnames(data1)
    t=ScaleData(Aa.Alldata,features = data1@assays$RNA@var.features)
    
    data1@assays$integrated@scale.data = t@assays$integrated@scale.data[,coi]
    rm(t)
    data1 <- RunPCA(data1, pcs.compute = 50,verbose = F) #this is now failing...
    # ElbowPlot(object = data1, ndims = 50)
    # Determine percent of variation associated with each PC
    pct <- data1[["pca"]]@stdev / sum(data1[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC: last point where change of % of variation is more than 0.1%
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    # choose the minimum of these two metrics as the PCs covering the majority of the variation in the data
    pcs <- min(co1, co2)
    d = 1:pcs
    data1 <- RunUMAP(data1, reduction ="pca", n.neighbors = 10L,spread =1,seed.use = 0, 
                     dims = d,reduction.name ='umap',
                     reduction.key ='umap',min.dist = 0.3,
                     local.connectivity = 10,verbose = F)#
    data1 <- FindNeighbors(object = data1,reduction ="pca",dims = d,
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
    data1@reductions=data1@reductions[c(1,2,5,4,3)]
    data1.subsets[[i]] = data1
  }#this sets coarse clustering as default
}
# DimPlot(data1,cols=clust.cp.separate)
##refine and label ----
####outer ----
{
  data1=data1.subsets$outer.epidermis
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'OuterClusters')
  
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1@active.assay='RNA'
  #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi,verbose = F) 
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName)#for long lists you might want to know what is missing
  clusterNames$ID[clName]
  # DimPlot(data1,cols=clust.cp.separate,label=T,reduction = 'umap')&NoAxes()
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='integrated'
  print(DimPlot(data1,group.by = c('IDs'),label=F,cols = (clust.cp.separate))+NoAxes())
  # data1@reductions=data1@reductions[c(1,3,2)]
  data1.subsets$outer.epidermis=data1
}
####inner ----
{
  data1=data1.subsets$inner.gastrodermis
  data1<-SetIdent(data1,value='IDs.coarse')
  # DimPlot(data1, cols = clust.cp.separate,label = T,reduction='umap')+
    NoAxes()
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'InnerClusters') 
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1@active.assay='RNA'
  #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi, split.by = 'orig.ident',verbose=F) #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
  
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName)#for long lists you might want to know what is missing
  clusterNames$ID[clName]
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  # data1@reductions=data1@reductions[c(1,3,2)]
  print(DimPlot(data1,group.by = c('IDs'),label=F,cols = (clust.cp.separate),reduction = 'umap')+NoAxes())
  
  data1.subsets$inner.gastrodermis=data1
}
####cnido ----
{
  
  data1=data1.subsets$cnidocyte #
  data1<-SetIdent(data1,value = 'IDs.coarse')
  
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'CnidoClusters')
  
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1@active.assay='RNA'
  #assign cluster ID to the individual libraries
  
  data1<-ScaleData(data1,features = goi, split.by='orig.ident',verbose=F) #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data',verbose=F,) #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  
  sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
  clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
  # plot(DimPlot(data1,label = T,cols=clust.cp.separate,reduction='umap')+NoAxes())
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident

  data1@active.assay='integrated'
  data1.subsets$cnidocyte=data1
  # print(DimPlot(data1,cols=clust.cp.separate,reduction='umap',group.by = 'IDs')+NoAxes())#+
}
####mucin ----
{
  data1=data1.subsets$gland.muc
  data1<-SetIdent(data1,value = 'IDs.coarse')#fine picks up small clusters but don't know what to do with these
  # print(DimPlot(data1,cols = clust.cp.separate,label = T)&NoAxes())
  # data1@reductions <- data1@reductions[c(1,3,2)]
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'MucinClusters')
  
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1@active.assay='RNA'
  #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi, split.by = 'orig.ident',verbose=F) #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
  
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
  clusterNames$ID[clName] #you can also look at the names to see where the problem might have been
  
  # DimPlot(data1,label = T,cols=clust.cp.separate,reduction='umap')+NoAxes()
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='integrated'
  data1.subsets$gland.muc =data1
  
  # print(DimPlot(data1,group.by = c('IDs'),label=F,cols = (clust.cp.separate))+NoAxes())
  
  
}
####dig.gland ----
{
  data1=data1.subsets$gland.dig
  data1<-SetIdent(data1,value = 'IDs.coarse')        
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'Dig.Gland')
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName)#for long lists you might want to know what is missing
  clusterNames$ID[clName]
  # DimPlot(data1,cols=clust.cp.separate)
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='integrated'
  # data1@reductions=data1@reductions[c(1,3,2)]
  # print(DimPlot(data1,group.by = c('IDs'),label=F,reduction = 'umap',
                # cols = (clust.cp.separate))+NoAxes())
  
  data1.subsets$gland.dig=data1
  
}
####neuron.1----
{
  data1=data1.subsets$neural.1
  clusterNames<- read_excel('DataS3.xlsx', sheet = 'Neuro')

  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1<-SetIdent(data1,value = 'IDs.fine')
  # print(DimPlot(data1,cols=clust.cp.separate,reduction='umap'))
  data1@active.assay='RNA'
  
  #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi,verbose=F)
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName)#for long lists you might want to know what is missing
  clusterNames$ID[clName]
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='integrated'
  # data1@reductions=data1@reductions[c(1,3,2)]
  data1.subsets$neural.1=data1
  
  # print(DimPlot(data1,group.by = c('IDs'),label=F,reduction = 'umap',
                # cols = (clust.cp.separate))+NoAxes()+labs(title = 'neural.1 population'))
}
####neuron.2 ----
{
  data1=data1.subsets$neural.2
  clusterNames<- read_excel('DataS3.xlsx', sheet = 'Neuro')
  
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1<-SetIdent(data1,value = 'IDs.fine')
  # print(DimPlot(data1,cols=clust.cp.separate,reduction='umap'))
  data1@active.assay='RNA'
  
  #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi,verbose=F)
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName)#for long lists you might want to know what is missing
  clusterNames$ID[clName]
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='integrated'
  # data1@reductions=data1@reductions[c(1,3,2)]
  data1.subsets$neural.2=data1
  
  # print(DimPlot(data1,group.by = c('IDs'),label=F,reduction = 'umap',
                # cols = (clust.cp.separate))+NoAxes()+labs(title = 'neural.2 partition'))
}
####st.muscle ----
{
  data1=data1.subsets$muscle.st
  #do-over clustering
  data1<-SetIdent(data1,value = 'IDs.coarse')
  data1 <- FindNeighbors(object = data1,reduction ="pca",dims = 1:16,nn.method = 'annoy',  annoy.metric = 'cosine',verbose = F,k.param = 20)
  data1 <- FindClusters(object = data1,resolution = 0.4,random.seed = 0,verbose = F)
  # # #look at relationship between clusters
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reorder.numeric = TRUE,verbose = F)#, dims = c(1:d))
  # print(DimPlot(data1,cols=clust.cp.separate,reduction='umap')+DimPlot(data1,cols=LibCP.stages,group.by = 'lifehistory'))
  
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'Muscle.ST') 
  
  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  data1@active.assay='RNA'
  
  #assign cluster ID to the individual libraries
  data1<-ScaleData(data1,features = goi,verbose=F)#, split.by = 'orig.ident') #works with scaled gene values, so must scale all of these first (in case they are not within the VarFeatures set!)
  
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName)#for long lists you might want to know what is missing
  clusterNames$ID[clName]
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay='integrated'
  # data1@reductions=data1@reductions[c(1,3,2)]
  # print(DimPlot(data1,group.by = c('IDs'),label=F,reduction = 'umap',
  #               cols = (clust.cp.separate))+NoAxes())
  # 
  data1.subsets[[1]]=data1
}

save(data1.subsets,file='Ac_manuscript_final/Aa.subsets.RObj')

{
  cl.order=NULL
  cl.ind = NULL
  # Aa.Alldata$IDs = Aa.Alldata@active.ident 
  Aa.Alldata$ID.separate = as.character(Aa.Alldata$IDs)
  for (i in 1:length(names (data1.subsets)))
  {
    coi=NULL
    coi=colnames(data1.subsets[[i]])
    Aa.Alldata$ID.separate[coi] = as.character(data1.subsets[[i]]@active.ident[coi])
    cl.order = c(cl.order,levels(data1.subsets[[i]]))
  }
  cl.ind = match(cl.order,levels(as.factor(Aa.Alldata$ID.separate))) 
  
  Aa.Alldata$ID.separate = as.factor(Aa.Alldata$ID.separate)
  Aa.Alldata$ID.separate = factor(Aa.Alldata$ID.separate,levels(Aa.Alldata$ID.separate)[cl.ind])
  
}
save(Aa.Alldata,file='Ac_manuscript_final/AaAlldata.Robj')
###all.neurons ----
{
  data1=merge(data1.subsets$neural.1,data1.subsets$neural.2,merge.data=T)
  data1@active.assay = 'RNA'
  data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
  
  # #use the full dataset scaling:
  Aa.Alldata@active.assay='integrated'
  coi = colnames(data1)
  t=ScaleData(Aa.Alldata,features = data1@assays$RNA@var.features)
  
  data1@assays$integrated@scale.data = t@assays$integrated@scale.data[,coi]
  rm(t)
  data1@active.assay = 'integrated'
  
  data1 <- RunPCA(data1, pcs.compute = 50,verbose = F,features = data1@assays$RNA@var.features) #this is now failing...
  # ElbowPlot(object = data1, ndims = 50)
  d = 1:20
  data1 <- RunUMAP(data1, reduction ="pca", n.neighbors = 10L,spread =0.6,seed.use = 0, 
                   dims = d,reduction.name ='umap',
                   reduction.key ='umap',min.dist = 0.2,
                   local.connectivity = 1,verbose = F)#
  
  # DimPlot(data1,cols=clust.cp.graded)
  all.neurons=data1
  save(all.neurons,file='Ac_manuscript_final/all.neurons.2.Robj')
  # data1=all.neurons
}
### all.muscle subset ----
{
allmuscle=subset(Aa.Alldata,idents=levels(Aa.Alldata)[c(1,2)])
data1=SetIdent(allmuscle,value = 'ID.separate')
data1@active.assay='RNA'
data1 <- FindVariableFeatures(data1,nfeatures = 2000,verbose = F)#
# #use the full dataset scaling:
Aa.Alldata@active.assay='integrated'
coi = colnames(data1)
t=ScaleData(Aa.Alldata,features = data1@assays$RNA@var.features)
data1@active.assay='integrated'
data1@assays$integrated@scale.data = t@assays$integrated@scale.data[,coi]
rm(t)
data1 <- RunPCA(data1, pcs.compute = 50,verbose = F) #this is now failing...
d = 1:20
data1 <- RunUMAP(data1, reduction ="pca", n.neighbors = 10L,spread =1,seed.use = 42, 
                 dims = d,reduction.name ='umap',
                 reduction.key ='umap',min.dist = 0.6,
                 local.connectivity = 1,verbose = F)#

allmuscle=data1
save(allmuscle,file='Ac_manuscript_final/New.allmuscle.RObj')

coi.st=WhichCells(Aa.Alldata,idents ='muscle.st')
coi.sm=WhichCells(allmuscle,idents ='outer.smooth.muscle')


data1=Aa.Alldata
Idents(data1) <- data1$IDs
data1$reduced=as.character(data1@active.ident)
data1$reduced[coi.sm]=as.character(Aa.Alldata$ID.separate[coi.sm])
data1<-SetIdent(data1,value='reduced')
data1@active.assay='RNA'

data1<-FindVariableFeatures(data1,nfeatures = 4000,verbose=F)
goi=data1@assays$RNA@var.features
#use dotplot to generate pct exp data:
expression=DotPlot(data1,'RNA',goi,group.by = 'reduced')
expression.data=expression$data
expression.st=expression.data[expression.data$id=='muscle.st',]
expression.sm=expression.data[expression.data$id=='outer.smooth.muscle',]
expression.other=expression.data[!expression.data$id=='muscle.st' & !expression.data$id=='outer.smooth.muscle',]

goi.st=goi[which(rowSums(data1@assays$RNA@counts[goi,WhichCells(data1,idents = levels(data1)[9])])>=3)] #list of all genes detectable in the st.muscle

goi.sm=goi[which(rowSums(data1@assays$RNA@counts[goi,WhichCells(data1,idents = levels(data1)[3])])>=3)] #list of all genes detectable in the sm.muscle

# use the dot plot data to filter the genes:

# filter for at least 10% of the cells in either cluster
muscle.st=expression.data[(expression.data$pct.exp >= 10 & expression.data$avg.exp.scaled >0 & expression.data$id == 'muscle.st'),]
muscle.sm=expression.data[(expression.data$pct.exp >= 10 & expression.data$avg.exp.scaled >0 & expression.data$id == 'outer.smooth.muscle'),]

#check these lists:
# DotPlot(data1,'RNA',features=as.character(unique(muscle.sm$features.plot)),cols=c('lightgrey','red'),col.min = 0)&RotatedAxis()&theme(panel.grid =element_line(colour = 'grey90',linewidth = 0.25))
# DotPlot(data1,'RNA',features=as.character(unique(muscle.st$features.plot)),cols=c('lightgrey','red'),col.min = 0)&RotatedAxis()&theme(panel.grid =element_line(colour = 'grey90',linewidth = 0.25))

#generate a specificity profile: not expressed in other clusters:
library(dplyr)
specificity.index.st=  expression.other[expression.other$features.plot %in% muscle.st$features.plot,] %>% group_by(features.plot) %>% count(pct.exp <=40 & avg.exp.scaled <=0.5) %>% arrange(desc(n))  
muscle.st.specific=as.character(specificity.index.st$features.plot[specificity.index.st$n>=7&specificity.index.st[,2]==TRUE])

specificity.index.sm=  expression.other[expression.other$features.plot %in% muscle.sm$features.plot,] %>% group_by(features.plot) %>% count(pct.exp <=40 & avg.exp.scaled <=0.5) %>% arrange(desc(n))  
muscle.sm.specific=as.character(specificity.index.sm$features.plot[specificity.index.sm$n>=7&specificity.index.sm[,2]==TRUE])

muscle.all=intersect(muscle.sm.specific,muscle.st.specific)

goi=unique(c(muscle.all,muscle.sm.specific,muscle.st.specific))
save(data1,muscle.all,muscle.sm.specific,muscle.st.specific,file='muscle.profiles.RData')

# DataS2 ----
DataS2.sh=Aa_genes[match(muscle.all,Aa_genes$gene_short_name),3:8] 
DataS2.sm=Aa_genes[match(setdiff(muscle.sm.specific,muscle.all),Aa_genes$gene_short_name),3:8] 
DataS2.st=Aa_genes[match(setdiff(muscle.st.specific,muscle.all),Aa_genes$gene_short_name),3:8] 

xlsx::write.xlsx(DataS2.st,'Ac_manuscript_final/DataS2.xlsx',sheetName='striated.specific',col.names = T,row.names = F)
xlsx::write.xlsx(DataS2.sm,'Ac_manuscript_final/DataS2.xlsx',sheetName='smooth.specific',col.names = T,row.names = F,append=T)
xlsx::write.xlsx(DataS2.sh,'Ac_manuscript_final/DataS2.xlsx',sheetName='shared',col.names = T,row.names = F,append=T)
}


#calculate medus-only ----
{
  data1=SetIdent(Aa.Alldata,value='ID.separate')
  
  #  use the idents table:
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory))
  colnames(ids.cluster.library) = c('ID','Lifehistory','CellCount')
 
  medusa.cl=as.character(ids.cluster.library$ID[which(ids.cluster.library$CellCount[ids.cluster.library$Lifehistory=='polyp'] < 1)])
  
  polyp.cl=as.character(ids.cluster.library$ID[which(ids.cluster.library$CellCount[ids.cluster.library$Lifehistory=='medusa'] <=4)])
}
}else{
  load (file='Ac_manuscript_final/Aa.subsets.RObj')
  load (file='Ac_manuscript_final/AaAlldata.Robj')}
##all.idents.DEGs ----
data1=SetIdent(Aa.Alldata,value = 'ID.separate')
data1@active.assay='RNA'
all.markers.allidents <- FindAllMarkers(data1,
                                        logfc.threshold = 1,
                                        return.thresh = 0.0001,
                                        min.pct = 0.2,
                                        only.pos = TRUE, 
                                        verbose=F)

all.markers.allidents[,8:11] = 'NA'
names(all.markers.allidents)[8:11]=names(Aa_genes)[c(3,6:8)]
# add GO terms associated with this list:
for (j in 1:length(levels(data1@active.ident))) # 
{
  x=all.markers.allidents[as.numeric(all.markers.allidents$cluster)==j,][1:length(which(as.numeric(all.markers.allidents$cluster)==j)),7]
  x2=Aa_genes[match(x,Aa_genes$gene_short_name),c(3,6:8)]
  all.markers.allidents[as.numeric(all.markers.allidents$cluster)==j,][1:length(which(as.numeric(all.markers.allidents$cluster)==j)),8:12]<-x2
}  

#also for transcription factors only:
all.markers.TF.allidents <- FindAllMarkers(data1,
                                           logfc.threshold = 0.2,
                                           features = intersect(AaTF_list,rownames(data1)),
                                           min.pct = 0.05,
                                           only.pos = TRUE,
                                           return.thresh = 0.001,verbose=F)

all.markers.TF.allidents[,8:11] = 'NA'
names(all.markers.TF.allidents)[8:11]=names(Aa_genes)[c(3,6:8)]
# add GO terms associated with this list:
for (j in 1:length(levels(data1@active.ident))) #
{
  x=all.markers.TF.allidents[as.numeric(all.markers.TF.allidents$cluster)==j,][1:length(which(as.numeric(all.markers.TF.allidents$cluster)==j)),7]
  x2=Aa_genes[match(x,Aa_genes$gene_short_name),c(3,6:8)]
  all.markers.TF.allidents[as.numeric(all.markers.TF.allidents$cluster)==j,][1:length(which(as.numeric(all.markers.TF.allidents$cluster)==j)),8:12]<-x2
}    
save (all.markers.TF.allidents,all.markers.allidents,file='Ac_manuscript_final/all.idents.markers.RData')

list = NULL
for (c in  1:length(levels(data1@active.ident)))
{
  x=all.markers.allidents[as.numeric(all.markers.allidents$cluster)==c,][1:min(5,length(which(as.numeric(all.markers.allidents$cluster)==c))),7]
  list=c(list,x)
}
list2 = NULL
for (c in  1:length(levels(data1@active.ident)))
{
  x=all.markers.TF.allidents[as.numeric(all.markers.TF.allidents$cluster)==c,][1:min(5,length(which(as.numeric(all.markers.TF.allidents$cluster)==c))),7]
  list2=c(list2,x)
}

# print(DotPlot(data1,'RNA', features = unique(list),
#               scale.by='radius' , col.min = 0, col.max = 4,
#               cols = c('lightgrey','darkred'))+
#         RotatedAxis() +FontSize(6,8) +labs(title='DEG')+theme(legend.position = 'bottom',legend.title = element_text(size=8)))+coord_flip()
# print(
#   DotPlot(data1,'RNA', features = unique(list2),
#           scale.by='radius' , col.min = 0, col.max = 4,
#           cols = c('lightgrey','darkred'))+
#     RotatedAxis() +FontSize(6,8) +labs(title ='DETFs')+theme(axis.text.y = element_text(colour=all.clusters.cp,size=8)))

# Generate gene lists subsets ----
{
  markers = NULL
  markers.TF = NULL# 
  for (i in 1:length(names (data1.subsets)))
  {
    data1=data1.subsets[[i]]
    data1@active.assay='RNA'
    all.markers <- FindAllMarkers(data1,
                                  logfc.threshold = 0.6,
                                  return.thresh = 0.001,
                                  min.pct = 0.2,
                                  only.pos = TRUE)

    all.markers[,8:11] = 'NA'
    names(all.markers)[8:11]=names(Aa_genes)[c(3,6,7,8)]
    # add GO terms associated with this list:
    for (j in 1:length(levels(data1@active.ident))) #
    {
      x=all.markers[as.numeric(all.markers$cluster)==j,][1:length(which(as.numeric(all.markers$cluster)==j)),7]
      x2=Aa_genes[match(x,Aa_genes$gene_short_name),c(3,6,7,8)]
      all.markers[as.numeric(all.markers$cluster)==j,][1:length(which(as.numeric(all.markers$cluster)==j)),8:11]<-x2
    }
    markers[[i]] = all.markers
    
    #also for transcription factors only: this did not work...
    all.markers_TF <- FindAllMarkers(data1,
                                     features = intersect(AaTF_list,rownames(data1)),
                                     min.pct = 0.01,
                                     only.pos = TRUE,
                                     return.thresh = 0.001,verbose=F)
    
    all.markers_TF[,8:11] = 'NA'
    names(all.markers_TF)[8:11]=names(Aa_genes)[c(3,6,7,8)]
    # add GO terms associated with this list:
    for (j in 1:length(levels(data1@active.ident))) #
    {
      x=all.markers_TF[as.numeric(all.markers_TF$cluster)==j,][1:length(which(as.numeric(all.markers_TF$cluster)==j)),7]
      x2=Aa_genes[match(x,Aa_genes$gene_short_name),c(3,6,7,8)]
      all.markers_TF[as.numeric(all.markers_TF$cluster)==j,][1:length(which(as.numeric(all.markers_TF$cluster)==j)),8:12]<-x2
    }
    markers.TF[[i]] = all.markers_TF
  }
  # View(markers[[i]])
  names(markers) = names(data1.subsets)
  names(markers.TF) = names(data1.subsets)
  
  save(markers,file = 'Ac_manuscript_final/AllData.integrated.subsets.markers.RObj')
  save(markers.TF,file = 'Ac_manuscript_final/AllData.integrated.subsets.TF.markers.RObj')
  for (i in 1:8)
  {
    write.csv(markers[[i]],file=paste('Ac_manuscript_final/DS1.',paste(i+4,".csv",sep=''),sep = ''))
    write.csv(markers.TF[[i]],file=paste('Ac_manuscript_final/DS1.',paste(i+4,"b.TF.csv",sep=''),sep = ''))
    
  }
  
}
save(Aa.Alldata,file='Ac_manuscript_final/AaAlldata.Robj')
