#OMA
.libPaths(c("/lisc/data/scratch/molevo/agcole/R/libs/seurat4/","/lisc/data/scratch/molevo/agcole/R/libs/course24/","/lisc/opt/sw/software/R/4.5.0/lib64/R/library"))
setwd("/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE_revisions/")

library(easypackages)
libraries("readxl","RColorBrewer",'ggplot2',
          'patchwork','dplyr','viridisLite','pals')
library(rliger)#,lib.loc ="/lisc/user/agcole/R/course24" )
library(SeuratObject)#,lib='/lisc/app/seurat5/5.1.0-4.4.0') 
library(Seurat)#,lib='/lisc/app/seurat4/4.4.0-4.3.3') #name updates don't work with Vs5
library(SeuratWrappers)#,lib.loc = "/lisc/app/seurat5/5.0.3-4.3.3")

## set colour palettes
gene.cp.dark=c('grey30',rev(brewer.pal(11 , "Spectral" )))
gene.cp=c('lightgrey',rev(brewer.pal(11 , "Spectral" )))
clust.cp.separate = unique (c(cols25(25),alphabet2(26),glasbey(32),alphabet(26)))
clust.cp.graded = unique(c(stepped3(16),stepped(20),stepped2(20)))

# load(file='temp.oma.RData')
######################
# load (file="/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/AcGenes.RData")
load (file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE.genes.RData')
load (file='/lisc/data/scratch/molevo/agcole/R/Nematostella_Nv2/Genes.Nv2.RData')

## OMA gene lists ----
lut=read.delim('/lisc/data/scratch/molevo/jmontenegro/alison/aaurita/results/annotation/kostya+/geneID_AAUR2.map.txt',header = F)
names(lut)=c('geneID','AAUR2')

OMA<- readxl::read_xlsx('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/OMAOrthologousGroupsAcNv.xlsx')
#Kostya Plus has different name files; First update AAUR2 with the gene_ model:
ind=match(OMA$Ac,lut$AAUR2,nomatch = 0)
OMA$Ac[1:42]=lut$geneID[ind]

ind2=match(OMA$Ac,genes.ac$cellranger,nomatch = '0')
which(ind2 == '0')
#everything is there... good.

#add gene name to the OMA file first:
OMA$nematostella.gene=OMA$Nv
ind1=match(OMA$Nv,genes$geneID,nomatch = '0')
which(ind1 == '0') #all good
OMA$nematostella.gene=genes$gene_short_name[ind1]

OMA$name=OMA$nematostella.gene
OMA$name=stringr::str_replace(OMA$name,"-like-","-OMA-")

ind2=match(OMA$Ac,genes.ac$cellranger,nomatch = '0')
which(ind2 == '0')

ind.nv=match(genes$geneID,OMA$Nv,nomatch = '0')
ind.nv=ind.nv[ind.nv>0]
genes$name = genes$gene_short_name
genes$name[ind1]=OMA$name
goi=genes$name[match(OMA$Nv,genes$geneID)]

genes.ac$name = genes.ac$gene_short_name
genes.ac$name[ind2]=OMA$name
goi.ac=genes.ac$name[match(OMA$Ac,genes.ac$geneID)]

# add ACOE column to OMA file:
OMA <- OMA %>%
  left_join(genes.ac %>% select(cellranger, ACO),
            by = c("Ac" = "cellranger"))

save (genes.ac, genes, OMA, file='OMA.genesAcNv.RData')

## cnidocytes ----
#(load(file='OMAmergeFull.Robj')
{
  # load('../../Nematostella_Nv2/data1.subsets.publish.Robj')
  # cnido.nem=data1.subsets$cnidocyte
  # save(cnido.nem,file='cnido.nv.RObj')
  load(file='../cnido.nv.RObj')
  load(file='OMA.genesAcNv.RData')
  DimPlot(cnido.nem)
temp.n=cnido.nem
rownames(temp.n@assays$RNA@counts)=genes$name
rownames(temp.n@assays$RNA@data)=genes$name
rownames(temp.n@assays$RNA@meta.features)=genes$name # otherwise fails with vs5
temp.n@assays$RNA@var.features=genes$name # otherwise fails with vs5
goi=genes$name[match(OMA$NV2,genes$geneID)]
temp.n<-ScaleData(temp.n,features = goi,split.by = 'orig.ident')
temp.n$species='nematostella'
#drop levels:
levels(temp.n)
temp.n = subset(temp.n,idents=levels(temp.n)[c(-2,-3,-4,-5)])


goi=genes.ac$name[match(OMA$Ac,genes.ac$cellranger)]

load(file='Ac.subsets.RObj')
temp=data1.subsets$cnidocyte#cnido.aur
rownames(temp@assays$RNA@counts)=genes.ac$name
rownames(temp@assays$RNA@data)=genes.ac$name
rownames(temp@assays$RNA@meta.features)=genes.ac$name # otherwise fails with vs5
temp@assays$RNA@var.features=genes.ac$name # otherwise fails with vs5
temp<-ScaleData(temp,features = goi,
                split.by = 'orig.ident')
temp$species='aurelia'

data1=merge(temp,temp.n) #good to go... run the pipeline
data1<-subset(data1,idents=levels(data1)[c(-3,-8)]) #drop medusa and aurelin
# use only the OMA:
data1=subset(data1,features = OMA$name)
data1 <- subset(x = data1, subset = nFeature_RNA > 350 & nCount_RNA > 500)
#filter out really poor-quality cells that no longer have enough information:
levels(data1)

# #put the cnidocytes in temporal order:
cl.order=c(7,3,1,2,4,6,5,14,8,11,12,13,15,16,10,9)
data1$IDs=as.factor(data1$IDs)
data1$IDs=factor(data1$IDs,levels(data1$IDs)[cl.order])
# data1$IDs=factor(data1$IDs,levels(data1$IDs)[c(1,7,8,9,2,3,10:15,4,5,17,16,18,6,19,20,21)])
data1@active.ident=data1$IDs
levels(data1)
}
## Seurat ----
{
  VlnPlot(data1, features = c('nFeature_RNA', 'nCount_RNA'))
  data1 <- NormalizeData(data1, scale.factor = 5000)
  # {
  #   list=  NULL
  #   vargenelist <- SplitObject(data1, split.by = "orig.ident")
  #   for (i in 1:length(vargenelist)) {
  #     vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 2000, verbose = FALSE)
  #   }
  #   for (i in 1:length(vargenelist)) {
  #     x <- vargenelist[[i]]@assays$RNA@var.features
  #     list=c(list,x)}
  #   list=unique(list)
  #   length(list)
  #   list=intersect(list,OMA$nematostella.gene)
  #   data1@assays$RNA@var.features = list
  #   
  # } #split by libraries failed because not enough cells from each sample
  data1 <- FindVariableFeatures(data1, nfeatures = 4000)

  data1 <- ScaleData(data1,split.by = 'orig.ident')
  data1 <- RunPCA(data1, pcs.compute = 50)
  d = 1:20
  
  ### harmony ----
  # insufficient when using all genes
 data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')
  ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
    data1 <- RunUMAP(data1, dims = 1:50,
                   reduction = 'harmony',
                   reduction.name ='umap.h',reduction.key ='umap.h', #this is default; can change name to save different maps
                   #the following parameters control your output.
                   n.neighbors = 10L, #how many similar cells do you expect
                   spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                   min.dist = 0.2, #how close to plot each cell higher=spread
                   local.connectivity = 10)#overall how connected is the graph
  DimPlot(data1, group.by = 'species',reduction = 'umap.h')+NoAxes()+labs(title='Harmony',subtitle = 'UMAP.alt | ID')#&
    DimPlot(data1, group.by = 'IDs',reduction = 'umap.h',
            cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Original IDs')
    # Use integrated data to calculate also clustering:
    data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                           # nn.method = 'annoy',  
                           #these parameters are same as UMAP uses. This is NOT seurat default
                           annoy.metric = 'cosine', k.param = 40)
    data1<- FindClusters(data1,resolution = 1)
    #calculate relationship between clusters; DISCUSS
    data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                              dims = 1:10,
                              reorder.numeric = T)
    DimPlot(data1,reduction = 'umap.h',
            cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Clusters')
}
  save(data1,file='OMAmergeCnido.Robj')
  save(temp,temp.n,data1,OMA,genes,genes.ac,file='oma.nv.ac.cnido.RData')
  
  ## Neurons ----
  {
    load('../../Nematostella_Nv2/data1.subsets.publish.Robj')
    neuro.nem=data1.subsets$neurogland.all
    # save(neuro.nem,file='neuro.nv.RObj')
    DimPlot(neuro.nem)
    temp.n=neuro.nem
    rownames(temp.n@assays$RNA@counts)=genes$name
    rownames(temp.n@assays$RNA@data)=genes$name
    rownames(temp.n@assays$RNA@meta.features)=genes$name # otherwise fails with vs5
    temp.n@assays$RNA@var.features=genes$name # otherwise fails with vs5
    goi=genes$name[match(OMA$NV2,genes$geneID)]
    temp.n<-ScaleData(temp.n,features = goi,split.by = 'orig.ident')
    temp.n$species='nematostella'
    #drop levels GD:
    levels(temp.n)
    temp.n = subset(temp.n,idents=levels(temp.n)[c(-24,-25,-26,-33,-34,-35,-36,-37,-47)])#maybe drop also S3 and 'early'... check first
    
    
    goi=genes.ac$gene_short_name[match(OMA$Ac_K,genes.ac$geneID)]

    load(file='../Ac_manuscript_revision/Ac.subsets.RObj')
    temp=data1.subsets$neural
    rownames(temp@assays$RNA@counts)=genes.ac$name
    rownames(temp@assays$RNA@data)=genes.ac$name
    rownames(temp@assays$RNA@meta.features)=genes.ac$name # otherwise fails with vs5
    temp@assays$RNA@var.features=genes.ac$name # otherwise fails with vs5
    temp<-ScaleData(temp,features = goi,split.by = 'orig.ident')
    temp$species='aurelia'
    
    data1=merge(temp,temp.n) #good to go... run the pipeline
    levels(data1)
    # data1<-subset(data1,idents=levels(data1)[c(1:7,9:17)])
    # use only the OMA:
    data1=subset(data1,features = OMA$name)
    data1 <- subset(x = data1, subset = nFeature_RNA > 350 & nCount_RNA > 500)
    #filter out really poor-quality cells that no longer have enough information:
    levels(data1)
  }
  ## Seurat ----
  {
    VlnPlot(data1, features = c('nFeature_RNA', 'nCount_RNA'))
    data1 <- NormalizeData(data1, scale.factor = 5000)
    data1 <- FindVariableFeatures(data1, nfeatures = 4000)
    
    data1 <- ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50)
    d = 1:20
    
    ### harmony ----
    # insufficient when using all genes
    data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')
    ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
    data1 <- RunUMAP(data1, dims = 1:50,
                     reduction = 'harmony',
                     reduction.name ='umap.h',reduction.key ='umap.h', #this is default; can change name to save different maps
                     #the following parameters control your output.
                     n.neighbors = 10L, #how many similar cells do you expect
                     spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                     min.dist = 0.2, #how close to plot each cell higher=spread
                     local.connectivity = 10)#overall how connected is the graph
    DimPlot(data1, group.by = 'species',reduction = 'umap.h')+NoAxes()+labs(title='Harmony',subtitle = 'UMAP.alt | ID')#&
    DimPlot(data1, group.by = 'IDs',reduction = 'umap.h',
            cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Original IDs')
    # # Use integrated data to calculate also clustering:
    # data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
    #                        # nn.method = 'annoy',  
    #                        #these parameters are same as UMAP uses. This is NOT seurat default
    #                        annoy.metric = 'cosine', k.param = 40)
    # data1<- FindClusters(data1,resolution = 1)
    # #calculate relationship between clusters; DISCUSS
    # data1 <- BuildClusterTree(object = data1, reorder = TRUE,
    #                           dims = 1:10,
    #                           reorder.numeric = T)
    # DimPlot(data1,reduction = 'umap.h',
    #         cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Clusters')
    
    #### Celltype tree
    
    data1<-SetIdent(data1,value = 'IDs')
    embeddings <- Embeddings(object = data1, reduction = 'harmony')[,1:50]
    data.dims=matrix('0',50L,length(levels(data1)))
    for (i in 1:length(levels(data1)))  
    {  cells <- WhichCells(object = data1, idents = levels(data1)[i])
    data.dims[,i] <- colMeans(embeddings[cells,])
    }
    library(ape)
    colnames(x = data.dims) <- levels(x = data1)
    data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
    nj.tree <- ape::nj(data.dist)
    nj.boot.tree <- ape::boot.phylo(nj.tree, t(data.dims), FUN = function(xx) nj(dist(xx)), B = 1000)
    
    # add bootstraps:
    nj.tree$node.label <- nj.boot.tree/10
    nj.tree$node.label <- round(nj.tree$node.label)
    
    node_col <- nj.tree$node.label
    node_col[nj.tree$node.label < 80] <- "slateblue"
    node_col[nj.tree$node.label < 50] <- "red3"
    node_col[nj.tree$node.label >= 80] <- "black"
    nj.tree$edge.length <- sqrt(nj.tree$edge.length)
    
    # plot tree
    
    plot.phylo(nj.tree, type = "tidy", label.offset = 0.2,show.node.label = T,use.edge.length = F,lab4ut = 'axial', cex = 1,no.margin = T, rotate.tree = 90)
    # ape::plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = T,use.edge.length = F,lab4ut = 'axial', cex = 1,no.margin = T, rotate.tree = 120)
    for (i in 1:length(node_col)) {
      nodelabels(node = length(nj.tree$tip.label)+i, pch=21, col="black", bg=node_col[i], cex=2)
    }
    #NOT ENOUGH INFO in the OMA 1:1...
  }
  save(data1,file='OMAmergeNeuro.Robj')
  
  ## 4.1 plots.summary ----
  #* generates UMAP plots and bar charts
  # set a maximally unique colour palette for clustering:
  clust.cp = clust.cp.separate
  {
    #summarize the data of interest:
    data1$orig.ident = droplevels(as.factor(data1$orig.ident))
    ids.cluster.sample = as.data.frame(table(Idents(data1), data1$IDs))
    colnames(ids.cluster.sample) = c('ID', 'sample', 'CellCount')
    
    #generate your sample and cluster UMAP plots:
    sample.plot =      DimPlot(
      data1,
      group.by = 'orig.ident',
      order = rev(levels(data1$orig.ident)),
      cols = rev(LibCP)
    ) + NoAxes() +
      labs(title = 'Time | sample origin')#+NoLegend()
    
    cluster.plot = DimPlot(
      data1,
      label = T,
      label.size = 4,
      repel = T,
      # reduction = 'umap',
      #order=rev(levels(data1@active.ident)),
      #ordering is optional
      cols = clust.cp
    ) + 
      NoAxes() +
      labs(title = 'Clusters | ID') + NoLegend()
    cluster.plot
    #barplot of cluster identities in each sample:
    
    dist.clust2 =
      ggplot(ids.cluster.sample, aes(fill = ID, y = CellCount,
                                     x = sample)) +
      geom_bar(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          x = (sample)
        ),
        position = "fill",
        stat = "identity",
        width = 0.5
      ) +
      scale_fill_manual(values = clust.cp) +
      theme(axis.text.x = element_text(
        #face="bold", color="#993333",
        size = 8,
        angle = -45,
        hjust = 0,
        vjust = 0.5
      )) +
      geom_area(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          x = as.integer(sample)
        ),
        position = "fill",
        stat = "identity",
        alpha = 0.2 ,
        size = .5,
        colour = "white"
      ) +
      geom_bar(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          #this re-plots the bars over the area
          x = (sample)
        ),
        position = "fill",
        stat = "identity",
        width = 0.5
      ) +
      ggtitle("Distribution of cell types in time and space")
    
    dist.clust =
      ggplot(ids.cluster.sample, aes(fill = ID, y = CellCount,
                                     x = sample)) +
      geom_bar(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          x = (sample)
        ),
        position = "stack",
        stat = "identity",
        width = 0.5
      ) +
      scale_fill_manual(values = clust.cp) +
      theme(axis.text.x = element_text(
        #face="bold", color="#993333",
        size = 8,
        angle = -45,
        hjust = 0,
        vjust = 0.5
      )) +
      geom_area(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          x = as.integer(sample)
        ),
        position = "stack",
        stat = "identity",
        alpha = 0.2 ,
        size = .5,
        colour = "white"
      ) +
      geom_bar(
        mapping = aes(
          fill = ID,
          y = (CellCount),
          #this re-plots the bars over the area
          x = (sample)
        ),
        position = "stack",
        stat = "identity",
        width = 0.5
      ) +
      ggtitle("Distribution of cell types in time and space")
    
    dist.lib.clust =
      ggplot(ids.cluster.sample, aes(fill = sample, y = CellCount,
                                     x = ID)) +
      geom_bar(
        mapping = aes(
          fill = sample,
          y = (CellCount),
          x = (ID)
        ),
        position = "fill",
        stat = "identity",
        width = 0.5
      ) +
      scale_fill_manual(values = LibCP) +
      theme(axis.text.x = element_text(
        #face="bold", color="#993333",
        size = 8,
        angle = -45,
        hjust = 0,
        vjust = 0.5
      )) +
      geom_area(
        mapping = aes(
          fill = sample,
          y = (CellCount),
          x = as.integer(ID)
        ),
        position = "fill",
        stat = "identity",
        alpha = 0.2 ,
        size = .5,
        colour = "white"
      ) +
      geom_bar(
        mapping = aes(
          fill = sample,
          y = (CellCount),
          #this re-plots the bars over the area
          x = (ID)
        ),
        position = "fill",
        stat = "identity",
        width = 0.5
      ) +
      ggtitle("Distribution of samples across cell states")
    
    
  }
  
  dist.lib.clust+cluster.plot
  FeaturePlot(data1,c('SoxC','SoxB.2','Myc2','myc1'),order=T,cols=gene.cp,split.by = 'species')&NoAxes()
  FeaturePlot(data1,c('c-jun','Calmodulin','CALB2-OMA-1','PKD2-OMA-2'),order=T,cols=gene.cp,split.by = 'species')&NoAxes()
  FeaturePlot(data1,c('GFI1B-like-1','GFI1B-like-1 : NV2.25410'),order=T,cols=gene.cp,split.by = 'species')&NoAxes()
  
    #run degs ----
    data1<-SetIdent(data1,value = 'tree.ident')
    all.markers <- FindAllMarkers(
      data1,
      # logfc.threshold = 1,
      features = data1@assays$RNA@var.features,
      #this is faster but restricted to variable gene set.
      return.thresh = 0.0001,
      # min.pct = 0.2,
      only.pos = TRUE,
    )
    # very nematostella heavy, no genes from aurelia included.
    list = NULL
    for (i in 1:length(levels(data1@active.ident)))
    {
      x = all.markers[as.numeric(all.markers$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers$cluster) == i))), 7]
      
      list = c(list, x)
    }
    DotPlot(
           data1, assay = 'RNA',
           features = unique(c(list)),
           scale.by = 'radius' ,
           col.min = 0,
           col.max = 2, 
           # split.by = 'species',
           cols = c('slateblue', 'darkorange'))&
           RotatedAxis()&  FontSize(6, 6)&  coord_flip()&
           labs(title = 'Top 5 DEGs', subtitle = 'p-val < 0.001')
    
  data1<-SetIdent(data1,value = 'IDs')
  all.markers.IDs <- FindAllMarkers(
    data1,
    logfc.threshold = 1,
    features = data1@assays$RNA@var.features,
    #this is faster but restricted to variable gene set.
    return.thresh = 0.001,
    min.pct = 0.2,
    only.pos = TRUE,
  )
  
  list.ids = NULL
  for (i in 1:length(levels(data1@active.ident)))
  {
    x = all.markers.IDs[as.numeric(all.markers.IDs$cluster) == i, ][1:min(5, length(which(as.numeric(all.markers.IDs$cluster) == i))), 7]
    
    list.ids = c(list.ids, x)
  }
  DotPlot(
    data1, assay = 'RNA',
    features = unique(c(list.ids)),
    scale.by = 'size' ,
    col.min = 0,
    col.max = 2, split.by = 'species',
    cols = c('slateblue', 'darkorange'))&
    RotatedAxis()&  FontSize(6, 6)&  coord_flip()&
    labs(title = 'Top 5 DEGs', subtitle = 'p-val < 0.001')
  
  #can plot a tree of cluster relationships:
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reduction = 'harmony',
dims = 1:50,
# features=OMA$nematostella.gene,
reorder.numeric = F)

  ape::plot.phylo(data1@tools$BuildClusterTree, direction = 'right')

  
  ### 3.3c CellPlots ----
  DimPlot(
    data1,
    label = T,
    label.size = 5,
    repel = T,
    cols = clust.cp.separate
  ) &
    NoLegend() & NoAxes() &
    DimPlot(
      data1,
      label = T,
      label.size = 5,
      repel = T,
      cols = clust.cp.separate,
      reduction = 'tsne'
    ) &
    NoLegend() & NoAxes()
  
  #*
    #plot to confirm all is as expected:
    DimPlot(data1,
            group.by = c('IDs'),
            label = T,
            cols = clust.cp.separate) + NoAxes()

  
# OMA reduced set: AllData NV ----
    #load Aurelia data:
    load(file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/Ac_manuscript_revision/Ac.Alldata.Robj')
    
# Reduced objects ---- 
    #* NOTE: these no longer have enough information to detect the same clusters.
    #* Especially bad for the Aurelia neurons. OK for cnidocytes; the maturing and mature profiles are present
    #* 
    oma.aur=genes.ac$gene_short_name[match(OMA$Aur2,genes.ac$geneID,nomatch = '0')]
    data1= subset(neur.aur,features = oma.aur)
    data1@active.assay='RNA'
# run Seurat:
    { 
      VlnPlot(data1, features = c('nFeature_RNA', 'nCount_RNA'))
      data1 <- NormalizeData(data1, scale.factor = 1000)
      data1 <- FindVariableFeatures(data1, nfeatures = 2000)
      
      data1 <- ScaleData(data1,split.by = 'orig.ident')
      data1 <- RunPCA(data1, pcs.compute = 50)
      print(ElbowPlot(object = data1, ndims = 50) +
              # can draw lines on graph to assist interpretation:
              geom_hline(yintercept = 2))
      pct <- data1[["pca"]]@stdev / sum(data1[["pca"]]@stdev) * 100
      # Calculate cumulative percents for each PC
      cumu <- cumsum(pct)
      # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
      co1 <- which(cumu > 90 & pct < 5)[1]
      # Determine the difference between variation of PC and subsequent PC: last point where change of % of variation is more than 0.1%
      co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
      # choose the minimum of these two metrics as the PCs covering the majority of the variation in the data
      pcs <- min(co1, co2)
      d = 1:20
      
      ### harmony ----
      # insufficient when using all genes
      data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident',)
      ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
      data1 <- RunUMAP(data1, dims = d,
                       reduction = 'harmony',
                       reduction.name ='umap.h',reduction.key ='umap.h', #this is default; can change name to save different maps
                       #the following parameters control your output.
                       n.neighbors = 20L, #how many similar cells do you expect
                       spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                       min.dist = 0.4, #how close to plot each cell higher=spread
                       local.connectivity = 10)#overall how connected is the graph
      DimPlot(data1, group.by = 'orig.ident',reduction = 'umap.h',
              cols = LibCP)+NoAxes()+labs(title='Harmony',subtitle = 'UMAP.alt | ID')&
      DimPlot(data1, group.by = 'ID.separate',reduction = 'umap.h',
              cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Original IDs')
    }
    cnido.aur.oma=data1
    oma.nem=genes$gene_short_name[match(OMA$Nv2,genes$geneID,nomatch = '0')]
    data1= subset(cnido.nem,features = oma.nem)
    # run Seurat:
    { 
      VlnPlot(data1, features = c('nFeature_RNA', 'nCount_RNA'))
      data1 <- NormalizeData(data1, scale.factor = 1000)
      data1 <- FindVariableFeatures(data1, nfeatures = 2000)
      
      data1 <- ScaleData(data1,split.by = 'orig.ident')
      data1 <- RunPCA(data1, pcs.compute = 50)
      print(ElbowPlot(object = data1, ndims = 50) +
              # can draw lines on graph to assist interpretation:
              geom_hline(yintercept = 2))
      pct <- data1[["pca"]]@stdev / sum(data1[["pca"]]@stdev) * 100
      # Calculate cumulative percents for each PC
      cumu <- cumsum(pct)
      # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
      co1 <- which(cumu > 90 & pct < 5)[1]
      # Determine the difference between variation of PC and subsequent PC: last point where change of % of variation is more than 0.1%
      co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
      # choose the minimum of these two metrics as the PCs covering the majority of the variation in the data
      pcs <- min(co1, co2)
      d = 1:20
      
      ### harmony ----
      # insufficient when using all genes
      data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')
      ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
      data1 <- RunUMAP(data1, dims = d,
                       reduction = 'harmony',
                       reduction.name ='umap.h',reduction.key ='umap.h', #this is default; can change name to save different maps
                       #the following parameters control your output.
                       n.neighbors = 10L, #how many similar cells do you expect
                       spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                       min.dist = 0.4, #how close to plot each cell higher=spread
                       local.connectivity = 1)#overall how connected is the graph
      DimPlot(data1, group.by = 'orig.ident',reduction = 'umap.h',
              cols = LibCP)+NoAxes()+labs(title='Harmony',subtitle = 'UMAP.alt | ID')#&
        DimPlot(data1, group.by = 'IDs',reduction = 'umap.h',
                cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Original IDs')
    }
        
# DEGS separate OMA ----
    load(file='../Nematostella_Nv2/Alldata.Nv2.publish.Robj')
    load (file='Aa.Alldata.newnames.RObj')
  
# generate DEG lists of alldata ID.separate
temp<-SetIdent(Alldata.Nv2,value='IDs')
rownames(temp@assays$RNA@counts)=genes$name
rownames(temp@assays$RNA@data)=genes$name
rownames(temp@assays$RNA@meta.features)=genes$name # otherwise fails with vs5
# temp@assays$RNA@var.features=genes$name # otherwise fails with vs5
temp$reduced=temp@active.ident
levels(temp$reduced)=c("pSC|PGCs", "pSC|PGCs", "neuronal", "gland.mucous", "unchar.immune", "gland.S1/S2/dig", "cnidocytes", "cnidocytes.mat", "pharyngeal", 
                           "epith.", "epith.", "retractor muscle", "gastroderm", 
                           "gastroderm")
temp$reduced=droplevels(temp$reduced)

temp<-SetIdent(temp,value='reduced')

# aurelia
temp2<-SetIdent(Aa.Alldata,value='ID.cluster')
rownames(temp2@assays$RNA@counts)=genes.ac$name
rownames(temp2@assays$RNA@data)=genes.ac$name
rownames(temp2@assays$RNA@meta.features)=genes.ac$name # otherwise fails with vs5
# temp@assays$RNA@var.features=genes$name # otherwise fails with vs5
temp2@active.assay='RNA'

levels(temp2@active.ident)=c("gland.S/dig","gland.muc",'striated muscle',"unchar.center", "gast.inner",  "neuron.N1","unchar.center",  "gland.S/dig","cnidocyte","epith.outer","cnidocyte.mat", "neuron.N2")
temp2@active.ident=factor(temp2@active.ident,levels(temp2)[c(4,6,10,2,1,7,9,8,3,5)])
# DimPlot(Aa.Alldata,group.by='orig.ident',cols=brewer.paired(10))&DarkTheme()&NoAxes()&NoLegend()

temp2@active.ident=factor(temp2@active.ident,levels(temp2)[c(6,7,2,3,5,4,9,10,8,1)])
temp@active.ident=factor(temp@active.ident,levels(temp)[c(6,7,2,5,3,10,11,9,8,4,1)])
temp2$species='aurelia'
temp$species='nematostella'

#* what if you include all OMA, also paralogs in both directions:
#* 
OMA.all<- readxl::read_xlsx('OMAnv2-ac.xlsx',sheet='nv2-ac.format')
g=unique(c(OMA.all$Nv2,OMA.all$Aur2))
g.aur=genes.ac$name[match(g,genes.ac$geneID,nomatch = 0)]
g.nem=genes$name[match(g,genes$geneID,nomatch = 0)]

data1<-subset(temp2,features = g.aur)#OMA$nematostella.gene)
dim(data1)
# data1 <- subset(x = data1, subset = nFeature_RNA > 250 & nCount_RNA > 400)
data1 <- subset(x = data1, subset = nFeature_RNA > 350 & nCount_RNA > 500)
dim(data1)

data2<-subset(temp,features = g.nem)#OMA$nematostella.gene)
dim(data2)
data2 <- subset(x = data2, subset = nFeature_RNA > 350 & nCount_RNA > 500)
dim(data2)

# DEG lists
DEG.OMA.aurelia<-FindAllMarkers(data1,logfc.threshold = 1,min.pct = 0.2,
                                only.pos = T,return.thresh = 0.01)
DEG.nv2<-FindAllMarkers(data2,logfc.threshold = 1,min.pct = 0.2,
                        only.pos = T,return.thresh = 0.01)

#generate the adjacency matrix:
M=matrix(nrow = length(unique(DEG.OMA.aurelia$cluster))+1, ncol = length(levels(data2@active.ident))+1)
for (r in 1:length(unique(DEG.OMA.aurelia$cluster)))
  for (c in 1:length(levels(data2@active.ident)))
    M[r,c]=length(intersect(DEG.nv2$gene[DEG.nv2$cluster==levels(data2@active.ident)[c]],
                            DEG.OMA.aurelia$gene[DEG.OMA.aurelia$cluster==unique(DEG.OMA.aurelia$cluster)[r]]))
rand = OMA$nematostella.gene[floor(runif(100, min=1, max=length(OMA$nematostella.gene)+1))]
# rand = OMA.TF[floor(runif(5, min=1, max=length(OMA.TF)+1))]

for (c in 1:(length(levels(data2@active.ident))+1))
  M[11,c]=length(intersect(rand,DEG.nv2$gene[DEG.nv2$cluster==levels(data2@active.ident)[c]]))

rand2 = OMA$nematostella.gene[floor(runif(100, min=1, max=length(OMA$nematostella.gene)+1))]
for (r in 1:(length(unique(DEG.OMA.aurelia$cluster))+1))
  M[r,12]=length(intersect(rand2,DEG.OMA.aurelia$gene[DEG.OMA.aurelia$cluster==unique(DEG.OMA.aurelia$cluster)[r]]))

colnames(M)=c(levels(data2),'random2')
rownames(M)=c(as.character(unique(DEG.OMA.aurelia$cluster)),'random')#[order.Aa]
pheatmap::pheatmap(M,cluster_rows = F,cluster_cols = F,scale = 'none')
image(M)
CLcp.Nv = c(clust.cp.graded[c(51,3,38,53,14,6,7,27,26,1,55)][c(6,7,2,5,3,10,11,9,8,4,1)],'black')
CLcp.Aa = c(clust.cp.graded[c(51,3,2,38,14,6,7,26,1,55)][c(6,7,2,3,5,4,9,10,8,1)],'black')
clust.cp=c((CLcp.Nv),(CLcp.Aa))
names(clust.cp) = c(colnames(M),rownames(M))
pal.bands(clust.cp)
# library(circlize)
circlize::chordDiagram(t(M),
                       symmetric = F,keep.diagonal = F, row.col = clust.cp,
                       grid.col = clust.cp,
                       # transparency = 0.8,
                       # order = c(names.Nv,'random',names.Aa),
                       link.largest.ontop = T,
                       directional=2,
                       big.gap = 10
)#[c('st-M','TRm'),c('st-M','TRm')]

cluster.plot.Nv =DimPlot(data2, cols = (CLcp.Nv))+DarkTheme()+NoAxes()+NoLegend()
  # labs(title = 'Clusters | Nv')
cluster.plot.Aa =DimPlot(data1, cols = CLcp.Aa,reduction = 'umap.int')+DarkTheme()+NoAxes()+NoLegend()
  # labs(title = 'Clusters | Aa')
cluster.plot.Aa+cluster.plot.Nv#+plot_layout(ncol = 1)

#still works ok. What if you merge these two?
data1.aur.all=data1

data1=merge(temp,temp2) # this is OK, but needs to maybe be filtered a little. anyhow without the OMA genes the species effects are large.
data1 <- subset(x = data1, subset = nFeature_RNA > 350 & nCount_RNA > 500)
# data1=merge(data1,data2)
#* Including all OMA genes is equivalent, BUT still using only 1:1 to integrate; try without.
{ 
  data1 <- NormalizeData(data1, scale.factor = 1000)
  {
    list=  NULL
    vargenelist <- SplitObject(data1, split.by = "orig.ident")
    for (i in 1:length(vargenelist)) {
      vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 2000, verbose = FALSE)
    }
    for (i in 1:length(vargenelist)) {
      x <- vargenelist[[i]]@assays$RNA@var.features
      list=c(list,x)}
    list=unique(list)
    length(list)
    list=intersect(list,OMA$nematostella.gene)
    data1@assays$RNA@var.features = list
    
}
  # data1 <- FindVariableFeatures(data1, nfeatures = 2000)
  
  data1 <- ScaleData(data1,split.by = 'orig.ident')
  data1 <- RunPCA(data1, pcs.compute = 50)
  d = 1:50
  
  ### harmony ----
  data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident',)
  ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
  data1 <- RunUMAP(data1, dims = d,
                   reduction = 'harmony',
                   reduction.name ='umap.h2',reduction.key ='umap.h', #this is default; can change name to save different maps
                   #the following parameters control your output.
                   n.neighbors = 20L, #how many similar cells do you expect
                   spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                   min.dist = 0.2, #how close to plot each cell higher=spread
                   local.connectivity = 1)#overall how connected is the graph
  DimPlot(data1,group.by = 'species',cols=c('orange','lightblue','darkblue'),reduction='umap.h2')&DarkTheme()&NoAxes()
  data1<-SetIdent(data1,value = 'ID.species')
  DimPlot(data1,cols=clust.cp.separate,pt.size = 1,reduction='umap.h2',order =sort(levels(data1)))&NoAxes()&DimPlot(data1,cols=rev(clust.cp.separate[1:length(levels(data1))]),pt.size = 1,order =rev(sort(levels(data1))))&NoAxes()&plot_layout(ncol = 1)
  
  DimPlot(data1,cells.highlight = WhichCells(data1,idents = 'neuronal'))
  
  FeaturePlot(data1,c('SoxC','BRN3-OMA-1','SoxB.2','Myc2','myc1','ID4-OMA-2','Elav1'),split.by = 'species',order=T)& NoAxes() & NoLegend() &scale_colour_gradientn(colours =gene.cp)
    
  data1$ID.species=data1@active.ident
  
  #can plot a tree of cluster relationships:
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            reduction = 'harmony',
                            dims = 1:20,
                            features=intersect(OMA$nematostella.gene,TF_list),
                            reorder.numeric = F)
  # g.merge=as.data.frame(rownames(data1))
  ape::plot.phylo(data1@tools$BuildClusterTree, direction = 'down')
  
  # Use integrated data to calculate also clustering:
  data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                         # nn.method = 'annoy',  
                         #these parameters are same as UMAP uses. This is NOT seurat default
                         annoy.metric = 'cosine', k.param = 40)
  data1<- FindClusters(data1,resolution = 1)
  #calculate relationship between clusters; DISCUSS
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            dims = 1:10,
                            reorder.numeric = T)
  DimPlot(data1,reduction = 'umap.h2',
          cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Clusters')
  
}
data1<-SetIdent(data1,value='ID.species')
data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                          reduction = 'harmony',
                          dims = 1:50,
                          features=intersect(AaTF_list,TF_list),#(OMA$nematostella.gene,TF_list),
                          reorder.numeric = F)
# g.merge=as.data.frame(rownames(data1))
ape::plot.phylo(ape::root(data1@tools$BuildClusterTree,'unchar.center'))
ape::plot.phylo(data1@tools$BuildClusterTree)

DimPlot(data1,group.by = 'IDs',cols=clust.cp.)

# save(data1,file='OMA.merge.all.RObj')
save(data1,file='OMA.merge.tissue.RObj')

#drop dev data from nematostella:
# data1=SetIdent(data1,value = 'lifehistory')
# data1.all=data1
# data1=subset(data1,idents = levels(data1)[2:6])
# data1=droplevels(data1)
# re-run pipeline...

### correlations ----

# not sure how to think about this..

## neurons ----
{
  load(file='nem.neur.Robj')
  load(file='all.neurons.2.Robj')
  all.neurons=data1
  temp.n=nem.neur
  rownames(temp.n@assays$RNA@counts)=genes$name
  rownames(temp.n@assays$RNA@data)=genes$name
  rownames(temp.n@assays$RNA@meta.features)=genes$name # otherwise fails with vs5
  temp.n@assays$RNA@var.features=as.character(genes$name) # otherwise fails with vs5
  goi=genes$name[match(OMA$Nv2,genes$geneID)]
  temp.n<-ScaleData(temp.n,features = goi,split.by = 'orig.ident')
  temp.n$species='nematostella'
  temp.n$ID.species=paste('Nv',temp.n$IDs,sep='.')
  temp.n<-SetIdent(temp.n,value = 'ID.species')
  levels(temp.n)
  
  goi=genes.ac$name[match(OMA$Aur2,genes.ac$geneID)]
  d1=all.neurons
  d1@active.assay='RNA'
  d2=data1.subsets$gland.dig
  d2@active.assay='RNA'
  temp=merge(d1,d2)
  rownames(temp@assays$RNA@counts)=genes.ac$name
  rownames(temp@assays$RNA@data)=genes.ac$name
  rownames(temp@assays$RNA@meta.features)=genes.ac$name # otherwise fails with vs5
  temp@assays$RNA@var.features=genes.ac$name # otherwise fails with vs5
  # temp<-ScaleData(temp,features = goi,split.by = 'orig.ident')
  temp$species='aurelia'
  temp$ID.species=paste('Ac',temp$IDs,sep='.')
  temp<-SetIdent(temp,value = 'ID.species')
  levels(temp)
  
  data1=merge(temp,temp.n) #good to go... run the pipeline
  # this is mediocre at best with all genes; strong library effects.
  
  # use only the OMA:
  data1=subset(data1,features = OMA$name)
  # data1 <- subset(x = data1, subset = nFeature_RNA > 350 & nCount_RNA > 500)

      unique(data1$ID.species)
      
      data1$orig.ident=as.factor(data1$orig.ident)
      
      data1$orig.ident=factor(data1$orig.ident,levels(data1$orig.ident)[c(8:12,1:7,14:25,13,26:33)])
      LibCP=c(rev(stepped(21)),stepped2(10))
  
}
## Seurat ----
{
  VlnPlot(data1, features = c('nFeature_RNA', 'nCount_RNA'))
  data1 <- NormalizeData(data1, scale.factor = 1000)
  {
    list=  NULL
    vargenelist <- SplitObject(data1, split.by = "orig.ident")
    for (i in 1:length(vargenelist)) {
      vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 2000, verbose = FALSE)
    }
    for (i in 1:length(vargenelist)) {
      x <- vargenelist[[i]]@assays$RNA@var.features
      list=c(list,x)}
    list=unique(list)
    length(list)
    list=intersect(list,OMA$name)
    data1@assays$RNA@var.features = list
    
  }
  data1 <- FindVariableFeatures(data1, nfeatures = 3000)
  
  data1 <- ScaleData(data1,split.by = 'orig.ident')
  data1 <- RunPCA(data1, pcs.compute = 50)
  d = 1:50
  
  ### harmony ----
  # insufficient when using all genes
  data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')

  data1 <- RunUMAP(data1, dims = d,
                   reduction = 'harmony',
                   reduction.name ='umap.h',reduction.key ='umap.h', #this is default; can change name to save different maps
                   #the following parameters control your output.
                   n.neighbors = 10L, #how many similar cells do you expect
                   spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                   min.dist = 0.2, #how close to plot each cell higher=spread
                   local.connectivity = 1)#overall how connected is the graph
  DimPlot(data1, group.by = 'orig.ident',reduction = 'umap.h',
          cols = LibCP,pt.size = 2)+NoAxes()+labs(title='Harmony',subtitle = 'UMAP.alt | ID')#&
  DimPlot(data1,reduction = 'umap.h',
          cols = c(clust.cp.graded,clust.cp.separate))+NoAxes()+labs(title='Harmony',subtitle = 'Original IDs')
  DimPlot(data1, group.by = 'species',reduction = 'umap.h',
          cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Original IDs')
  
  #can plot a tree of cluster relationships:
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            # reduction = 'harmony',
                            # dims = d,
                            # features=intersect(OMA$nematostella.gene,TF_list),
                            # graph = 'RNA_nn',
                            reorder.numeric = F)
  
  ape::plot.phylo(data1@tools$BuildClusterTree, direction = 'downwards')
  { 
    embeddings <- Embeddings(object = data1, reduction = 'pca')[,1:50]
    data.dims=matrix('0',20L,length(levels(data1)))
    for (i in 1:length(levels(data1)))  
    {  cells <- WhichCells(object = data1, idents = levels(data1)[i])
    data.dims[,i] <- colMeans(embeddings[cells,])
    }
    colnames(x = data.dims) <- levels(x = data1)
    data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
    nj.tree <- ape::nj(data.dist)
    ape::plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = F,use.edge.length = F,lab4ut = 'axial', cex = 0.8,no.margin = T)
    
    }
   
  # Use integrated data to calculate also clustering:
  data1 <- FindNeighbors(object = data1,reduction ="harmony",dims = d,
                         # nn.method = 'annoy',  
                         #these parameters are same as UMAP uses. This is NOT seurat default
                         annoy.metric = 'cosine', k.param = 40)
  data1<- FindClusters(data1,resolution = 1)
  #calculate relationship between clusters; DISCUSS
  data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                            dims = 1:10,
                            reorder.numeric = T)
  DimPlot(data1,reduction = 'umap.h',
          cols = clust.cp.separate)+NoAxes()+labs(title='Harmony',subtitle = 'Clusters')
}

data1 <- BuildClusterTree(object = data1, reorder = TRUE,
                          reduction = 'harmony',
                          dims = 1:50,
                          features=intersect(AaTF_list,TF_list),#(OMA$nematostella.gene,TF_list),
                          reorder.numeric = F)
ape::plot.phylo(data1@tools$BuildClusterTree)

save(data1,file='OMAmergeFullneural.Robj')

# new cell type tree ----
load(file='OMAmergeFullneural.Robj')

#drop the early states first:
levels(data1)
all.neural=data1
data1<-subset(data1,idents=levels(data1)[c(1:9,13:33,35:40,44:55,57:64,67:74)])
              # [c(3:10,12:16,18:20,22:23,25:45,47:52,54:55,57:59,61:75)])
data1<-droplevels(data1)

#use only OMA TFs...
OMA.tf=OMA$name[match(TF_list,OMA$nematostella.gene,nomatch=0)]
all.tf=rownames(data1)[match(unique(c(OMA.tf,TF_list,AaTF_list)),rownames(data1))]
data1 <- ScaleData(data1,split.by = 'orig.ident',features =all.tf )
data1 <- RunPCA(data1, pcs.compute = 50,features = all.tf)

### harmony ----
# insufficient when using all genes
data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')

library(stats)
library(ape)

dim(data1)
#use PCs instead to build the tree:
embeddings <- Embeddings(object = data1, reduction = 'harmony')[,1:50]
data.dims=matrix(nrow = 50L,ncol =  length(levels(data1)))
for (i in 1:length(levels(data1)))  
{  cells <- WhichCells(object = data1, idents = levels(data1)[i])
data.dims[,i] <- colMeans(embeddings[cells,])
}

colnames(x = data.dims) <- levels(x = data1)
data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
nj.tree <- nj(data.dist)

nj.boot.tree <- boot.phylo(nj.tree, t(data.dims), FUN = function(xx) nj(dist(xx)), B = 1000)

# add bootstraps:
nj.tree$node.label <- nj.boot.tree/10
nj.tree$node.label <- round(nj.tree$node.label)

# plot tree
plot.phylo(nj.tree, type = "unrooted", label.offset = 0.2,show.node.label = F,use.edge.length = F,lab4ut = 'axial', cex = 1,no.margin = T, rotate.tree = 90)

plot.phylo(nj.tree, type = "tidy", label.offset = 0.2,show.node.label = T,use.edge.length = F,lab4ut = 'axial', cex = 1,no.margin = T, rotate.tree = 90)

all.markers_TF <- FindAllMarkers(
  data1,
  logfc.threshold = 0.4,
  features = OMA.tf, #include only TFs that are present
  min.pct = 0.05,
  only.pos = TRUE,
  return.thresh = 0.001,
)
View(all.markers_TF)

-------------------------------------------------
### seurat ----
#very time and memory intensive; does not necessarily perform better than harmony, which is very quick.
{
  # split the dataset into a list of separate seurat objects
  data1<- SetIdent(data1,value='orig.ident')
  #this was failing without SeuratObject vs5
  data1@assays$integrated=NULL
  data1@assays$alra=NULL
  data1.list <- SplitObject(data1, split.by = "orig.ident")
  
  #need to drop the mesentery libraries here; there are no cnidocytes in these so it throws an error:
  data1.list=data1.list[c(1:25,27:29)]
  # normalize and identify variable features for each dataset independently
  data1.list <- lapply(X = data1.list, FUN = function(x) {
    x <- NormalizeData(x,scale.factor = 5000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = data1.list,nfeatures = 5000)  
  data1.anchors <- FindIntegrationAnchors(object.list = data1.list, anchor.features = features)
  
  # this command creates an 'integrated' data assay
  data1.combined <- IntegrateData(anchorset = data1.anchors,k.weight = 20)
  #re-run analysis on new assay:

  data1.combined <- ScaleData(data1.combined, verbose = FALSE)
  data1.combined <- RunPCA(data1.combined, npcs = 50, verbose = FALSE)
  
  ### UMAP ----
  data1.combined <- RunUMAP(
    data1.combined,
    dims = d,
    reduction = 'pca',
    reduction.name = 'umap',
    reduction.key = 'umap',
    #this is default; can change name to save different maps
    #the following parameters control your output.
    seed.use = 42,
    #random starting point; without your maps will differ every time!
    n.neighbors = 10L,
    #how many similar cells do you expect
    spread = 1,
    #how big is the axis space: 0.2 for lineages; 1 for separation
    min.dist = 0.2,
    #how close to plot each cell higher=spread
    local.connectivity = 1
  )#overall how connected is the graph
  
  #check that output is not based only on information content:
  FeaturePlot(
    data1.combined,
    'nCount_RNA',
    reduction = 'umap',order=T,
    cols = c('lightgrey', 'darkred')
  ) + NoAxes() + NoLegend() +
    ggplot2::labs(title = 'nFeatures', subtitle = 'umap') 
  DimPlot(data1.combined,cols=clust.cp.separate,group.by = 'IDs',pt.size = 2)&NoAxes()
  
  data1.combined <- FindNeighbors(data1.combined, reduction = "pca", dims = 1:20,k.param = 45)
  data1.combined <- FindClusters(data1.combined, resolution = 0.2)
  data1.combined <- BuildClusterTree(data1.combined,reorder = T,reorder.numeric = T)
  # DimPlot(data1.combined,cols=clust.cp.separate,group.by = 'orig.ident')&NoAxes()&
  DimPlot(data1.combined,cols=clust.cp.separate)&NoAxes()
  
  DimPlot(data1.combined, group.by = 'orig.ident',reduction = 'umap',
          cols = LibCP)+NoAxes()+ggplot2::labs(title='SeuratInt',subtitle = 'UMAP.alt | ID')&
  DimPlot(data1.combined, group.by = 'IDs',reduction = 'umap',
          cols = clust.cp.separate)+NoAxes()+ggplot2::labs(title='SeuratInt',subtitle = 'Original IDs')
  
  #can plot a tree of cluster relationships:
  data1.combined<-SetIdent(data1.combined,value = 'IDs')
  data1.combined <- BuildClusterTree(object = data1.combined, reorder = TRUE,
                            dims = 1:50,
                            reorder.numeric = F)

  ape::plot.phylo(data1.combined@tools$BuildClusterTree, direction = 'right')+ape::plot.phylo(data1@tools$BuildClusterTree, direction = 'right')
}

#EED Helsinki ----
aurelia<-SetIdent(Aa.Alldata,value='ID.cluster')
aurelia@active.assay='RNA'
levels(aurelia@active.ident)=c("neuroglandular","gland.mucus",'striated muscle',"unchar.center", "gastrodermis",  "neuronal","unchar.center",  "neuroglandular","cnidocyte","epithelia.outer","cnidocyte.mat", "neuronal.N2")
aurelia@active.ident=factor(aurelia@active.ident,levels(aurelia)[c(4,10,6,2,1,7,9,8,3,5)])


DotPlot(aurelia,'RNA',features=c('cpath','otx2','otx1a','otx1b','otx1c','otx1d'),cols=c('grey40','skyblue'),scale.by = 'size',col.min = 0,idents = levels(aurelia)[c(1:2,4:10)])&#&DarkTheme()&RotatedAxis()&theme(legend.position = 'bottom')

DotPlot(Alldata.Nv2,features=c('Nem64','OtxA','OtxB','OtxC'),cols=c('grey40','orange'),scale.by = 'size',col.min = 0,idents = levels(Alldata.Nv2)[c(1,3,4,6:8,10,12:13)])&RotatedAxis()& theme(legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6,color = 'white'))&DarkTheme()

DotPlot(aurelia,'RNA',features=c('cpath','otx2','otx1a','otx1b','otx1c','otx1d'),cols=c('grey40','skyblue'),scale.by = 'size',col.min = 0,group.by = 'ID.separate')&RotatedAxis()& theme(legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6,color = 'white'))&DarkTheme()&coord_flip()

DotPlot(Alldata.Nv2,'RNA',features=c('Nem64','OtxA','OtxB','OtxC'),cols=c('grey40','skyblue'),scale.by = 'size',col.min = 0,group.by = 'ID.separate')&RotatedAxis()& theme(legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6,color = 'white'))&DarkTheme()&coord_flip()

FeaturePlot(Alldata.Nv2,c('OtxA','OtxB','OtxC'),order=T,min.cutoff = 1)&DarkTheme()&NoLegend()& NoAxes() &scale_colour_gradientn(colours =gene.cp.dark)
FeaturePlot(Aa.Alldata,c('otx2','otx1a','otx1b','otx1c'),order=T,min.cutoff = 1)&DarkTheme()& NoAxes() &scale_colour_gradientn(colours =gene.cp.dark)&NoLegend()

Aa.Alldata@active.assay='alra'
FeaturePlot(Aa.Alldata,c('cpath'),order=T,max.cutoff = )&DarkTheme()& NoAxes() & NoLegend() &scale_colour_gradientn(colours =gene.cp.dark)

DotPlot(Aa.Alldata,'RNA',features=c('soxc','myc1.1','myc3','myc4','myc2','myc5','myc.1','myc.2'),cols=c('grey40','skyblue'),scale.by = 'size',col.min = 0,group.by = 'ID.separate')&
  DotPlot(Alldata.Nv2,'RNA',features=c('SoxC','myc1','Myc3','myc4','Myc2','myc5','myc'),cols=c('grey40','skyblue'),scale.by = 'size',col.min = 0,group.by = 'ID.separate')&
  RotatedAxis()&DarkTheme()& theme(panel.grid = element_line(size=0.2,colour = 'grey'),legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6,color = 'white'))

DotPlot(Aa.Alldata,'RNA',features=c('un130-like-2','al-like-1','tead3-like-1','soxc',	
'soxb1b','soxb1a','prd13-like-1','prd13-like-2','myc1.1','myc.1','myc.2','myc3','myc4','myc2','myc5'),cols=c('grey40','skyblue'),scale.by = 'size',col.min = 0,group.by = 'ID.separate')&RotatedAxis()& theme(legend.title = element_text(size = 8),legend.position = 'bottom',legend.text = element_text(size=6,color = 'white'))&DarkTheme()

#* Can you estimate how many times TF expansion is related to putatively novel cell types?
View(as.data.frame(AaTF_list))

#* what if you run clustering on GENES and not cells:
{#* first: make an aggregate matrix for simplicity
data.g=AggregateExpression(Aa.Alldata,group.by = 'ID.separate',return.seurat = F, scale.factor = 5000,assays = 'RNA',features=AaTF_list)
# create seurat where genes are columns
data1=CreateSeuratObject(t(data.g$RNA),)
#filter genes to have at least 10 reads
data1 <- subset(x = data1, subset = nCount_RNA > 100)
#maybe should separate by sample not just cluster...
#* to do.
dim(data1)
data1 <- FindNeighbors(
  object = data1,
  features=rownames(data1),#this parameter is what UMAP uses:
  annoy.metric = 'cosine',
  #default = 'euclidean'
  # see https://www.nature.com/articles/s41592-019-0372-4
  # k.param = 10
)
}
#This doesn't produce anything very sensible.

#*
AaTF_list.filter=AaTF_list[which(rowSums(Aa.Alldata[AaTF_list,])>=50)]
DotPlot(Aa.Alldata,features= sort(AaTF_list.filter)[1:50],group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'))&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))


eprotein.aur=genes.ac$gene_short_name[match(c('AAUR2.3994','AAUR2.7419','AAUR2.7418','AAUR2.3995','AAUR2.42838','AAUR2.7416','AAUR2.7417','AAUR2.42837'),genes.ac$geneID)]
DotPlot(Aa.Alldata,features= eprotein.aur,group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'))&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

mad.aur=genes.ac$gene_short_name[match(c('AAUR2.10536','AAUR2.10535','AAUR2.47648','AAUR2.22471','AAUR2.20688'),genes.ac$geneID)]
DotPlot(Aa.Alldata,features= mad.aur,group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'))&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

ATOH8.aur=genes.ac$gene_short_name[match(c('AAUR2.42412','AAUR2.21664','AAUR2.45801','AAUR2.45802','AAUR2.19288','AAUR2.17097','AAUR2.10736','AAUR2.15768','AAUR2.25654','AAUR2.7216'),genes.ac$geneID)]
DotPlot(Aa.Alldata,features= ATOH8.aur,group.by = 'ID.separate',col.min = 0,split.by = 'lifehistory',cols=c('blue','green','purple','red'),scale.by = 'size',scale.max = 10)&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

DIM.aur=genes.ac$gene_short_name[match(c('AAUR2.14962','AAUR2.30528','AAUR2.42531','AAUR2.17881','AAUR2.30061'),genes.ac$geneID)]

DIM.nem=genes$gene_short_name[match(c('NV2.6620','NV2.12812','NV2.6622','NV2.22218'),genes$geneID)]

DotPlot(Aa.Alldata,features= DIM.aur,group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'),scale.by = 'size',scale.max = 10)&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))&
  
  DotPlot(Alldata.Nv2,features= DIM.nem,group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'))&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

HOXE.aur=genes.ac$gene_short_name[match(c('AAUR2.18593','AAUR2.31977','AAUR2.12126','AAUR2.19806','AAUR2.14283','AAUR2.11997','AAUR2.41983','AAUR2.20449'),genes.ac$geneID)]
DotPlot(Aa.Alldata,features= HOXE.aur,group.by = 'ID.separate',col.min = 0,
        split.by = 'lifehistory',cols=c('blue','green','purple','red'),
        scale.by = 'size',scale.max = 10)&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

Aa.Alldata<-SetIdent(Aa.Alldata,value='ID.separate')
DEG.TF.aurelia<-FindAllMarkers(Aa.Alldata,features=AaTF_list.filter,logfc.threshold = 1,min.pct = 0.2,only.pos = T,return.thresh = 0.01)
DotPlot(Aa.Alldata,features= unique(DEG.TF.aurelia$gene),group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'))&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

prdm6.nem=genes$gene_short_name[match(c('NV2.919','Nv2.2247','NV2.22966','NV2.22967','NV2.22968','NV2.23254','NV2.23255','NV2.23257','NV2.23258'),genes$geneID)]

DotPlot(Alldata.Nv2,features= prdm6.nem,group.by = 'ID.separate',col.min = 0,cols=c('lightgrey','red'))&RotatedAxis()&theme(panel.grid = element_line(size=0.2,colour = 'grey90'))

# OMA.all ----
# use all OMA genes: combine expression of paralogs
# first: test on smaller subset: cnidocytes:


##rLiger integration FAIL ----
#fails. bugs in installation.
data1[["RNA"]] <- split(data1[["RNA"]], f = data1$orig.ident)

data1 <- data1 %>%
  normalize() %>%
  selectGenes() %>%
  scaleNotCenter()

data1 <- data1 %>%
  runINMF(k = 20) %>%
  quantileNorm()


#harmony v5 ----
#update to 'layered' v5 object
obj=data1
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
obj@reductions$pca=NULL

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,split.by = 'orig.ident')
obj <- RunPCA(obj)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE)

obj <- RunUMAP(obj, dims = d,
                 reduction = 'harmony',
                 reduction.name ='umap.h2',reduction.key ='umap.h', #this is default; can change name to save different maps
                 #the following parameters control your output.
                 n.neighbors = 10L, #how many similar cells do you expect
                 spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                 min.dist = 0.2, #how close to plot each cell higher=spread
                 local.connectivity = 10)#overall how connected is the graph
DimPlot(data1, group.by = 'species',reduction = 'umap.h',
        cols = LibCP)+NoAxes()+labs(title='Harmony',subtitle = 'original')&DimPlot(data1, group.by = 'species',reduction = 'umap.h2',cols = LibCP)+NoAxes()+labs(title='Harmony',subtitle = 'IntegrateLayers')

#this worked fine. removed and not saved.

cnido.nem=subset(AllData.Nv2$IDs.separate)

# fuck off python
# SAMap ----
mappings=read.csv(file = '/lisc/user/agcole/jupnotebooks/nvac.cnido.mapping.cvs',row.names = 1)
mappings=read.csv(file = '/lisc/user/agcole/jupnotebooks/nvac.cnido.mapping.cvs',row.names = 1)
pheatmap::pheatmap(mappings,cluster_rows = T)

# this is OK; doesn't add anything really to the cell type tree...

# OMA 2 with Hydra ----
load (file='/lisc/data/scratch/molevo/blickenstorfer/OrthoMAP_jmontenegro/results/OrthoMAP_Seurat_Objects/oma/result.Robj')
load (file='/lisc/data/scratch/molevo/blickenstorfer/OrthoMAP_acole/results/OrthoMAP_Seurat_Objects/oma/result.Robj')

data1<-SetIdent(orthomap_obj,value='ID.separate')
coi=WhichCells(data1,ident=levels(data1)[c(2,11,39,41,43:47,50)])
all.cnido=subset(data1,idents = levels(data1)[c(c(2,11,21,39,41,43:47,50))])
DimPlot(data1,cells.highlight =coi ) #seurat v4 fails with v5 object

data1<-SetIdent(all.cnido,value='ID.separate') #this is missing from this object...
# DimPlot(all.cnido,cols=clust.cp.separate)
# load (file='/lisc/data/scratch/molevo/blickenstorfer/OrthoMAP_acole/results/OrthoMAP_Seurat_Objects/oma/result.Robj') #these don't match
# 
# data1<-SetIdent(orthomap_obj,value='ID.separate')
# coi.h=WhichCells(SetIdent(all.cnido,value = 'species'),ident='Hv')
# coi=colnames(all.cnido)
# all.cnido$ID.separate = all.cnido$IDs
# all.cnido$ID.separate = as.character(all.cnido$ID.separate)
# all.cnido$ID.separate[coi] = as.character(data1$ID.separate)[coi]
# all.cnido$ID.separate[coi.h]=as.character(all.cnido$IDs)[coi.h]
# levels(as.factor(all.cnido$ID.separate))
# all.cnido<-SetIdent(all.cnido,value="ID.separate")
# DimPlot(all.cnido,cols=clust.cp.separate)
