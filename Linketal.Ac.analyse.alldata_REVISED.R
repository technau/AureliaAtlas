#setup ----
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
library(dplyr)

load(file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE.genes.RData') 

gene.cp = c('lightgrey', rev(brewer.pal(11 , "Spectral")))
clust.cp.separate = unique (c(cols25(25), glasbey(32), alphabet(26)))
clust.cp.graded = unique(c(stepped3(20), stepped2(20), stepped(20)))
LibCP =brewer.paired(12)
LibCP.stages=LibCP[c(2,7,10,12)]
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
do.over=F
run.seurat=F
if (do.over)
{ load(file='AcAlldata.raw.16.Robj')
  load(file='single.Alldata.16.ROBJ')
Ac.Alldata=Ac.Alldata.unprocessed}else{load(file='Ac.Alldata.Robj')}

if (do.over)
{
  # Generate Integrated Dataset ----
  # normalize and identify variable features for each dataset independently
  single.AllData <- lapply(single.AllData, function(x) {
    x <- NormalizeData(x, scale.factor = 5000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x)
    x <- RunPCA(x)
    x <- RunUMAP(x,dims=1:20)
  })

  #make a list of ALL variable features
  #* this was necessary for Seurat integration, but you don't use these again with harmony
  vargenelist= NULL
  for (i in 1:length(single.AllData)) {
    x <- single.AllData[[i]]@assays$RNA@var.features
    vargenelist=c(vargenelist,x)}
  vargenelist=unique(vargenelist)
  length(vargenelist)

  #save and clean up the workspace
  save(single.AllData,file='single.Alldata.16.Robj')
  rm(single.AllData, Ac.Alldata.unprocessed)
  gc()
  
  data1=SetIdent(Ac.Alldata,value='orig.ident')
  data1$orig.ident = as.factor(data1$orig.ident)
  lib.names.original=dput(levels((data1$orig.ident)))
  lib.names.original
  order = c(8:16,1:7)
  data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,
                                       levels(data1@meta.data$orig.ident)[order])
  data1=  SetIdent(data1,value = 'orig.ident')
  Ac.Alldata=data1
  save(Ac.Alldata,file='Ac.Alldata.Robj')
}  

{
 data1=SetIdent(Ac.Alldata,value='orig.ident')
  data1@active.assay='RNA'
  polyp.all=WhichCells(data1,idents = levels(data1)[1:7])
  length(polyp.all)
  medusa.all = WhichCells(data1,idents = levels(data1)[12:16])
  length(medusa.all)
  strobila.all=WhichCells(data1,idents = levels(data1)[8:9])
  length(strobila.all)
  ephyra.all = WhichCells(data1,idents = levels(data1)[10:11])
  length(ephyra.all)
  data1$lifehistory= data1$orig.ident
  data1$lifehistory=as.character(data1$lifehistory)
  data1$lifehistory[polyp.all]='polyp'
  data1$lifehistory[strobila.all]='strobila'
  data1$lifehistory[ephyra.all]='ephyra'
  data1$lifehistory[medusa.all]='medusa'
  
  data1$lifehistory.tissue= data1$orig.ident
  data1$lifehistory.tissue=as.character(data1$lifehistory.tissue)
  data1$lifehistory.tissue[polyp.all]='polyp'
  data1$lifehistory.tissue[strobila.all]='strobila'
  data1$lifehistory.tissue[ephyra.all]='ephyra'
  data1$lifehistory.tissue=as.factor(data1$lifehistory.tissue)
  data1$lifehistory.tissue<-factor(data1$lifehistory.tissue,levels=levels(data1$lifehistory.tissue)[c(7,8,6,1:5)])
  data1$lifehistory.tissue<-factor(data1$lifehistory.tissue,levels=levels(data1$lifehistory.tissue)[c(1,2,4,5,6:8,3)])
  
  data1=SetIdent(data1,value='lifehistory')
  Ac.Alldata$lifehistory = data1@active.ident #saved
  data1=SetIdent(data1,value='lifehistory.tissue')
  Ac.Alldata$lifehistory.tissue = data1@active.ident #saved
  
  data1=SetIdent(Ac.Alldata,value='lifehistory')
  data1@active.assay='RNA'
  
  VlnPlot(data1,c('nFeature_RNA','nCount_RNA'),cols=LibCP.stages)
  
  keep = which(rowSums(data1@assays$RNA@counts) >= 10) #17560 !! :)
  length(keep)

  coi.8000=WhichCells(SetIdent(Ac.Alldata,value='lifehistory'),
                      downsample = 8000,seed = NULL)
  downsampled.data1=subset(data1,cells = coi.8000,features = keep) #memory issue here.
  genes=as.data.frame(rownames(Ac.Alldata@assays$RNA))
  colnames(genes)='gene_short_name'
  genes=genes[keep,]
  genes=as.data.frame(genes)
  colnames(genes)='gene_short_name'
  data1=downsampled.data1
  
  #replicates
  data1.temp=subset(SetIdent(Ac.Alldata,value='lifehistory'),features = keep)
  data1.temp@active.assay='RNA'
  {
    nr=100 #number of replicates
    replicates = NULL
    coi=NULL
    for (r in 1:nr)
    {
      coi[[r]]=WhichCells(data1.temp,downsample = 8000,seed = NULL)
    }
    
    
    goi.lifehistory = NULL
    num.gene.all.rep = matrix (0L,nr,4)#4
    # data.temp=data1
    for (r in 1:nr)
    {
      data1=subset(data1.temp,cells=coi[[r]])
      data1@active.assay='RNA'
      polyp=WhichCells(data1,idents = levels(data1)[1])
      medusa = WhichCells(data1,idents = levels(data1)[4])
      strobila=WhichCells(data1,idents = levels(data1)[2])
      ephyra = WhichCells(data1,idents = levels(data1)[3])
      polyp.genes.ac=NULL
      polyp.genes.ac = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[
        ,polyp])>=3)] 
      medusa.genes.ac=NULL
      medusa.genes.ac = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,medusa])>=3)] 
      strobila.genes.ac = NULL
      strobila.genes.ac = genes.ac$gene_short_name[which(rowSums(data1@assays$RNA@counts[,strobila])>=3)]
      ephyra.genes.ac = NULL
      ephyra.genes.ac = genes.ac$gene_short_name[which(rowSums(data1@assays$RNA@counts[,ephyra])>=3)] 
      total.counts=NULL  
      mean.total.counts = NULL
      num.gene.all = matrix (0L,1,4)#4
      for (l in 1:4)
      {total.counts[[l]]=colSums(data1@assays$RNA@data[,WhichCells(data1,idents=levels(data1)[l])])
      mean.total.counts[l]=mean(total.counts[[l]])}
      a1=(polyp.genes.ac)
      a3=(strobila.genes.ac)
      a4=(ephyra.genes.ac)
      a2=(medusa.genes.ac)
      
      num.gene.all.rep[r,]=c(length(a1),length(a3),length(a4),length(a2))
      ST.venn.plot <- VennDiagram::venn.diagram(list(a1,a2,a3,a4),filename = NULL,category = c('polyp','medusa','strobila','ephyra'), fill=LibCP[c(2,8,4,5)],#margin=1,
    cex=2,print.mode = 'percent',sigdigs=2,disable.logging=F)
      
      colnames(num.gene.all.rep) =c('polyp','strobila','ephyra','medusa')
      # barplot(num.gene.all.rep)
      replicates[[r]] = ST.venn.plot
      goi.lifehistory[[r]]=VennDiagram::calculate.overlap(list(a1,a2,a3,a4))
    }
    
    # grid::grid.newpage()
    # grid::grid.draw(replicates[[1]])
    # boxplot(num.gene.all.rep,col=LibCP[c(2,4,5,8)])+title('Number of Expressed genes.ac')
    save(num.gene.all.rep,file='replicates8000.kostya.RObj')
  }   


  ### gene usage ----
  data1=SetIdent(Ac.Alldata,value='lifehistory')
  polyp=WhichCells(data1,idents = levels(data1)[1])
  strobila=WhichCells(data1,idents = levels(data1)[2])
  ephyra=WhichCells(data1,idents = levels(data1)[3])
  medusa = WhichCells(data1,idents = levels(data1)[4])
  {
    n=1
    polyp.genes.ac=NULL
    polyp.genes.ac = genes.ac$gene_short_name[which(rowSums(data1@assays$RNA@counts[,polyp])>=n)] #31103 c | 29885 d
    medusa.genes.ac=NULL
    medusa.genes.ac = genes.ac$gene_short_name[which(rowSums(data1@assays$RNA@counts[,medusa])>=n)] #40555 | 37903 d
    strobila.genes.ac = NULL
    strobila.genes.ac = genes.ac$gene_short_name[which(rowSums(data1@assays$RNA@counts[,strobila])>=n)] #34476 | 31176 d
    ephyra.genes.ac = NULL
    ephyra.genes.ac = genes.ac$gene_short_name[which(rowSums(data1@assays$RNA@counts[,ephyra])>=n)] #36893 | 33921 d
    total.counts=NULL  
    mean.total.counts = NULL
    num.gene.all = matrix (0L,1,4)#4
    for (l in 1:4)
    {total.counts[[l]]=colSums(data1@assays$RNA@data[,WhichCells(data1,idents=levels(data1)[l])])
    mean.total.counts[l]=mean(total.counts[[l]])}
    a1=(polyp.genes.ac)
    a3=(strobila.genes.ac)
    a4=(ephyra.genes.ac)
    a2=(medusa.genes.ac)
    num.gene.all=c(length(a1),length(a3),length(a4),length(a2))
    ST.venn.plot <- VennDiagram::venn.diagram(list(a1,a2,a3,a4),filename = NULL,
                                              category = c('polyp','medusa','strobila','ephyra'),
                                              fill=LibCP[c(2,8,4,5)],
                                              cex=2,disable.logging=F,sigdigs = 0,print.mode = 'percent')
    ST.venn.plot.raw <- VennDiagram::venn.diagram(list(a1,a2,a3,a4),filename = NULL,
                                                  category = c('polyp','medusa','strobila','ephyra'),
                                                  fill=LibCP[c(2,8,4,5)],
                                                  cex=2,disable.logging=F,sigdigs = 0,print.mode = 'raw')
    
    goi.lifehistory.all=VennDiagram::calculate.overlap(list(a1,a2,a3,a4))
    
  }
  
  medusa.specific=setdiff(c(medusa.genes.ac,ephyra.genes.ac),c(polyp.genes.ac,strobila.genes.ac))
  medusa.specific=genes.ac[match(medusa.specific,genes.ac$gene_short_name),] #1172
  View(medusa.specific)
  # sanity check:
 DotPlot(data1,features=medusa.specific$gene_short_name,col.min = 0)&RotatedAxis()&coord_flip()
  write.csv(medusa.specific,file='DS1.1a.csv')
  
  ### TopGO medusa ----
  run.topGO=T
  if (run.topGO)
  {
  {
    library(topGO)
    
    Unfold <- annot.red.Plus %>% 
      dplyr::mutate(`go` = strsplit(as.character(`go`), ",")) %>% 
      tidyr::unnest(`go`) 
    geneID2GO <- Unfold %>% split(x = .$`go`, f = .$gene_short_name.x)
    geneNames <- names(geneID2GO)
    
    myInterestingGenes = medusa.specific$gene_short_name
    geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
    names(geneList) <- geneNames
    
    GOdata.BP <- new("topGOdata",ontology = "BP", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
    GOdata.MF <- new("topGOdata",ontology = "MF", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
  }
  {
    p=NULL
    resultFis=NULL
    topGO.medusa=NULL
    go.type=c("Biological Process",'Molecular Function')
    for (i in 1:2)
    {    if(i==1)
    {GOdata=GOdata.BP}else if(i==2)
    {GOdata=GOdata.MF}  
      #run the test...
      resultFis[[i]] <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      pvalFis <- score(resultFis[[i]])
      #filter for only >0.05
      pvalFis = pvalFis[pvalFis>=0.05]
      allRes_intgenes<- GenTable(GOdata, pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=30)
      allRes_intgenes$pvalues<-as.numeric(allRes_intgenes$pvalues)
      #convert NA values to zero; only if p value is so small it cannot be displayed in r
      allRes_intgenes[is.na(allRes_intgenes)]<-0.00000001
      
      #plot GOenrichment 
      topGO.medusa[[i]]=
        ggplot2::ggplot(allRes_intgenes, ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
        ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
                              col='#FF7F00',fill='#FF7F00') +
        ggplot2::coord_flip()+
        ggplot2::xlab(go.type[[i]]) +
        ggplot2::ylab("Enrichment -log10 p-value") +
        ggplot2::labs(title="GO enrichment",subtitle='medusa specific genes.ac')
      
      save(topGO.medusa,file='topGO.medusa')
    }
  }
  }
  ### DEGs ----
  data1=SetIdent(Ac.Alldata,value='lifehistory')
  data1@active.assay='RNA'
  # all markers
  {
    life.history.markers=FindAllMarkers(data1,logfc.threshold = 0.6, return.thresh = 0.001,  min.pct = 0.3, only.pos = T, verbose=F)
    
    life.history.markers = update_marker_list(life.history.markers)
    save(life.history.markers,file='markers.lifehistory.RObj')
    write.csv(life.history.markers,file='DS1.2.csv')
  }
  # TFs only
  {

    life.history.markers.TF=FindAllMarkers(data1,
                                           return.thresh = 0.01,min.pct = 0.05,only.pos = T,features = intersect(AcTF_list,rownames(data1)),verbose=F)
    
    life.history.markers.TF = update_marker_list(life.history.markers.TF)
    save(life.history.markers.TF,file='TF.markers.lifehistory.RObj')
    write.csv(life.history.markers.TF,file='DS1.2b.csv')
  }
  ### TopGO Life history ----
  {
    library(topGO )
    
    Unfold <- annot.red.Plus %>% 
      dplyr::mutate(`go` = strsplit(as.character(`go`), ",")) %>% 
      tidyr::unnest(`go`) 
    geneID2GO <- Unfold %>% split(x = .$`go`, f = .$gene_short_name.x)
    geneNames <- names(geneID2GO)
    #set dataset:
    all.markers =life.history.markers
    #process data
    GOdata = NULL
    myInterestingGenes = NULL
    geneList = NULL
    for(i in 1:length(unique(all.markers$cluster)))
    {
      #identify the GOI from the DE list
      myInterestingGenes[[i]] <- all.markers$gene[all.markers$cluster==unique(all.markers$cluster)[i]] #list of genes.ac you want to perform GO enrichment for
      geneList[[i]] <- factor(as.integer(geneNames %in% myInterestingGenes[[i]]))
      names(geneList[[i]]) <- geneNames
      GOdata[[i]] <- new("topGOdata",ontology = "BP", allGenes = geneList[[i]],annot = annFUN.gene2GO, gene2GO = geneID2GO)
    }
    # filter and generate figures.
    p=NULL
    resultFis=NULL
    topGO.lifehistory=NULL
    for(i in 1:length(unique(all.markers$cluster)))
    {
      #run the test...
      resultFis[[i]] <- runTest(GOdata[[i]], algorithm = "classic", statistic = "fisher") 
      pvalFis <- score(resultFis[[i]])
      #filter for only >0.05
      pvalFis = pvalFis[pvalFis>=0.5]
      allRes_intgenes<- GenTable(GOdata[[i]], pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=30)
      allRes_intgenes$pvalues<-as.numeric(allRes_intgenes$pvalues)
      go_ids <- allRes_intgenes$GO.ID
      # library(simplifyEnrichment)
      
      # go_scores <- allRes_intgenes$pvalues
      # names(go_scores) <- allRes_intgenes$GO.ID
      # simplifyGO(
      #   go_scores,
      #   ont = "BP",
      #   # organism = "human",   # change if needed (e.g., "mouse")
      #   method = "binary_cut"
      # )
      #convert NA values to zero; only if p value is so small it cannot be displayed in r
      allRes_intgenes[is.na(allRes_intgenes)]<-0.00000001
      
      #plot GOenrichment 
      
      david_col=LibCP.stages
      topGO.lifehistory[[i]]=
        ggplot2::ggplot(allRes_intgenes, ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
        ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
                              col=david_col[i],fill=david_col[i]) +
        ggplot2::coord_flip()+
        ggplot2::xlab("Biological Process") +
        ggplot2::ylab("Enrichment -log10 p-value") +
        ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
      
      # print(topGO.lifehistory[[i]])
    }
    
    save(topGO.lifehistory,file='lifehistory.topGO.RObj')
  } 
}    
if (run.seurat)
  # Run seurat ----
{
  ##harmony integration ----
  {    
    data1<-Ac.Alldata
    data1@active.assay='RNA'
    ### 3.1a normalize ----
    data1 <- NormalizeData(data1, scale.factor =5000)
    ### 3.1b calculate variable genes.ac ----
    data1 <- FindVariableFeatures(data1, nfeatures = 2000)
    VariableFeaturePlot(data1)
    
    list=  NULL
    vargenelist <- SplitObject(data1, split.by = "orig.ident")
    for (i in 1:length(vargenelist)) {
      vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 1000, verbose = FALSE)
      x <- vargenelist[[i]]@assays$RNA@var.features
      list=c(list,x)
    }
      list=unique(list)
      length(list)
      data1@assays$RNA@var.features = list
    data1 <- ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50)
    data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident',)
    ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
    pct <- data1[["harmony"]]@stdev / sum(data1[["harmony"]]@stdev) * 100
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    d=1:co1
    d
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    d = 1:co2
    
    data1 <- RunUMAP(data1, dims = d,
                     reduction = 'harmony',
                     reduction.name ='umap.h',reduction.key ='umap.h',
                     n.neighbors = 30L, #how many similar cells do you expect
                     spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                     min.dist = 0.4, #how close to plot each cell higher=spread
                     local.connectivity = 10,#overall how connected is the graph
                     seed.use = 42)
    
    data1 <- FindNeighbors(data1, 
                           reduction = "harmony", 
                           dims = d,
                           annoy.metric = 'cosine',
                           verbose=F,
                           k.param = 30) 
    
    data1 <- FindClusters(data1, verbose=F,
                          resolution = 0.1) 
    data1 <- BuildClusterTree(data1,
                              reorder = TRUE,
                              reorder.numeric = T,verbose = F)
    data1$h.cluster=data1@active.ident
      }
  Ac.Alldata=data1
  data1@active.assay = 'RNA'
  data1<-SetIdent(Ac.Alldata,value = 'h.cluster')
}
#name clusters ID ----
{
  clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                            sheet = 'AGC.markers.integrated')
  

  goi = genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)] #account for the updated names; this also doesn't work!!!
  goi
  
  DotPlot(data1,'RNA',goi,col.min = 0)&RotatedAxis()&DimPlot(data1,cols=clust.cp.separate,label=T)+NoAxes()+NoLegend()+coord_flip()
  #assign cluster ID to the individual libraries
  data1@active.assay='RNA'
  data1<-ScaleData(data1,features = goi,verbose=F)#
  clName = vector()
  
  #generate a matrix of values of each cluster for each gene:
  m=AverageExpression(data1,features=goi,slot='scale.data') #Seurat has a function for this
  
  for (j in 1:length(levels(data1@active.ident))) #for each cluster set
  {
    clName[j]=as.integer(which.max(m$RNA[,j])) #choose the highest value
  }
  sort(clName) #this prints which IDs were selected; sanity check that it did what you wanted... should have no unwanted replicated IDs
  clusterNames$ID[clName]
  #Note: here we merge the cnidocytes, dig. gland, and neural
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay = 'RNA'
  DimPlot(data1,cols=clust.cp.separate,label=T)&NoLegend()
  Ac.Alldata=data1
  save(Ac.Alldata,file='Ac.Alldata.Robj')
}

# DEGs ----
{  
  data1=SetIdent(Ac.Alldata,value='IDs')
  {
    data1@active.assay='RNA'
    all.markers <- FindAllMarkers(data1,
                                  return.thresh = 0.001,
                                  min.pct = 0.2,
                                  only.pos = TRUE,
                                  verbose = F)
    
    all.markers<-update_marker_list(all.markers)
    
  }
  all.markers.Alldata=all.markers 
  #TFs
  {all.markers <- FindAllMarkers(data1,
                                 return.thresh = 0.001,
                                 features = intersect(AcTF_list,rownames(data1)),
                                 min.pct = 0.05,
                                 only.pos = TRUE,
                                 verbose = F)
    
    all.markers<-update_marker_list(all.markers)
    }  
  all.markers.Alldata.TF=all.markers
  save(all.markers.Alldata,all.markers.Alldata.TF,file = 'AllData.ID.markers.RData')
  write.csv(all.markers.Alldata,file='DS1.3.csv')
  write.csv(all.markers.Alldata.TF,file='DS1.3b.csv')
}

# topGO ----
Unfold <- annot.red.Plus %>% 
  dplyr::mutate(`go` = strsplit(as.character(`go`), ",")) %>% 
  tidyr::unnest(`go`) 
geneID2GO <- Unfold %>% split(x = .$`go`, f = .$gene_short_name.x)
geneNames <- names(geneID2GO)
#set dataset:
all.markers=all.markers.Alldata
#process data
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
  
  david_col=clust.cp.separate
  for(i in 1:length(unique(all.markers$cluster)))  
  {
  topGO.IDs[[i]]=
    ggplot2::ggplot(allRes_intgenes[[i]], ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
    ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
                          col=david_col[i],fill=david_col[i]) +
    ggplot2::coord_flip()+
    ggplot2::xlab("Biological Process") +
    ggplot2::ylab("Enrichment -log10 p-value") +
    ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
}  
  for(i in 1:length(topGO.IDs))  
   print(topGO.IDs[[i]])


save(topGO.IDs,allRes_intgenes,file='topGO.IDs.RData')

