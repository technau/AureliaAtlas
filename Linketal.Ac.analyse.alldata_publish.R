

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
do.over=T
run.seurat=T
if (do.over)
{load(file='Ac_manuscript_final/AaAlldata.raw.new.Robj')
Aa.Alldata=Aa.Alldata.unprocessed}else{load(file='Ac_manuscript_final/AaAlldata.Robj')}

if (do.over)
{
  # normalize and identify variable features for each dataset independently
  single.AllData=SplitObject(Aa.Alldata,split.by = 'orig.ident')
  single.AllData <- lapply(single.AllData, function(x) {
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
  
  single.AllData <- lapply(single.AllData, function(x) {
    x <- ScaleData(x,features=features.red)
    x <- RunPCA(x)
  })
  save(single.AllData,file='Ac_manuscript_final/single.Alldata.new.ROBJ')
  #* Perform integration
  anchors.rpca <- FindIntegrationAnchors(object.list = single.AllData,
                                         scale=F,
                                         normalization.method = 'LogNormalize',
                                         reduction='rpca',
                                         anchor.features = features.red)
  save(anchors.rpca,file='Ac_manuscript_final/anchors.single.AAUR2.strobila_ephyra2.RObj')
  
  #if run in parts, this would be needed to reload anchors and integrate
  # load(file='Ac_manuscript_final/anchors.single.AAUR2.strobila_ephyra2.RObj')
  
  vargenelist= NULL
  for (i in 1:12) {
    x <- anchors.rpca@object.list[[i]]@assays$RNA@var.features
    vargenelist=c(vargenelist,x)}
  vargenelist=unique(vargenelist)
  length(vargenelist)

  goi.integrate = unique(c(vargenelist,#all variable genes from the separate datasets
                           AaTF_list,cell.cycle.goi))#,#all expressed TFs
  
  Aa.Alldata <- IntegrateData(anchorset = anchors.rpca, 
                              features.to.integrate = goi.integrate)
  DefaultAssay(Aa.Alldata) <- "integrated"
  Aa.Alldata@assays$integrated@var.features=(Aa.Alldata@assays$integrated@var.features)
  
  data1=SetIdent(Aa.Alldata,value='orig.ident')
  data1$orig.ident = as.factor(data1$orig.ident)
  lib.names.original=dput(levels((data1$orig.ident)))
  order = c(8:12,1:7)
  data1@meta.data$orig.ident = factor(data1@meta.data$orig.ident,
                                      levels(data1@meta.data$orig.ident)[order])
  data1=  SetIdent(data1,value = 'orig.ident')
  Aa.Alldata=data1
  save(Aa.Alldata,file='Ac_manuscript_final/AaAlldata.Robj')
}  

{
  
  data1=SetIdent(Aa.Alldata,value='orig.ident')
  data1@active.assay='RNA'
  polyp.all=WhichCells(data1,idents = levels(data1)[1:3])
  medusa.all = WhichCells(data1,idents = levels(data1)[8:12])
  strobila.all=WhichCells(data1,idents = levels(data1)[4:5])
  ephyra.all = WhichCells(data1,idents = levels(data1)[6:7])
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
  
  data1=SetIdent(data1,value='lifehistory')
  Aa.Alldata$lifehistory = data1@active.ident #saved
  data1=SetIdent(data1,value='lifehistory.tissue')
  Aa.Alldata$lifehistory.tissue = data1@active.ident #saved
  
  data1=SetIdent(Aa.Alldata,value='lifehistory')
  data1@active.assay='RNA'
  
  # VlnPlot(data1,c('nFeature_RNA','nCount_RNA'),cols=LibCP.stages)
  
  keep = which(rowSums(data1@assays$RNA@counts) >= 10)
  coi.8000=WhichCells(SetIdent(Aa.Alldata,value='lifehistory'),
                      downsample = 8000,seed = NULL)
  downsampled.data1=subset(data1,cells = coi.8000,features = keep) #memory issue here.
  genes=Aa_genes[keep,]
  data1=downsampled.data1
  
  #replicates
  data1.temp=subset(SetIdent(Aa.Alldata,value='lifehistory'),features = keep)
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
      polyp.genes=NULL
      polyp.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[
        ,polyp])>=3)] 
      medusa.genes=NULL
      medusa.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,medusa])>=3)] 
      strobila.genes = NULL
      strobila.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,strobila])>=3)]
      ephyra.genes = NULL
      ephyra.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,ephyra])>=3)] 
      total.counts=NULL  
      mean.total.counts = NULL
      num.gene.all = matrix (0L,1,4)#4
      for (l in 1:4)
      {total.counts[[l]]=colSums(data1@assays$RNA@data[,WhichCells(data1,idents=levels(data1)[l])])
      mean.total.counts[l]=mean(total.counts[[l]])}
      a1=(polyp.genes)
      a3=(strobila.genes)
      a4=(ephyra.genes)
      a2=(medusa.genes)
      
      num.gene.all.rep[r,]=c(length(a1),length(a3),length(a4),length(a2))
      ST.venn.plot <- VennDiagram::venn.diagram(list(a1,a2,a3,a4),filename = NULL,
                                                category = c('polyp','medusa','strobila','ephyra'),
                                                fill=LibCP[c(2,8,4,5)],#margin=1,
                                                cex=2,print.mode = 'percent',sigdigs=2,disable.logging=F)
      colnames(num.gene.all.rep) =c('polyp','strobila','ephyra','medusa')
      # barplot(num.gene.all.rep)
      replicates[[r]] = ST.venn.plot
      goi.lifehistory[[r]]=VennDiagram::calculate.overlap(list(a1,a2,a3,a4))
    }
    
    # grid::grid.newpage()
    # grid::grid.draw(replicates[[1]])
    # boxplot(num.gene.all.rep,col=LibCP[c(2,4,5,8)])#+title('Number of Expressed Genes')
    save(num.gene.all.rep,file='Ac_manuscript_final/replicates8000.StrEph.RObj')
  }   
  
  ### gene usage ----
  data1=downsampled.data1
  polyp=WhichCells(data1,idents = levels(data1)[1])
  strobila=WhichCells(data1,idents = levels(data1)[2])
  ephyra=WhichCells(data1,idents = levels(data1)[3])
  medusa = WhichCells(data1,idents = levels(data1)[4])
  {
    n=3
    polyp.genes=NULL
    polyp.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,polyp])>=n)] #31103 c | 29885 d
    medusa.genes=NULL
    medusa.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,medusa])>=n)] #40555 | 37903 d
    strobila.genes = NULL
    strobila.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,strobila])>=n)] #34476 | 31176 d
    ephyra.genes = NULL
    ephyra.genes = genes$gene_short_name[which(rowSums(data1@assays$RNA@counts[,ephyra])>=n)] #36893 | 33921 d
    total.counts=NULL  
    mean.total.counts = NULL
    num.gene.all = matrix (0L,1,4)#4
    for (l in 1:4)
    {total.counts[[l]]=colSums(data1@assays$RNA@data[,WhichCells(data1,idents=levels(data1)[l])])
    mean.total.counts[l]=mean(total.counts[[l]])}
    a1=(polyp.genes)
    a3=(strobila.genes)
    a4=(ephyra.genes)
    a2=(medusa.genes)
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
  
  medusa.specific=setdiff(c(medusa.genes,ephyra.genes),c(polyp.genes,strobila.genes))
  medusa.specific=Aa_genes[match(medusa.specific,Aa_genes$gene_short_name),c(3:8)]
  

  ### TopGO medusa ----
  run.topGO=T
  if (run.topGO)
  {
  {
    library(topGO)
    
    Unfold <- Aa_annotations %>% 
      dplyr::mutate(`GOs` = strsplit(as.character(`GOs`), ",")) %>% 
      tidyr::unnest(`GOs`) 
    geneID2GO <- Unfold %>% split(x = .$`GOs`, f = .$name)
    geneNames <- names(geneID2GO)
    
    myInterestingGenes = medusa.specific$name
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
        ggplot2::labs(title="GO enrichment",subtitle='medusa specific genes')
      
      save(topGO.medusa,file='Ac_manuscript_final/topGO.medusa')
    }
  }
  }
  ### DEGs ----
  {
    data1=SetIdent(Aa.Alldata,value='lifehistory')
    data1@active.assay='RNA'
    life.history.markers=FindAllMarkers(data1,logfc.threshold = 0.6, return.thresh = 0.0001,  min.pct = 0.3,       only.pos = T, verbose=F)
    
    life.history.markers[,8:12] = 'NA'
    names(life.history.markers)[8:12]=names(Aa_genes)[c(3,4,6,7,8)]
    # add GO terms associated with this list:
    for (i in 1:length(levels(data1@active.ident))) # 
    {
      x=life.history.markers[as.numeric(life.history.markers$cluster)==i,][1:length(which(as.numeric(life.history.markers$cluster)==i)),7]
      x2=Aa_genes[match(x,Aa_genes$gene_short_name),c(3,4,6,7,8)]
      life.history.markers[as.numeric(life.history.markers$cluster)==i,][1:length(which(as.numeric(life.history.markers$cluster)==i)),8:12]<-x2
    }
    save(life.history.markers,file='Ac_manuscript_final/markers.lifehistory.RObj')
    write.csv(life.history.markers,file='Ac_manuscript_final/DS1.2.csv')
  }
  {
    data1=SetIdent(Aa.Alldata,value='lifehistory')
    data1@active.assay='RNA'
    life.history.markers.TF=FindAllMarkers(data1,
                                           return.thresh = 0.0001,min.pct = 0.05,only.pos = T,features = intersect(AaTF_list,rownames(data1)),verbose=F)
    
    life.history.markers.TF[,8:12] = 'NA'
    names(life.history.markers.TF)[8:12]=names(Aa_genes)[c(1,3:6)]
    # add GO terms associated with this list:
    for (i in 1:length(levels(data1@active.ident))) # 
    {
      x=life.history.markers.TF[as.numeric(life.history.markers.TF$cluster)==i,][1:length(which(as.numeric(life.history.markers.TF$cluster)==i)),7]
      x2=Aa_genes[match(x,Aa_genes$gene_short_name),c(1,3:6)]
      life.history.markers.TF[as.numeric(life.history.markers.TF$cluster)==i,][1:length(which(as.numeric(life.history.markers.TF$cluster)==i)),8:12]<-x2
    }
    save(life.history.markers.TF,file='Ac_manuscript_final/TF.markers.lifehistory.RObj')
    write.csv(life.history.markers.TF,file='Ac_manuscript_final/DS1.2b.csv')
  }
  ### TopGO Life history ----
  {
    library(topGO)
    
    Unfold <- Aa_annotations %>% 
      dplyr::mutate(`GOs` = strsplit(as.character(`GOs`), ",")) %>% 
      tidyr::unnest(`GOs`) 
    geneID2GO <- Unfold %>% split(x = .$`GOs`, f = .$gene_short_name)
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
      myInterestingGenes[[i]] <- all.markers$gene[all.markers$cluster==unique(all.markers$cluster)[i]] #list of genes you want to perform GO enrichment for
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
      pvalFis = pvalFis[pvalFis>=0.05]
      allRes_intgenes<- GenTable(GOdata[[i]], pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=30)
      allRes_intgenes$pvalues<-as.numeric(allRes_intgenes$pvalues)
      #convert NA values to zero; only if p value is so small it cannot be displayed in r
      allRes_intgenes[is.na(allRes_intgenes)]<-0.00000001
      
      #plot GOenrichment 
      
      david_col=LibCP[c(2,10,4,8)]
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
    
    save(topGO.lifehistory,file='Ac_manuscript_final/lifehistory.topGO.StrEph.RObj')
  }
}    
if (run.seurat)
  # Run seurat ----
{
  ##harmony integration ----
  {    
    data1<-Aa.Alldata
    data1@active.assay='RNA'
    ### 3.1a normalize ----
    data1 <- NormalizeData(data1, scale.factor = 5000)
    ### 3.1b calculate variable genes ----
    data1 <- FindVariableFeatures(data1, nfeatures = 2000)
    VariableFeaturePlot(data1)
    
    list=  NULL
    vargenelist <- SplitObject(data1, split.by = "orig.ident")
    for (i in 1:length(vargenelist)) {
      vargenelist[[i]] <- FindVariableFeatures(vargenelist[[i]],nfeatures = 1000, verbose = FALSE)
    }
    {
      for (i in 1:length(vargenelist)) {
        x <- vargenelist[[i]]@assays$RNA@var.features
        list=c(list,x)}
      list=unique(list)
      length(list)
      data1@assays$RNA@var.features = list
      
    }
    data1 <- ScaleData(data1,split.by = 'orig.ident')
    data1 <- RunPCA(data1, pcs.compute = 50)
    data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident',)
    # ElbowPlot(data1,reduction = 'harmony',ndims = 50)+ElbowPlot(data1,ndims=50)
    d = 1:20
    
    data1 <- RunUMAP(data1, dims = d,
                     reduction = 'harmony',
                     reduction.name ='umap.h',reduction.key ='umap.h',
                     n.neighbors = 20L, #how many similar cells do you expect
                     spread =1, #how big is the axis space: 0.2 for lineages; 1 for separation
                     min.dist = 0.2, #how close to plot each cell higher=spread
                     local.connectivity = 100,
                     seed.use = 42)#overall how connected is the graph
    
  }
  ##Seurat integration ---- 
  {
    data1@active.assay='integrated'
    data1 <- ScaleData(data1, verbose = FALSE)
    data1 <- RunPCA(data1,
                    features = data1@assays$integrated@var.features,
                    npcs = 50, verbose = FALSE)
    PCAPlot(data1)
    
    # Determine percent of variation associated with each PC
    pct <- data1[["pca"]]@stdev / sum(data1[["pca"]]@stdev) * 100
    
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    
    # choose the minimum of these two metrics as the PCs covering the majority of the variation in the data
    pcs <- min(co1, co2)
    d = 1:pcs
    print(pcs)
    data1 <- RunUMAP(data1, reduction = "pca", 
                     dims = d,
                     n.neighbors = 20L, #how many similar cells do you expect save==10
                     spread =0.6, #how big is the axis space: 0.2 for lineages; 1 for separation
                     min.dist = 0.2, #how close to plot each cell higher=spread
                     local.connectivity = 10,#overall how connected is the graph save==1
                     reduction.name = 'umap.int',reduction.key = 'umap.int',
                     verbose = F,
                     seed.use = 10)
    
    # DimPlot(data1,group.by = 'orig.ident',cols=alpha(LibCP,1))&NoAxes()+NoLegend()&labs(title='medusa top')&
    #   DimPlot(data1,group.by = 'orig.ident',cols=alpha(rev(LibCP),1),
    #           order=levels(as.factor(data1$orig.ident)))+labs(title='polyp top')&NoAxes()
    # # 
    data1 <- FindNeighbors(data1, 
                           reduction = "pca", 
                           dims = 1:20,
                           annoy.metric = 'cosine',
                           verbose=F,
                           k.param = 20) 
    
    data1 <- FindClusters(data1, verbose=F,
                          resolution = 0.1) 
    data1 <- BuildClusterTree(data1,
                              reorder = TRUE,
                              reorder.numeric = T,verbose = F)
    data1$ID.cluster=data1@active.ident
  }
  
  Aa.Alldata=data1
  data1@active.assay = 'RNA'
  data1<-SetIdent(data1,value = 'ID.cluster')
}
#name clusters ID ----
{
  clusterNames<- read_excel('DataS3.xlsx',
                            sheet = 'AGC.markers.integrated')

  goi = Aa_genes$gene_short_name[match(clusterNames$marker,Aa_genes$name)]
  #assign cluster ID to the individual libraries
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
  #Note: here we merge the cnidocytes
  
  #use the wanted order from the spreadsheet to re-order the clusters:
  #first order the identities..
  data1@active.ident = factor(data1@active.ident,
                              levels(data1@active.ident)[order(clName)])
  #set the names to your IDs; once the gene list works, it is really easy to then update your cluster names with your excel spreadsheet.
  levels(data1@active.ident) = clusterNames$ID[clName][order(clName)]
  #save the IDs in metadata:
  data1@meta.data$IDs = data1@active.ident
  data1@active.assay = 'integrated'

  Aa.Alldata=data1
  save(Aa.Alldata,file='Ac_manuscript_final/AaAlldata.Robj')
}

# DEGs ----
{  
  data1=SetIdent(Aa.Alldata,value='IDs')
  {
    data1@active.assay='RNA'
    genes=Aa_genes
    all.markers <- FindAllMarkers(data1,
                                  return.thresh = 0.001,
                                  min.pct = 0.2,
                                  only.pos = TRUE,
                                  verbose = F)
    
    all.markers[,8:11] = 'NA'
    names(all.markers)[8:11]=names(genes)[c(3,6:8)]
    # add GO terms associated with this list:
    for (i in 1:length(levels(data1@active.ident))) # 
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),7]
      x2=genes[match(x,genes$gene_short_name),c(3,6:8)]
      all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),8:11]<-x2
    }  
    
  }
  all.markers.Alldata=all.markers 
  #TFs
  {all.markers <- FindAllMarkers(data1,
                                 return.thresh = 0.001,
                                 features = intersect(AaTF_list,rownames(data1)),
                                 min.pct = 0.05,
                                 only.pos = TRUE,
                                 verbose = F)
    
    all.markers[,8:11] = 'NA'
    names(all.markers)[8:11]=names(genes)[c(3,6:8)]
    # add GO terms associated with this list:
    for (i in 1:length(levels(data1@active.ident))) # 
    {
      x=all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),7]
      x2=genes[match(x,genes$gene_short_name),c(3,6:8)]
      all.markers[as.numeric(all.markers$cluster)==i,][1:length(which(as.numeric(all.markers$cluster)==i)),8:11]<-x2
    }  }
  all.markers.Alldata.TF=all.markers
  save(all.markers.Alldata,file = 'Ac_manuscript_final/AllData.ID.markers.RObj')
  write.csv(all.markers.Alldata,file='Ac_manuscript_final/DS1.3.csv')
  save(all.markers.Alldata.TF,file = 'Ac_manuscript_final/AllData.ID.markers.TF.RObj')
  write.csv(all.markers.Alldata.TF,file='Ac_manuscript_final/DS1.3b.csv')
}

# topGO ----
library(topGO)

Unfold <- Aa_annotations %>% 
  dplyr::mutate(`GOs` = strsplit(as.character(`GOs`), ",")) %>% 
  tidyr::unnest(`GOs`) 
geneID2GO <- Unfold %>% split(x = .$`GOs`, f = .$gene_short_name)
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
  myInterestingGenes[[i]] <- all.markers$gene[all.markers$cluster==unique(all.markers$cluster)[i]] #list of genes you want to perform GO enrichment for
  geneList[[i]] <- factor(as.integer(geneNames %in% myInterestingGenes[[i]]))
  names(geneList[[i]]) <- geneNames
  GOdata[[i]] <- new("topGOdata",ontology = "BP", allGenes = geneList[[i]],annot = annFUN.gene2GO, gene2GO = geneID2GO)
}
# filter and generate figures.
p=NULL
resultFis=NULL
topGO.IDs=NULL
for(i in 1:length(unique(all.markers$cluster)))
{
  #run the test...
  resultFis[[i]] <- runTest(GOdata[[i]], algorithm = "classic", statistic = "fisher") 
  pvalFis <- score(resultFis[[i]])
  #filter for only >0.05
  pvalFis = pvalFis[pvalFis>=0.05]
  allRes_intgenes<- GenTable(GOdata[[i]], pvalues = resultFis[[i]], orderBy = "pvalues", topNodes=30)
  allRes_intgenes$pvalues<-as.numeric(allRes_intgenes$pvalues)
  #convert NA values to zero; only if p value is so small it cannot be displayed in r
  allRes_intgenes[is.na(allRes_intgenes)]<-0.00000001
  
  #plot GOenrichment 
  
  david_col=LibCP[c(2,10,4,8)]
  topGO.IDs[[i]]=
    ggplot2::ggplot(allRes_intgenes, ggplot2::aes(x=(reorder(Term,(-log10(pvalues)))), y=(-log10(pvalues)))) +
    ggplot2::stat_summary(geom = "bar", fun = mean, position = "dodge",
                          col=david_col[i],fill=david_col[i]) +
    ggplot2::coord_flip()+
    ggplot2::xlab("Biological Process") +
    ggplot2::ylab("Enrichment -log10 p-value") +
    ggplot2::labs(title="GO enrichment",subtitle=levels(all.markers$cluster)[i])
  
  # print(topGO.IDs[[i]])
}

save(topGO.IDs,file='Ac_manuscript_final/topGO.IDs.RObj')
