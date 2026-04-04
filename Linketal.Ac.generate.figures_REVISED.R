# R script for generating figures in updated BioRxiv submission:
#Setup ----
#set library paths and load libraries:
.libPaths(c("/lisc/data/scratch/molevo/agcole/R/libs/seurat4/","/lisc/data/scratch/molevo/agcole/R/libs/course24/"))

setwd("/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/Ac_manuscript_revision_ACOE")

library(Seurat,quietly=T)
packageVersion('Seurat')
packageVersion('SeuratObject') #vs5 causes problems with everything!
library(RColorBrewer,quietly=T)
library(patchwork,quietly=T)
library(ggplot2,quietly=T)
library(pals,quietly=T)
library(readxl,quietly=T)
library(SeuratWrappers,quietly=T, lib.loc = "/lisc/data/scratch/molevo/agcole/R/libs/course24")
library(Nebulosa, lib.loc = "/lisc/data/scratch/molevo/agcole/R/libs/course24")
library(tidyr,quietly=T)
library(dplyr)
load (file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE.genes.RData')

if(!exists('Ac.Alldata'))
load(file='Ac.Alldata.Robj') 
if(!exists('data1.subsets'))
load(file='Ac.subsets.RObj')

genes=genes.ac

# how many genes have no data:
remove=genes.ac$gene_short_name[which(rowSums(Ac.Alldata@assays$RNA@counts)==0)]
length(genes.ac$gene_short_name)-length(remove)

# how many gene have at least 10 reads
keep=genes.ac$gene_short_name[which(rowSums(Ac.Alldata@assays$RNA@counts)>=10)]
length(keep) #17783
# how many genes have no annotation information?
genes=genes[keep%in%genes$gene_short_name,]
length(grep('ACOE',genes$gene_short_name))

# additional colour palettes: #check these with new data

gene.cp = c('lightgrey', rev(brewer.pal(11 , "Spectral")))
clust.cp.separate = unique (c(cols25(25), glasbey(32), alphabet(26)))
clust.cp.graded = unique(c(stepped3(20), stepped2(20), stepped(20)))
LibCP = c(stepped3(3),stepped(20)[c(13:16)],stepped(20)[c(11,9)],stepped3(20)[c(5,6)],'red',stepped(20)[c(1,2)],stepped2(20)[c(17,18)])
Ac.Alldata$orig.ident=as.factor(Ac.Alldata$orig.ident)
# Ac.Alldata$orig.ident=factor(Ac.Alldata$orig.ident,levels = levels(Ac.Alldata$orig.ident)[c(1,2,3,6,4,7,5,8:16)])
names(LibCP)=levels(as.factor(Ac.Alldata$orig.ident))
#brewer.paired(16)
LibCP.stages=LibCP[c(1,9,11,16)]
names(LibCP.stages)=c('polyp','strobila','ephyra','medusa')
neur.cp=rev(c('wheat3','burlywood',stepped(8)[7:8],"#ff7f00",stepped3(8)[5:7],'orange','chocolate','darkkhaki',stepped(12)[9:11],"#33a02c",'green3','forestgreen','darkgreen','darkolivegreen',(stepped2(8)[5:8])))

all.clusters.cp=c("skyblue","#1F78C8","darkblue","seagreen3","#FFD700" ,'goldenrod3','goldenrod1','yellow','goldenrod4',"grey",brewer.greys(14)[6:9],"#ff0000",brewer.reds(6)[3:6],neur.cp,"#FB6496",brewer.rdpu(20)[c(6,8,10,14,17,12,18)],'violet',"#6A33C2", brewer.purples(10)[4:10],'black')

names(all.clusters.cp)=levels(SetIdent(Ac.Alldata,value='ID.separate'))
names(neur.cp)=levels(data1.subsets$neural)
# 

#pal.bands(neur.cp,all.clusters.cp,LibCP)

# Figure 1 ----
## Fig. 1A insert ----
data1<- SetIdent(Ac.Alldata,value='orig.ident')
#generate plot with all cells, and only life history highlighted:
coi.samples=NULL
for(i in 1:length(levels(data1)))
{coi.samples[[i]]=WhichCells(data1,idents=levels(data1)[i])}
names(coi.samples)=levels(data1)

p1=DimPlot(data1,cells.highlight = coi.samples[1:7], cols.highlight=alpha((LibCP[1:7]),0.75),reduction = 'umap.h')+NoAxes()+NoLegend()
p2=DimPlot(data1,cells.highlight = coi.samples[8:9], cols.highlight=alpha((LibCP[8:9]),0.75),reduction = 'umap.h')+NoAxes()+NoLegend()
p3=DimPlot(data1,cells.highlight = coi.samples[10:11], cols.highlight=alpha((LibCP[10:11]),0.75),reduction = 'umap.h')+NoAxes()+NoLegend()
p4=DimPlot(data1,cells.highlight = coi.samples[12:16], cols.highlight=alpha((LibCP[12:16]),0.75),reduction = 'umap.h',order=levels(as.factor(Ac.Alldata$orig.ident)))+NoAxes()+NoLegend()

p3+p4+p2+p1+patchwork::plot_layout(ncol=2)

## Fig 1 C
# Calculate all genes for all cells
data1=SetIdent(Ac.Alldata,value = 'lifehistory')
keep =which(rowSums(data1@assays$RNA@counts) >= 10)
# coi.8000=WhichCells(SetIdent(Ac.Alldata,value='lifehistory'),
                    # downsample = 8000,seed = NULL)
data1@active.assay='RNA'
data1=subset(data1,features = genes.ac$gene_short_name[keep])#,cells = coi.8000)
genes=genes.ac[keep,]
polyp=WhichCells(data1,idents = levels(data1)[1])
strobila=WhichCells(data1,idents = levels(data1)[2])
ephyra=WhichCells(data1,idents = levels(data1)[3])
medusa = WhichCells(data1,idents = levels(data1)[4])
## Fig. 1C ----
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
                                            fill=LibCP.stages[c(1,4,2,3)],
                                            cex=2,disable.logging=F,sigdigs = 0,print.mode = 'percent')
  plot.new()
  grid::grid.draw(ST.venn.plot)
  ST.venn.plot.raw <- VennDiagram::venn.diagram(list(a1,a2,a3,a4),filename = NULL,
                                                category = c('polyp','medusa','strobila','ephyra'),
                                                fill=LibCP.stages[c(1,4,2,3)],
                                                cex=2,disable.logging=F,sigdigs = 0,print.mode = 'raw')

## Fig.1C ----  
  grid::grid.draw(ST.venn.plot.raw)
  grid::grid.draw(ST.venn.plot)
  goi.lifehistory.all=VennDiagram::calculate.overlap(list(a1,a2,a3,a4))
}

load(file='replicates8000.kostya.RObj') #not done yet

pdf('Fig1b.1.pdf')
boxplot(num.gene.all.rep,col=LibCP.stages)#+title('Number of Expressed Genes')
dev.off()
## Fig. 1B ----
boxplot(num.gene.all.rep,col=LibCP.stages)#+title('Number of Expressed Genes')

# everything is significant; but is this the right way to test this?
rm(replicates,num.gene.all.rep)

## Fig.1D medusa genes ----
medusa.specific=setdiff(unique(c(medusa.genes,ephyra.genes)),unique(c(polyp.genes,strobila.genes)))
medusa.specific=genes.ac[match(medusa.specific,genes.ac$gene_short_name),]
write.csv(medusa.specific,file='DS1.1a.csv')

medusa.specific = read.csv(file='DS1.1a.csv')
g=medusa.specific$gene_short_name
# DotPlot(Ac.Alldata,'RNA',features = g,group.by = 'lifehistory',col.min = 0) #sanity check; these ARE medusa-specific?
x.all=DotPlot(Ac.Alldata,'RNA',features = g,group.by = 'ID.separate')

specificity.index.filter=  x.all$data %>%
  group_by(features.plot) %>% 
  filter(pct.exp >=5) %>%
  count(pct.exp >=5)
gene.plot=as.character(specificity.index.filter$features.plot)
x=DotPlot(Ac.Alldata,assay='RNA',features=gene.plot,group.by = 'ID.separate')
# sort gene.plot according to clusters with dist and hclust to get gene orders:
data1<-ScaleData(Ac.Alldata,features=gene.plot,split.by = 'lifehistory',assay='RNA')
m=AverageExpression(data1,slot = 'scale.data',assay='RNA',features=gene.plot)
h=hclust(dist(m$RNA),method='ward.D2')
#Fig1E ----
DotPlot(Ac.Alldata,assay='RNA',features=rev(gene.plot[h$order]),scale.by='size' ,dot.min = 0.01,group.by = 'ID.separate', col.min = 0,cols = c('lightgrey','red'))+RotatedAxis() +FontSize(12,8)+ labs(title='Medusa-specific genes',subtitle='up-regulated in >=5% of cell ID')+theme(legend.position = 'bottom',legend.text = (element_text(size=12)),panel.grid = element_line('grey80',linewidth = 0.2))+coord_flip()
medusa.genes.f=rev(genes.ac[match(gene.plot[h$order],genes.ac$gene_short_name),])
library(dplyr)
specificity.index2=  x.all$data %>% #select either x.all for all 672 or x for just the filtered
  group_by(id) %>% 
  filter(pct.exp >0 & avg.exp.scaled >=0)%>% #& features.plot %in% goi1) 
  count(pct.exp >0) 

summary.medusa.exp=as.data.frame(levels(Ac.Alldata$ID.separate))
summary.medusa.exp$n = 0L
names(summary.medusa.exp)[1] = 'id'
summary.medusa.exp$n[summary.medusa.exp$id %in% specificity.index2$id]=specificity.index2$n
#Fig1D ----

specificity.index2=  x$data %>% #select either x.all for all 672 or x for just the filtered
  group_by(id) %>% 
  filter(pct.exp >0 & avg.exp.scaled >=0)%>% #& features.plot %in% goi1) 
  count(pct.exp >0) 

summary.medusa.exp=as.data.frame(levels(Ac.Alldata$ID.separate))
summary.medusa.exp$n = 0L
names(summary.medusa.exp)[1] = 'id'
summary.medusa.exp$n[summary.medusa.exp$id %in% specificity.index2$id]=specificity.index2$n
#normalize this by cluster size?
coi=NULL
coi.l=NULL
data1<-SetIdent(Ac.Alldata,value='ID.separate')
for(i in 1:length(levels(data1)))
{
  coi[[i]]=WhichCells(data1,idents=levels(data1)[i])
  coi.l[i]=length(coi[[i]])
}
barplot((summary.medusa.exp$n/coi.l)*45,names.arg = summary.medusa.exp$id,las=2,col=all.clusters.cp[summary.medusa.exp$id],ylim =c(0,10)  )+title('Filtered.medusa.genes, detection / cluster normalized')
x=DotPlot(Ac.Alldata,assay='RNA',features=gene.plot[h$order],group.by = 'ID.separate')

medusa.specific.filtered=x$data %>% group_by(id) %>%  filter(pct.exp >5 & avg.exp.scaled >=0.1)
names(medusa.specific.filtered)[3]='gene_short_name'
medusa.specific.filtered.list=merge(medusa.specific.filtered,genes.ac,all.x=T,all.y=F)

write.csv(medusa.specific.filtered.list,file='DS1.1b.csv')

## Fig. 2A dimplot ----
#highlight smooth muscle on main
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
#' **update** the script... saved into AcAlldata...
clust.cp.reduced=c('goldenrod2','#FB6496','seagreen3','purple3','grey','black',"#ff7f00",'red2','#1F78C8')
names(clust.cp.reduced)=levels(SetIdent(Ac.Alldata,value='reduced'))

DimPlot(SetIdent(Ac.Alldata,value='reduced'),cols=clust.cp.reduced,pt.size = 0.75,order=rev(levels(data1)))+NoAxes()+NoLegend()
## Fig. 2C.distplots ----
{
  #summarize the data of interest:
  ids.cluster.sample = as.data.frame(table(Ac.Alldata@meta.data$IDs, Ac.Alldata@meta.data$lifehistory.tissue))
  colnames(ids.cluster.sample) = c('ID', 'sample', 'CellCount')
  
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
    ),legend.position = 'right') +
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
  
  
  
}
dist.clust2

{
  #summarize the data of interest:
  ids.cluster.sample = as.data.frame(table(Ac.Alldata@meta.data$ID.separate, Ac.Alldata@meta.data$lifehistory.tissue))
  colnames(ids.cluster.sample) = c('ID', 'sample', 'CellCount')
  
  dist.clust =
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
    scale_fill_manual(values = all.clusters.cp) +
    theme(axis.text.x = element_text(
      #face="bold", color="#993333",
      size = 8,
      angle = -45,
      hjust = 0,
      vjust = 0.5
    ),legend.position = 'right') +
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
  
  
  
}
dist.clust
# Figure 2B ----
#use PCs to build the tree:

embeddings <- Embeddings(object = Ac.Alldata, reduction = 'harmony')[,1:50]
Ac.Alldata<-SetIdent(Ac.Alldata,value='ID.separate')
data.dims=matrix('0',50L,length(levels(Ac.Alldata)))
for (i in 1:length(levels(Ac.Alldata)))  
{  cells <- WhichCells(object = Ac.Alldata, idents = levels(Ac.Alldata)[i])
data.dims[,i] <- colMeans(embeddings[cells,])
}
colnames(x = data.dims) <- levels(x = Ac.Alldata)
data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
nj.tree <- ape::nj(data.dist)
ape::plot.phylo(nj.tree, type = "fan", label.offset = 0.1,show.node.label = F,use.edge.length = F,lab4ut = 'axial',tip.color = all.clusters.cp, cex = 1,no.margin = T)

# Figure 3 cnido ----
  ## Fig3A,b ----
clust.cp=c("#FB6496",brewer.rdpu(20)[c(8:11,20,17,12,16)])
clust.cp[2]='lightpink'
clust.cp[3]='green4'
clust.cp[4]='slateblue'
clust.cp[5]='slateblue4'
clust.cp[8]='black'
clust.cp[9]='grey'
names(clust.cp)=levels(data1.subsets$cnidocyte)
cnido.clust=DimPlot(data1.subsets$cnidocyte,cols=clust.cp)+NoAxes()+NoLegend()+scale_x_reverse()
data1=data1.subsets$cnidocyte
ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory.tissue))
colnames(ids.cluster.library) = c('ID','Stage','CellCount')

print(cnido.clust+
        ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                        x=Stage)) +
        geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                              x=(Stage)),
                 position="fill", stat="identity", width = 0.5)+
        scale_fill_manual(values = clust.cp)+
        theme(axis.text.x = element_text(#face="bold", color="#993333", 
          size=8, angle=-45,hjust=0,vjust = 0.5))+
        geom_area(mapping =aes(fill=ID, y= log(CellCount),
                               x=as.integer(Stage)),
                  position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
        geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                              x=(Stage)),
                 position="fill", stat="identity", width = 0.5)
)
coi=NULL
for(c in 1:length(levels(data1.subsets$cnidocyte)))
  coi[[c]]=WhichCells(data1.subsets$cnidocyte,idents=levels(data1.subsets$cnidocyte)[c])
names(coi)=levels(data1.subsets$cnidocyte)
DimPlot(Ac.Alldata,cells.highlight = coi,cols.highlight = rev(clust.cp[sort(names(coi))]))&NoAxes()&NoLegend()

##Fig 3C ----
# load('../OMA.genes.RData')
{
  # load(file='../cnido.nv.ac.RData')
  load('../ac.kostyaPlus/oma.nv.ac.cnido.RData')
  load (file='/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/ACOE.genes.RData')
  data1<-SetIdent(data1,value = 'IDs')
  embeddings <- Embeddings(object = data1, reduction = 'harmony')[,1:20]
  data.dims=matrix('0',20L,length(levels(data1)))
  for (i in 1:length(levels(data1)))  
  {  cells <- WhichCells(object = data1, idents = levels(data1)[i])
  data.dims[,i] <- colMeans(embeddings[cells,])
  }
   library(ape)
  colnames(x = data.dims) <- levels(x = data1)
  data.dist <- dist(x = t(x = data.dims),method = 'maximum')
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
  ape::plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = T,use.edge.length = F,lab4ut = 'axial', cex = 1,no.margin = T, rotate.tree = 120)
  for (i in 1:length(node_col)) {
    nodelabels(node = length(nj.tree$tip.label)+i, pch=21, col="black", bg=node_col[i], cex=2)
  }
  
}
## Fig3 E ----
library("Nebulosa")
load('../../Genes.Nv2.RData')
temp@active.assay='RNA'
goi.a=c("myc1.1","soxc", "prd13-like-1", "zn431-like-1", 'pax2a-like-1','myc3.1', "jun-like3",'fos-like-4', "gfi1b-like-1", "foxn4-like-1", 'sry-like-1','mycn.1')#
myc.a=genes.ac[match(c('scaffold17.g336','scaffold17.g338','scaffold17.g606',
                                       'scaffold17.g351','scaffold21.g619','scaffold14.g112'),genes.ac$cellranger),]
check=genes.ac[match(goi.a,genes.ac$gene_short_name.y),]
check$gene_short_name[6]='mycn'
  plot_density(data1.subsets$cnidocyte,check,reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')&scale_x_reverse()
  plot_density(data1.subsets$cnidocyte, (myc.a$gene_short_name),reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')&scale_x_reverse()
  
  myc.n=genes$gene_short_name[match(c('myc','myc5','myc1','Myc2','Myc3','myc4'),genes$gene_short_name)]
  goi.n=c('myc1',
          'SoxC',
          'Prdm13','ZNF845','PAX2A-like-1','JUN-like-11','FOS-like-4','GFI1B-like-1','FoxN1','SoxA','myc','Myc3')
  goi.n=genes$name[match(goi.n,genes$gene_short_name)]
  plot_density(temp.n, goi.n,reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')
  plot_density(temp.n, myc.n,reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')
  
# Fig4 neurons----
#all.neuron
{
  coi=NULL
  for(c in 1:length(levels(data1.subsets$neural)))
    coi[[c]]=WhichCells(data1.subsets$neural,idents=levels(data1.subsets$neural)[c])
  names(coi)=levels(data1.subsets$neural)
  DimPlot(Ac.Alldata,cells.highlight = coi,cols.highlight = rev(neur.cp[sort(names(coi))]))&NoAxes()&NoLegend()
 
  data1<-data1.subsets$neural
  data1$orig.ident=as.factor(data1$orig.ident)
  data1$lifehistory.tissue=as.factor(data1$lifehistory.tissue)
  # data1$lifehistory.tissue=factor(data1$lifehistory.tissue,levels(data1$lifehistory.tissue)[c(1,2,4,5:8,3)])
  DimPlot(data1,cols=neur.cp,label=F)&NoAxes()&NoLegend()
    ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$orig.ident))#lifehistory.tissue))
    colnames(ids.cluster.library) = c('ID','Library','CellCount')
    ids.cluster.library$ID<-factor(ids.cluster.library$ID,levels(ids.cluster.library$ID)[rev(1:24)])
    
    clust.cp=rev(neur.cp)
    print(
      ggplot(ids.cluster.library, aes(fill=ID, y= (CellCount),
                                      x=Library)) +
        geom_bar(mapping =aes(fill=ID, y= (CellCount),
                              x=(Library)),
                 position="fill", stat="identity", width = 0.5)+
        scale_fill_manual(values = clust.cp)+
        theme(axis.text.x = element_text(#face="bold", color="#993333", 
          size=8, angle=-45,hjust=0,vjust = 0.5))+
        geom_area(mapping =aes(fill=ID, y= (CellCount),
                               x=as.integer(Library)),
                  position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
        geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                              x=(Library)),
                 position="fill", stat="identity", width = 0.5)+labs(title = 'neurons')
    )

all.neurons=data1.subsets$neural      
    coi=NULL
    for(c in 1:length(levels(all.neurons)))
      coi[[c]]=WhichCells(all.neurons,idents=levels(all.neurons)[c])
    names(coi)=levels(all.neurons)
    DimPlot(Ac.Alldata,cells.highlight = coi,cols.highlight = rev(neur.cp[sort(names(coi))]))&NoAxes()&NoLegend()
    all.neurons$IDs=as.factor(all.neurons$IDs)
    goi=c('ashA','pou4','ins1a-like-1','soxc')
    ### new voltage figure ----
    voltage.juan<-read.delim('../ac.kostyaPlus/Voltage-dependent_channel_IPR003938_IPR027359_IPR014743_IPR016449 (1).lst',header = F)
    voltage.juan<-genes.ac[match(voltage.juan$V1,genes.ac$cellranger),]
    voltage.juan.goi<-genes.ac$gene_short_name[match(voltage.juan$cellranger,genes.ac$cellranger)]

    #filter by expression
    voltage.juan.goi<-voltage.juan.goi[rowSums(data1.subsets$neural@assays$RNA@counts[voltage.juan.goi,])>=10]
    voltage.juan.goi.f<-voltage.juan.goi[rowSums(data1.subsets$neural@assays$RNA@counts[voltage.juan.goi,])>=250]
    voltage<-genes.ac[match(voltage.juan.goi,genes.ac$gene_short_name),][c(17,31,32,39,42,53),]
    DotPlot(data1,'RNA',voltage.juan.goi[c(17,31,32,39,42,53)],cols=c('lightgrey','red'),scale.by = 'size',dot.min = 0.25,col.min = 0,scale = F)&RotatedAxis()&theme(legend.title = element_text(size = 8),legend.text = element_text(size=8),panel.grid = element_line(linewidth = 0.2,colour = 'grey90'))
    DotPlot(Ac.Alldata,'RNA',c(voltage.juan.goi[c(17,31,32,39,42,53)]),group.by='ID.separate',cols=c('black','red'),scale.by = 'size',dot.min = 0.25,scale = T,col.min = 0)&RotatedAxis()&theme(legend.title = element_text(size = 8),legend.text = element_text(size=8),panel.grid = element_line(linewidth = 0.2,colour = 'grey90'), axis.text.x = element_text(size = 8))&coord_flip()
    DotPlot(data1,'RNA',c('ins1a-like-1','soxc'),cols=c('lightgrey','darkorange3'),scale.by = 'size',dot.min = 0.01,col.min = 0,scale = F)&RotatedAxis()&
    DotPlot(data1,'RNA',c('ashA','pou4'),cols=c('lightgrey','darkgreen'),scale.by = 'radius',dot.min = 0.05,col.min = 0,scale = F)&RotatedAxis()&theme(legend.title = element_text(size = 8),legend.text = element_text(size=8),panel.grid = element_line(linewidth = 0.2,colour = 'grey90'))
    
### tree ----
    embeddings <- Embeddings(object = data1, reduction = 'harmony')[,1:50]
    data.dims=matrix('0',50L,length(levels(data1)))
    for (i in 1:length(levels(data1)))  
    {  cells <- WhichCells(object = data1, idents = levels(data1)[i])
    data.dims[,i] <- colMeans(embeddings[cells,])
    }
    colnames(x = data.dims) <- levels(x = data1)
    data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
    nj.tree <- ape::nj(data.dist)
    ape::plot.phylo(ape::root(nj.tree,outgroup = 'n.early.state'), type = "unrooted", label.offset = 0.5,show.node.label = F,use.edge.length = F,lab4ut = 'axial',tip.color = rev(clust.cp), cex = 0.8,no.margin = T,rotate.tree = 180)
    
    
goi=c('tba2-like-1',
      'prfa-like-1')#,'tbb-like-2','tbb-like-7')#'prfa-like-1','prfa-like-2','prfa4-like'

DotPlot(data1,'RNA',goi,group.by='IDs',cols=c('white','black'),scale.by = 'size',dot.min = 0.25,scale = F)&RotatedAxis()&theme(legend.title = element_text(size = 8),legend.text = element_text(size=8),panel.grid = element_line(linewidth = 0.2,colour = 'grey90'), axis.text.x = element_text(size = 1))&coord_flip()
# Figure 6 muscle ----
load(file='New.allmuscle.RObj')
features = genes.ac$gene_short_name[match(c("ACOE011G013993","ACOE016G019902","ACOE013G016778","ACOE013G016139"),genes.ac$ACO)]
##Fig 6Q ----
FeaturePlot(allmuscle, features = features, reduction = "umap", pt.size = 2, cols = c("lightgrey", "blue"),order = T,max.cutoff = 3,min.cutoff = 1,ncol=4)&NoAxes()&NoLegend()

#####change the order of the stages 
allmuscle$stage <- factor(
  allmuscle$lifehistory,
  levels = c("polyp", "strobila", "ephyra", "medusa"))

df <- as.data.frame(table(
  allmuscle@meta.data$lifehistory.tissue,
  allmuscle@meta.data$IDs
))
colnames(df) <- c("Library", "CellType", "CellCount")

colors <- c("muscle.st" = "#1F78C8", "smooth.muscle" = "seagreen3")

## Fig6RS ----
coi[[1]]=WhichCells(allmuscle,idents = levels(allmuscle)[1])
coi[[2]]=WhichCells(allmuscle,idents = levels(allmuscle)[2])
names(coi)=levels(allmuscle)
DimPlot(Ac.Alldata,cells.highlight = coi,cols.highlight =(c("seagreen3","#1F78C8")))&NoLegend()&NoAxes()
  (DimPlot(allmuscle,cols=(c("#1F78C8","seagreen3")),reduction = 'umap')&NoAxes() |
ggplot(df, aes(x = Library, y = CellCount, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
    legend.title = element_blank()
  )  +theme(legend.position = "none")+ylab('Percent library contribution')+xlab('Life cycle stage'))+plot_layout(widths = c(1, 3))

  
  # Figure 7 muscle ----
Ac.Alldata<-SetIdent(Ac.Alldata,value='IDs')
coi.st=WhichCells(Ac.Alldata,idents ='muscle.st')

Ac.Alldata<-SetIdent(Ac.Alldata,value='ID.separate')
coi.sm=WhichCells(Ac.Alldata,idents ='smooth.muscle')

data1=Ac.Alldata
Idents(data1) <- data1$IDs
data1$reduced=as.character(data1@active.ident)
data1$reduced[coi.sm]=as.character(Ac.Alldata$ID.separate[coi.sm])
data1<-SetIdent(data1,value='reduced')

#further reduce:
levels(as.factor(data1$reduced))
new.names=c('non-muscle','non-muscle','non-muscle','non-muscle','non-muscle','muscle.st','non-muscle','muscle.sm','non-muscle')
data1$reduced.muscle=as.factor(data1$reduced)
levels(data1$reduced.muscle)=new.names
data1$reduced.muscle=factor(data1$reduced.muscle,levels(data1$reduced.muscle)[c(2,3,1)])

#generate muscle protein gene lists:
calm=genes.ac$gene_short_name[grep('calm',genes.ac$gene_short_name)]
actin=genes.ac$gene_short_name[grep('Actin',genes.ac$PFAMs)]
myosin<-genes.ac$gene_short_name[grep('Myosin_',genes.ac$PFAMs)]
myos.l<-genes.ac$gene_short_name[grep('Myosin, light',genes.ac$Description)]
melc=c(myos.l[-6],'mlc2-like-5','mlc2-like-1')
mlck<-genes.ac$gene_short_name[grep('light chain kinase',genes.ac$Description)]
mlck=c(mlck,'scaffold16.g60','scaffold16.g61')#,'obscurin-like'
# 
melcPh<-genes.ac$gene_short_name[grep('myosin-light-chain-phosphatase',genes.ac$Description)]
melcPh=c(melcPh,'mypt2-like-1')
mylia<-genes.ac$gene_short_name[grep('myosin regulatory',genes.ac$Description)]

tpm = genes.ac$gene_short_name[grep('Tropomyosin',genes.ac$PFAMs)]
myoph<-genes.ac$gene_short_name[grep('Calponin',genes.ac$PFAMs)]
myoph<-c(myoph,'smtl1-like')
mrlc<-genes.ac$gene_short_name[grep('GKAP',genes.ac$PFAMs)]

goi.all=c(actin,myosin,melc,mlck,mrlc,melcPh,mylia,myoph,tpm,calm)
goi1=unique(goi.all)
#generate the plot:
x=DotPlot(data1,features=unique(goi1),group.by = 'reduced')
x=x$data
x=x[!is.nan(x$avg.exp.scaled),]
#subset the muscle:
x2a=x[x$id == 'muscle.st' | x$id == 'smooth.muscle',]
#filter for expression min. cells
x2=x2a[x2a$pct.exp>=6,]#was 10; tmp3 is only 6

#filter genes for up-regulation in muscle
goi=x2$features.plot[x2$avg.exp.scaled>0.5]
goi=unique(as.character(goi))
#filter gene sets:

goi.plot<-goi1[goi1 %in% goi]
check<-genes.ac[genes.ac$gene_short_name %in% goi.plot,]

#drop the low expressed actins:
goi.plot=goi.plot[c(-1,-2,-5,-6)]
goi.plot=goi.plot[c(1:2,6:10,3:5,11:33)]
#Plot with names
p=DotPlot(data1,features=rev(goi.plot),scale.by='size' ,dot.min = 0.01,group.by = 'reduced.muscle', col.min = -0.5,
          cols = c('grey90','red'))+RotatedAxis() +FontSize(12,10)+ labs(title='PFAM/Orthogroup')+theme(legend.position = 'bottom',legend.text = (element_text(size=6)),panel.grid = element_line('grey',linewidth = 0.1))+coord_flip()

#update names to gene models:
gene.names=p$data$features.plot
for (i in 1:length(levels(gene.names)))
  levels(gene.names)[i]=genes.ac$ACO[match(as.character(levels(gene.names)[i]), genes.ac$gene_short_name)]
p2=p
p2$data$features.plot=gene.names
#sanity check that re-naming worked correctly:
p+p2 
# plot Figure8:
p2

# Generate CellBrowser Files ----
Ac.Alldata@active.assay='RNA'
AllData <- DietSeurat(Ac.Alldata,data=T,assays='RNA',dimreducs = c('umap.int'))
AllData@meta.data = AllData@meta.data[c(1:5,10,11)]
head(AllData@meta.data)
DimPlot(AllData,label=T,cols=brewer.paired(14),label.size = 4)&NoAxes()&NoLegend()
save(AllData,file = 'TEMP_CellBrowser/AllData.Robj')

my.list=dput(names(data1.subsets))

for (i in 1:length(names(data1.subsets)))
{
  data1=data1.subsets[[i]]
  data1@active.assay='RNA'
  data1 <- DietSeurat(data1,data=T,assays='RNA',dimreducs = c('umap'))
  data1@meta.data = data1@meta.data[,c(1:5,9,10)]
  # head(data1@meta.data)
  data1<-SetIdent(data1,value = 'IDs')
  assign(my.list[i],data1) 
  save (list=my.list[i],file = paste('TEMP_CellBrowser/',my.list[i],'.RObj',sep='')) 
  #clean up
  rm (list=my.list[i])
  #works!!!
}


for (i in 1:length(my.list))
{
  jpeg(file=paste('TEMP_CellBrowser/',my.list[i],'.jpg',sep=''))
  print(DimPlot(Ac.Alldata,cells.highlight = colnames(data1.subsets[[i]]),cols.highlight = 'red',reduction='umap.h')&NoAxes()&NoLegend()&labs(title=my.list[i]))
  dev.off()
}

jpeg(file=paste('TEMP_CellBrowser/SeriesOverview.JPEG'),)
print(DimPlot(Ac.Alldata,cols=LibCP,group.by = 'orig.ident',pt.size = 0.5)&NoAxes()&labs(title = 'All Data'))
dev.off()

write.csv(genes.ac[,1:2],'TEMP_CellBrowser/gene.table.csv')

# generate Supplement S1 ----
medusa.specific=read.csv(file='DS1.1a.csv')
xlsx::write.xlsx(medusa.specific,'DataS1.xlsx',sheetName="DS1.1a.medusa.specific",col.names = T,row.names = F)
medusa.specific.filtered.list=read.csv(file='DS1.1b.csv')
xlsx::write.xlsx(medusa.specific.filtered.list,'DataS1.xlsx',sheetName="DS1.1b.medusa.specific.filtered",col.names = T,row.names = F,append=T)
load(file='markers.lifehistory.RObj')
xlsx::write.xlsx(life.history.markers,'DataS1.xlsx',sheetName="DS1.2.lifecycle",col.names = T,row.names = F,append=T)
load(file='TF.markers.lifehistory.RObj')
xlsx::write.xlsx(life.history.markers.TF,'DataS1.xlsx',sheetName="DS1.2b.lifecycle.TFs",col.names = T,row.names = F,append=T)
load(file = 'AllData.ID.markers.RData')
xlsx::write.xlsx(all.markers.Alldata,'DataS1.xlsx',sheetName="DS1.3.IDs",col.names = T,row.names = F,append=T)
xlsx::write.xlsx(all.markers.Alldata.TF,'DataS1.xlsx',sheetName="DS1.3b.IDs.TFs",col.names = T,row.names = F,append=T)
load(file = 'Ac.AllData.subsets.markers.RObj')
load(file = 'Ac.AllData.subsets.TF.markers.RObj')
for (i in 1:7)
{
  xlsx::write.xlsx(markers[[i]],'DataS1.xlsx',sheetName=paste('DS1',i+3,names(markers)[i],sep='.'),col.names = T,row.names = F,append=T)
  xlsx::write.xlsx(markers.TF[[i]],'DataS1.xlsx',sheetName=paste('DS1',i+3,"b",names(markers)[i],'TFs',sep='.'),col.names = T,row.names = F,append=T)
}
load('DataS2.RData')
DataS2<-rbind(DataS2.sh,DataS2.sm,DataS2.st)
DataS2=DataS2[,c(-7,-8,-9,-10,-11)]
xlsx::write.xlsx(DataS2,'DataS1.xlsx',sheetName='DataS1.13.muscle.specific',col.names = T,row.names = F,append=T)

##DataS2 ----
{
  genes.file=genes.ac[,c(-7,-8,-11)]
xlsx::write.xlsx(genes.file,'DataS2raw.xlsx',sheetName="DS2.2.GeneAnnotations",col.names = T,row.names = F)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet = 'AGC.markers.integrated')
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[,c(2,4,6)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName="DS2.3.Partitions",col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet = 'OuterClusters')
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1:6),c(2,3,5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName="DS2.4.Epidermis",col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet = 'InnerClusters') 
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1,2,3,4,7),c(2,3,5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName='DS2.5.Gastrodermis',col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet='Dig.Gland') 
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1:8),c(2,3,5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName='DS2.6.Digestive.Glands',col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet='MucinClusters') 
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1,4:7),c(2,3,5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName='DS2.7.Mucin.Gland',col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet='CnidoClusters') 
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1:9),c(2,3,5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName='DS2.8.Cnidocytes',col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet='Neuro') 
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1:23),c(3:5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName='DS2.9.All.Neuroglandular',col.names = T,row.names = T,append=T)

clusterNames<- read_excel('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/DataS3.kostyaPlus.xlsx',
                          sheet = 'Muscle.ST')
clusterNames$marker=genes.ac$gene_short_name[match(clusterNames$marker,genes.ac$gene_short_name.x)]
clusterNames=clusterNames[c(1,2,4),c(2,3,5)]
xlsx::write.xlsx(clusterNames,'DataS2raw.xlsx',sheetName='DS2.10.Muscle.ST',col.names = T,row.names = T,append=T)

OMA<- readxl::read_xlsx('/lisc/data/scratch/molevo/agcole/R/Aurelia_51k/ac.kostyaPlus/OMAOrthologousGroupsAcNv.xlsx')
lut=read.delim('/lisc/data/scratch/molevo/jmontenegro/alison/aaurita/results/annotation/kostya+/geneID_AAUR2.map.txt',header = F)
names(lut)=c('geneID','AAUR2')

#update with the gene_ model because the names had changed:
ind=match(OMA$Ac,lut$AAUR2,nomatch = 0)
OMA$Ac[1:42]=lut$geneID[ind]
OMA <- OMA %>%
  dplyr::left_join(genes.ac %>% dplyr::select(cellranger, ACO,Description,gene_short_name),
            by = c("Ac" = "cellranger"))
OMA=OMA[,c(4,6,5,3,1)]
xlsx::write.xlsx(OMA,'DataS2raw.xlsx',sheetName='DS2.11.Ac_Nv.OMA',col.names = T,row.names = T,append=T)
}
