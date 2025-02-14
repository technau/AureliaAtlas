# R script for generating figures in updated BioRxiv submission:

library(Seurat,quietly=T)
packageVersion('Seurat')
library(RColorBrewer,quietly=T)
library(patchwork,quietly=T)
library(ggplot2,quietly=T)
library(pals,quietly=T)
library(readxl,quietly=T)
library(SeuratWrappers,quietly=T)
library(tidyr,quietly=T)
library(Matrix)
library(dplyr)

load (file="Aaur2.newnames.RData")
LibCP =brewer.paired(12)
LibCP.stages=LibCP[c(2,4,5,8)]

if(!exists('Aa.Alldata'))
load(file='Ac_manuscript_final/AaAlldata.Robj')
if(!exists('data1.subsets'))
load(file='Ac_manuscript_final/Aa.subsets.RObj')

genes=Aa_genes

# how many genes have no data:
remove=Aa_genes$Cellranger.Aaur2[which(rowSums(Aa.Alldata@assays$RNA@counts)==0)]
length(Aa_genes$Cellranger.Aaur2)-length(remove)

# how many gene have at least 10 reads
keep=Aa_genes$Cellranger.Aaur2[which(rowSums(Aa.Alldata@assays$RNA@counts)>=10)]
length(keep) #34301
# how many genes have no annotation information?
genes=genes[keep%in%genes$Cellranger.Aaur2,]
length(grep('gene.',genes$gene_short_name)) #23959

# additional colour palettes: #check these with new data
neur.cp=c("#ff7f00",'orange','chocolate',stepped3(8)[7:8],'orange3',stepped(8)[5:7],'darkkhaki',stepped2(8)[5:8],stepped(12)[10:12],"#33a02c",'green3')

all.clusters.cp=c("#1F78C8","deepskyblue2","#a6cee3",'lightblue3','aquamarine4','darkblue',"darkolivegreen4","#FFD700" ,'lightgoldenrod1','#e7a300','goldenrod1','black',"#565656",brewer.greys(20)[c(8,10,12,14)],"#ff0000",brewer.reds(20)[c(5,8,10,15)],"#6A33C2", brewer.purples(15)[6:12],neur.cp,"#FB6496",brewer.rdpu(20)[c(8:11,20,17,12,16)])

names(all.clusters.cp)=levels(SetIdent(Aa.Alldata,value='ID.separate'))
names(neur.cp)=c(levels(data1.subsets$neural.1),levels(data1.subsets$neural.2))

# Figure 1 ----
## Fig. 1A insert ----
DimPlot(Aa.Alldata,group.by = 'orig.ident',cols=alpha(LibCP,1))+NoAxes()+theme(legend.position = 'bottom')&labs(title='medusa top')&
  DimPlot(Aa.Alldata,group.by = 'orig.ident',
          order=levels(as.factor(Aa.Alldata$orig.ident)),cols=alpha(rev(LibCP),1))+labs(title='polyp top')+NoAxes()+NoLegend()

## Fig 1 C
# Calculate all genes for all cells
data1=SetIdent(Aa.Alldata,value = 'lifehistory')
keep =which(rowSums(data1@assays$RNA@counts) >= 10)

data1@active.assay='RNA'
data1=subset(data1,features = Aa_genes$gene_short_name[keep])#,cells = coi.8000)
genes=Aa_genes[keep,]
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
                                            fill=LibCP[c(2,8,4,5)],
                                            cex=2,disable.logging=F,sigdigs = 0,print.mode = 'percent')
  plot.new()
  grid::grid.draw(ST.venn.plot)
  ST.venn.plot.raw <- VennDiagram::venn.diagram(list(a1,a2,a3,a4),filename = NULL,
                                                category = c('polyp','medusa','strobila','ephyra'),
                                                fill=LibCP[c(2,8,4,5)],
                                                cex=2,disable.logging=F,sigdigs = 0,print.mode = 'raw')

## Fig.1C ----  
  grid::grid.draw(ST.venn.plot.raw)
  grid::grid.draw(ST.venn.plot)
  goi.lifehistory.all=VennDiagram::calculate.overlap(list(a1,a2,a3,a4))
}

load(file='Ac_manuscript/replicates8000.StrEph.RObj')

pdf('Ac_manuscript/Fig1b.1.pdf')
boxplot(num.gene.all.rep,col=LibCP[c(2,4,5,8)])#+title('Number of Expressed Genes')
dev.off()
## Fig. 1B ----
boxplot(num.gene.all.rep,col=LibCP[c(2,4,5,8)])#+title('Number of Expressed Genes')
rm(replicates,num.gene.all.rep)


## Fig. 1D ----
{
load(file='/lisc/scratch/molevo/agcole/R/Aurelia_51k/Ac_manuscript/markers.lifehistory.StrEph.RObj')
View(life.history.markers)
DotPlot(SetIdent(Aa.Alldata,value='lifehistory'),'RNA',rev(unique(life.history.markers$gene)),split.by = 'lifehistory',cols = LibCP[c(2,4,5,8)],scale.by = 'radius',
        col.min = 0,
        dot.scale = 5)&
  RotatedAxis()+FontSize(10,5)&theme(panel.grid = element_line(linewidth = 0.2,colour = 'grey90'),legend.title = element_text(size = 8),legend.text = element_text(size=6),legend.position = 'bottom')&coord_flip()

#generate the legends:
goi=(unique(life.history.markers$gene))
L_P=FeaturePlot(Aa.Alldata,goi[1],cols=c('lightgrey',LibCP[c(2)]))
L_S=FeaturePlot(Aa.Alldata,goi[1],cols=c('lightgrey',LibCP[c(4)]))
L_E=FeaturePlot(Aa.Alldata,goi[1],cols=c('lightgrey',LibCP[c(5)]))
L_M=FeaturePlot(Aa.Alldata,goi[1],cols=c('lightgrey',LibCP[c(8)]))
legP <- ggpubr::get_legend(L_P+theme(legend.position = 'top',legend.text.align = 1))
legS <- ggpubr::get_legend(L_S+theme(legend.position = 'top',legend.text.align = 1))
legE <- ggpubr::get_legend(L_E+theme(legend.position = 'top',legend.text.align = 1))
legm <- ggpubr::get_legend(L_M+theme(legend.position = 'top',legend.text.align = 1))

leg1=ggpubr::as_ggplot(legP)
leg2=ggpubr::as_ggplot(legS)
leg3=ggpubr::as_ggplot(legE)
leg4=ggpubr::as_ggplot(legm)
colourscales=leg1+leg2+leg3+leg4+plot_layout(ncol=1)
colourscales


}
## Fig. 1E.GO ----

load(file='Ac_manuscript_final/lifehistory.topGO.StrEph.RObj')
library(ggside)
# pdf('Fig1E.pdf',height = 8,width = 8)
# print(topGO.lifehistory[[1]]+topGO.lifehistory[[4]])
# dev.off()
topGO.lifehistory[[1]]+topGO.lifehistory[[4]]
# rm(topGO.lifehistory)

## Fig. 1e dimplot
sample.plot=DimPlot(SetIdent(Aa.Alldata,value='IDs'),cols=c(clust.cp.separate[c(1,7,13,2,4,5,3,9,15:30)]),pt.size = 0.75)+NoAxes()+NoLegend()

## Fig. 1f.barplots ----

clust.cp =c(clust.cp.separate[c(1,7,13,2,4,5,3,9,15:30)])
{
  #summarize the data of interest:
  ids.cluster.sample = as.data.frame(table(Aa.Alldata@meta.data$IDs, Aa.Alldata@meta.data$orig.ident))
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
    ),legend.position = 'bottom') +
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
sample.plot+dist.clust2

## Fig.1G medusa genes ----
medusa.specific=setdiff(c(medusa.genes,ephyra.genes),c(polyp.genes,strobila.genes))

medusa.specific=Aa_genes[match(medusa.specific,Aa_genes$gene_short_name),c(3:8)]
write.csv(medusa.specific,file='Ac_manuscript_final/DS1.1a.csv')

g=medusa.specific$gene_short_name
# DotPlot(Aa.Alldata,'RNA',features = g,group.by = 'lifehistory') #sanity check; these ARE medusa-specific?
x.all=DotPlot(Aa.Alldata,'RNA',features = g,group.by = 'ID.separate')

specificity.index.filter=  x.all$data %>%
  group_by(features.plot) %>% 
  filter(pct.exp >5) %>%
  count(pct.exp >5)
gene.plot=as.character(specificity.index.filter$features.plot)
x=DotPlot(Aa.Alldata,assay='RNA',features=gene.plot,group.by = 'ID.separate')
# sort gene.plot according to clusters with dist and hclust to get gene orders:
data1<-ScaleData(Aa.Alldata,features=gene.plot,split.by = 'lifehistory',assay='RNA')
m=AverageExpression(data1,slot = 'scale.data',assay='RNA',features=gene.plot)
h=hclust(dist(m$RNA),method='ward.D2')
#Fig1G ----
DotPlot(Aa.Alldata,assay='RNA',features=(gene.plot[h$order]),scale.by='size' ,dot.min = 0.05,group.by = 'ID.separate', col.min = 0,cols = c('lightgrey','red'))+RotatedAxis() +FontSize(12,8)+ labs(title='Medusa-specific genes',subtitle='up-regulated in >=5% of cell ID')+theme(legend.position = 'bottom',legend.text = (element_text(size=12)),panel.grid = element_line('grey90',linewidth = 0.1))+coord_flip()


specificity.index2=  x$data %>% 
  group_by(id) %>% 
  filter(pct.exp >0 & avg.exp.scaled >=0.1)%>% #& features.plot %in% goi1) 
  count(pct.exp >0) 

summary.medusa.exp=as.data.frame(levels(Aa.Alldata$ID.separate))
summary.medusa.exp$n = 0L
names(summary.medusa.exp)[1] = 'id'
summary.medusa.exp$n[summary.medusa.exp$id %in% specificity.index2$id]=specificity.index2$n
#Fig1E upper ----
barplot(summary.medusa.exp$n,names.arg = summary.medusa.exp$id,las=2,col=all.clusters.cp[summary.medusa.exp$id])+title('Filtered.medusa.genes, detection / cluster')

x=DotPlot(Aa.Alldata,assay='RNA',features=gene.plot[h$order],group.by = 'ID.separate')

medusa.specific.filtered=x$data %>% group_by(id) %>%  filter(pct.exp >5 & avg.exp.scaled >=0.1)
names(medusa.specific.filtered)[3]='gene_short_name'
medusa.specific.filtered.list=merge(medusa.specific.filtered,Aa_genes,all.x=T,all.y=F)

write.csv(medusa.specific.filtered.list,file='Ac_manuscript_final/DS1.1b.csv')


# Figure S1 ----
{
load (file='Ac_manuscript/topGO.medusa')
FS1A=topGO.medusa[[1]]+topGO.medusa[[2]]
FS1A
data1=SetIdent(Aa.Alldata,value='orig.ident')
cellID=NULL
for (i in 1:length(levels(data1)))
  cellID[[i]]=WhichCells(data1,idents=levels(data1)[i])
names(cellID)=levels(data1)


DimPlot(data1,cells.highlight = cellID[1],cols.highlight = LibCP[1])+
  DimPlot(data1,cells.highlight = cellID[2],cols.highlight = LibCP[2])+
  DimPlot(data1,cells.highlight = cellID[3],cols.highlight = LibCP[3])+
  DimPlot(data1,cells.highlight = cellID[4],cols.highlight = LibCP[4])+
  DimPlot(data1,cells.highlight = cellID[5],cols.highlight = LibCP[5])+
  DimPlot(data1,cells.highlight = cellID[6],cols.highlight = LibCP[6])+
  DimPlot(data1,cells.highlight = cellID[7],cols.highlight = LibCP[7])+
  DimPlot(data1,cells.highlight = cellID[8],cols.highlight = LibCP[8])+
  DimPlot(data1,cells.highlight = cellID[9],cols.highlight = LibCP[9])+
  DimPlot(data1,cells.highlight = cellID[10],cols.highlight = LibCP[10])+
  DimPlot(data1,cells.highlight = cellID[11],cols.highlight = LibCP[11])+
  DimPlot(data1,cells.highlight = cellID[12],cols.highlight = LibCP[12])+plot_layout(ncol=3)
}

# Figure 2 ----
## Fig.2a
{
  #PLOT WITHOUT THE EPITHELIA 
  clust.cp=all.clusters.cp
  Aa.Alldata<-SetIdent(Aa.Alldata,value='IDs')
  data1<-subset(Aa.Alldata,idents=levels(Aa.Alldata)[c(1,4:8)])
  data1$ID.separate<-droplevels(data1$ID.separate)
  
  data1=SetIdent(data1,value = 'ID.separate')
  levels(data1)
  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory.tissue))
  colnames(ids.cluster.library) = c('ID','Stage','CellCount')
  
  print(
    ggplot(ids.cluster.library, aes(fill=ID, y= CellCount,
                                    x=Stage)) +
      geom_bar(mapping =aes(fill=ID, y= (CellCount),
                            x=(Stage)),
               position="fill", stat="identity", width = 0.5)+
      scale_fill_manual(values = clust.cp)+
      theme(axis.text.x = element_text(#face="bold", color="#993333", 
        size=8, angle=-45,hjust=0,vjust = 0.5))+
      geom_area(mapping =aes(fill=ID, y= (CellCount),
                             x=as.integer(Stage)),
                position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
      geom_bar(mapping =aes(fill=ID, y= (CellCount),#this re-plots the bars over the area
                            x=(Stage)),
               position="fill", stat="identity", width = 0.5)+labs(title = 'All cluster identities')+theme(legend.position = 'none')
  )
}
# Figure 2b ----
#use PCs to build the tree:
embeddings <- Embeddings(object = Aa.Alldata, reduction = 'pca')[,1:50]
Aa.Alldata<-SetIdent(Aa.Alldata,value='ID.separate')
data.dims=matrix('0',50L,length(levels(Aa.Alldata)))
for (i in 1:length(levels(Aa.Alldata)))  
{  cells <- WhichCells(object = Aa.Alldata, idents = levels(Aa.Alldata)[i])
data.dims[,i] <- colMeans(embeddings[cells,])
}
colnames(x = data.dims) <- levels(x = Aa.Alldata)
data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
nj.tree <- ape::nj(data.dist)
ape::plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = F,use.edge.length = F,lab4ut = 'axial',tip.color = all.clusters.cp, cex = 0.8,no.margin = T)

## Fig.2C ----
#colour TFs by family
{
  load(file='Ac_manuscript_final/all.idents.markers.RData')
write.csv(all.markers.allidents,file='Ac_manuscript_final/DS1.4.csv')
write.csv(all.markers.TF.allidents,file='Ac_manuscript_final/DS1.4b.csv')
Aa.Alldata=SetIdent(Aa.Alldata,value='ID.separate')
all.markers.all.celltype.filter=all.markers.TF.allidents
all.markers.all.celltype.filter = all.markers.TF.allidents[all.markers.TF.allidents$pct.2<=0.5,]
list = NULL
for (c in (1:length(levels(Aa.Alldata@active.ident))))
{
  x=all.markers.all.celltype.filter[as.numeric(all.markers.all.celltype.filter$cluster)==c,][1:min(3,length(which(as.numeric(all.markers.all.celltype.filter$cluster)==c))),7]
  list=c(list,x)
}
list=unique(list)

DotPlot(SetIdent(Aa.Alldata,value='ID.separate'),'RNA', features = (list),scale.by='size' , col.min = 0, cols = c('lightgrey','darkred'),dot.min = 0.02)+
  RotatedAxis() +labs(title='Top 3 DETFs / cell state',subtitle = 'adj.p-val <= 0.001')+theme(legend.position = 'bottom',legend.title = element_text(size=8),legend.text = element_text(size=8),panel.grid = element_line(color='grey90',linewidth = 0.25))&coord_flip()

}
#colour TFs by family
{
list.DEG=unique(all.markers.all.celltype.filter$gene)
list.short=list #plotted above
list.short=as.data.frame(list.short)
# import TF family from annotations file
x=Aa_genes$TF_Fam[match(list.short$list,Aa_genes$gene_short_name)]
list.short$family=x

TF.list.filtered=AaTF_list[rowSums(Aa.Alldata@assays$RNA@counts[AaTF_list,])>=1]
list.ALL=TF.list.filtered #expressed TFs
#colour palette for TFs:
x=Aa_genes$TF_Fam[match(list.ALL,Aa_genes$gene_short_name)]
g.cp=clust.cp.separate[1:length(unique(x))]
names(g.cp)=unique(x)

list.ALL=as.data.frame(list.ALL)
list.ALL$family=x
#filter for at least 4 members
tfs=as.data.frame(table(list.ALL$family))
tfs=as.character(tfs$Var1[tfs$Freq>3])
x=Aa_genes$TF_Fam[match(list.ALL$list,Aa_genes$gene_short_name)]
x=x[tfs%in%x]

list=list.DEG#list.short
list=as.data.frame(list)
x=Aa_genes$TF_Fam[match(list$list,Aa_genes$gene_short_name)]
list$family=x
ind=sort(unique(list$family),index.return=T)
list.DEG=list

list.ratio=table(list.ALL$family)
list.ratio=as.data.frame(list.ratio)
row.names(list.ratio)=list.ratio$Var1
list.ratio$DEG=as.numeric(0)
colnames(list.ratio)=c('FAM','ALL','DEG')
x=as.data.frame(table(list.DEG$family))
rownames(x)=x$Var1
colnames(x)=c('FAM','DEG')
list.ratio$DEG[match(x$FAM,list.ratio$FAM)]=as.numeric(x$DEG)
list.ratio$ratio=as.numeric(list.ratio$DEG)/as.numeric(list.ratio$ALL)*100

list.ratio$fraction.DEG=list.ratio$DEG/sum(list.ratio$DEG)*100
list.ratio$fraction.ALL=list.ratio$ALL /sum(list.ratio$ALL)*100
list.ratio$dif=list.ratio$fraction.DEG-list.ratio$fraction.ALL
names(g.cp)=list.ratio$FAM

# family percent of all DEG - family percent of all TFs
ind=which(abs(list.ratio$dif)>=2.5) # Only consider plus or minus 5%
ind
barplot(list.ratio$dif[ind],col=g.cp[ind],horiz = T,xlim = c(-20,10))
pal.bands(g.cp[list.short$family]) #plotted for the figure
pal.bands(g.cp[unique(list.short$family)]) 
pal.bands(g.cp[unique(ind)]) #all TFs
}

# Figure SX --- all DEGs
{  
  load(file='Ac_manuscript/all.idents.markers.RData')
  Aa.Alldata=SetIdent(Aa.Alldata,value='ID.separate')
  list = NULL
  for (c in (1:length(levels(Aa.Alldata@active.ident))))
  {
    x=all.markers.allidents[as.numeric(all.markers.allidents$cluster)==5,][1:min(50,length(which(as.numeric(all.markers.allidents$cluster)==c))),7]
    list=c(list,x)
  }
  list=unique(list)
  # list=list[!is.na(list)]
  
  DotPlot(SetIdent(Aa.Alldata,value='ID.separate'),'RNA', features = (list),scale.by='size' , col.min = 0, cols = c('lightgrey','darkred'),dot.min = 0.02)+
    RotatedAxis() +labs(title='Top 3 DEGs / cell state')+theme(legend.position = 'bottom',legend.title = element_text(size=8),legend.text = element_text(size=8),panel.grid = element_line(color='grey90',linewidth = 0.25))&coord_flip()
  
}
# Figure S3 ----
{
load (file = 'Ac_manuscript/AllData.integrated.subsets.markers.RObj')
list = NULL
for (c in (1:length(levels(data1@active.ident))))
{
  x=markers$gland.dig[as.numeric(markers$gland.dig$cluster)==c,][1:min(10,length(which(as.numeric(markers$gland.dig$cluster)==c))),7]
  list=c(list,x)
}
list=unique(list)
list=list[!is.na(list)]

data1=data1.subsets$gland.dig
ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory.tissue))
colnames(ids.cluster.library) = c('ID','Stage','CellCount')
ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                x=Stage)) +
  geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                        x=(Stage)),
           position="fill", stat="identity", width = 0.5)+
  scale_fill_manual(values = clust.cp.separate)+
  theme(axis.text.x = element_text(#face="bold", color="#993333", 
    size=8, angle=-45,hjust=0,vjust = 0.5))+
  geom_area(mapping =aes(fill=ID, y= log(CellCount),
                         x=as.integer(Stage)),
            position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
  geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                        x=(Stage)),
           position="fill", stat="identity", width = 0.5)&
  DimPlot(data1.subsets$gland.dig,cols=clust.cp.separate)+NoLegend()+NoAxes()

  DotPlot(data1.subsets$gland.dig,'RNA', features = (list),scale.by='radius' , col.min = 0, split.by = 'lifehistory',cols = LibCP.stages,# c('lightgrey','darkred'),
          dot.min = 0.02)+
  RotatedAxis() +FontSize(6,8) +labs(title='Top 3 DETFs / cell sub-type',subtitle = 'adj.p-val <= 0.001')+theme(legend.position = 'bottom',legend.title = element_text(size=8),legend.text = element_text(size=8),panel.grid = element_line(color='grey90',linewidth = 0.1))+coord_flip()

  # Figure S4 ----
    data1=data1.subsets$gland.muc
  list = NULL
  for (c in (1:length(levels(data1@active.ident))))
  {
    x=markers$gland.muc[as.numeric(markers$gland.muc$cluster)==c,][1:min(3,length(which(as.numeric(markers$gland.muc$cluster)==c))),7]
    list=c(list,x)
  }
  list=unique(list)
  list=list[!is.na(list)]
  x=markers$gland.muc[(markers$gland.muc$cluster)=='mu.unchar.medusa',7]

  ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory.tissue))
  colnames(ids.cluster.library) = c('ID','Stage','CellCount')
  ggplot(ids.cluster.library, aes(fill=ID, y= log(CellCount),
                                  x=Stage)) +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),
                          x=(Stage)),
             position="fill", stat="identity", width = 0.5)+
    scale_fill_manual(values = clust.cp.separate)+
    theme(axis.text.x = element_text(#face="bold", color="#993333", 
      size=8, angle=-45,hjust=0,vjust = 0.5))+
    geom_area(mapping =aes(fill=ID, y= log(CellCount),
                           x=as.integer(Stage)),
              position="fill", stat="identity",alpha=0.2 , size=.5, colour="white") +
    geom_bar(mapping =aes(fill=ID, y= log(CellCount),#this re-plots the bars over the area
                          x=(Stage)),
             position="fill", stat="identity", width = 0.5)&
    DimPlot(data1.subsets$gland.muc,cols=clust.cp.separate)+NoLegend()+NoAxes()
  
  DotPlot(data1.subsets$gland.muc,'RNA', features = (list),scale.by='radius' , col.min = 0, split.by = 'lifehistory',cols = LibCP.stages,# c('lightgrey','darkred'),
          dot.min = 0.02)+
    RotatedAxis() +FontSize(6,8) +labs(title='Top 3 DEGs / cell sub-type',subtitle = 'adj.p-val <= 0.001')+theme(legend.position = 'bottom',legend.title = element_text(size=8),legend.text = element_text(size=8),panel.grid = element_line(color='grey90',linewidth = 0.1))+coord_flip()
}

# Figure S5 ----
  load(file='Ac_manuscript/all.idents.markers.RData')
  data1=SetIdent(Aa.Alldata,value = 'ID.separate')
data1@active.assay='RNA'

list = NULL
for (c in  1:length(levels(data1@active.ident)))
{
x=all.markers.allidents[as.numeric(all.markers.allidents$cluster)==c,][1:min(4,length(which(as.numeric(all.markers.allidents$cluster)==c))),7]
list=c(list,x)
}
  
  DotPlot(data1,'RNA', features = unique(list),
          scale.by='radius' , col.min = 0, col.max = 4,  
          cols = c('lightgrey','darkred'))+
    RotatedAxis() +FontSize(6,8) +labs(title='All Cell States: Top 4 DEGs')+theme(legend.position = 'bottom',legend.title = element_text(size=6),panel.grid = element_line(color='grey90',linewidth = 0.1))+coord_flip()
  
# Figure 3 ----
  ## Fig3A,b ----
clust.cp=c("#FB6496",brewer.rdpu(20)[c(8:11,20,17,12,16)])
clust.cp[2]='lightpink'
clust.cp[3]='green4'
clust.cp[4]='slateblue'
clust.cp[5]='slateblue4'
clust.cp[8]='black'
clust.cp[9]='orange'
names(clust.cp)=levels(data1.subsets$cnidocyte)
cnido.clust=DimPlot(data1.subsets$cnidocyte,cols=clust.cp)+NoAxes()+NoLegend()
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
DimPlot(Aa.Alldata,cells.highlight = coi,cols.highlight = rev(clust.cp[sort(names(coi))]))&NoAxes()&NoLegend()

##Fig 3C ----
load('OMA.genes.RData')
{
  load(file='cnido.nv.ac.RData')
  temp.n=cnido.nv
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
  
  goi=Aa_genes$gene_short_name[match(OMA$Aur2,Aa_genes$geneID)]
  temp=cnido.ac
  temp@active.assay='RNA'
  rownames(temp@assays$RNA@counts)=Aa_genes$gene_short_name
  rownames(temp@assays$RNA@data)=Aa_genes$gene_short_name
  rownames(temp@assays$RNA@meta.features)=Aa_genes$gene_short_name # otherwise fails with vs5
  temp@assays$RNA@var.features=Aa_genes$gene_short_name # otherwise fails with vs5
  temp$species='aurelia'
  temp$ID.species=paste('Ac',temp$IDs,sep='.')
  temp<-SetIdent(temp,value = 'ID.species')
  temp<-ScaleData(temp,features = goi,split.by = 'orig.ident')
  levels(temp)
  
  all.cnido=merge(temp.n,temp) #good to go... run the pipeline
  dim(all.cnido)
  # use only the OMA:
  data1=subset(all.cnido,features = OMA$name)
  #drop the planula and gastrula clusters:
  data1<-subset(data1,idents=levels(data1)[c(3:5,7:9,12:16,20:21)])
  #filter out really poor-quality cells that no longer have enough information:
  data1 <- subset(x = data1, subset = nFeature_RNA > 350 & nCount_RNA > 500)
  data1 <- NormalizeData(data1, scale.factor = 5000)
    data1 <- FindVariableFeatures(data1, nfeatures = 5000)
    data1 <- ScaleData(data1,split.by = 'orig.ident')
  data1 <- RunPCA(data1, npcs = 50)

  
  ### harmony ----
  data1<-harmony::RunHarmony(data1,group.by.vars = 'orig.ident')
  embeddings <- Embeddings(object = data1, reduction = 'harmony')[,1:20]
  data.dims=matrix('0',20L,length(levels(data1)))
  for (i in 1:length(levels(data1)))  
  {  cells <- WhichCells(object = data1, idents = levels(data1)[i])
  data.dims[,i] <- colMeans(embeddings[cells,])
  }
  # library(ape)
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
  ape::plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = T,use.edge.length = T,lab4ut = 'axial', cex = 1,no.margin = T)
  for (i in 1:length(node_col)) {
    nodelabels(node = length(nj.tree$tip.label)+i, pch=21, col="black", bg=node_col[i], cex=2)
  }
  
}
## Fig3 E ----
library("Nebulosa")
temp@active.assay='RNA'
goi.a=c("myc4","SoxC", "prd13-like-1", "zn431-like-1", 'pax2a-like-1', "jun-like3",'FOS-OMA-1', "gfi1b-like-1", "FoxN1", 'sx19b-like-2','myc.1','myc5')#

  plot_density(temp, goi.a,reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')

  plot_density(temp, 'myc3',reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')

  goi.n=c('myc4','SoxC','Prdm13','ZNF845','PAX2A-like-1','JUN-like-11','FOS-like-4','GFI1B-like-1','FoxN1','SoxA','myc','myc5')
  goi.n=genes$name[match(goi.n,genes$gene_short_name)]
  plot_density(temp.n, goi.n,reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')
  plot_density(temp.n, 'Myc3 : NV2.18825',reduction = 'umap')&NoAxes()&NoLegend()&scale_colour_gradient(low = 'grey90',high = 'black',na.value = 'grey90')


  

## Fig4 ----
#all.neuron
{
load(file='Ac_manuscript/all.neurons.2.Robj')
  load(file = 'Ac_manuscript/AllData.integrated.subsets.markers.RObj')
data1<-all.neurons  
  data1$orig.ident=as.factor(data1$orig.ident)
  data1$orig.ident=factor(data1$orig.ident,levels(data1$orig.ident)[c(8:10,12,11,1:7)])
  data1$lifehistory.tissue=as.factor(data1$lifehistory.tissue)
  data1$lifehistory.tissue=factor(data1$lifehistory.tissue,levels(data1$lifehistory.tissue)[c(7,8,6,1:5)])
  DimPlot(data1,cols=neur.cp)&NoAxes()&NoLegend()
    ids.cluster.library = as.data.frame(table(Idents(data1), data1@meta.data$lifehistory.tissue))
    colnames(ids.cluster.library) = c('ID','Library','CellCount')
    ids.cluster.library$ID<-factor(ids.cluster.library$ID,levels(ids.cluster.library$ID)[rev(1:20)])
    
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
  
    coi=NULL
    for(c in 1:length(levels(all.neurons)))
      coi[[c]]=WhichCells(all.neurons,idents=levels(all.neurons)[c])
    names(coi)=levels(all.neurons)
    DimPlot(Aa.Alldata,cells.highlight = coi,cols.highlight = rev(neur.cp[sort(names(coi))]))&NoAxes()&NoLegend()
    all.neurons$IDs=as.factor(all.neurons$IDs)
    all.neurons$IDs=factor(all.neurons$IDs,levels(all.neurons$IDs)[rev(c(5,1:4,6,17:20,7:9,14:16,10:13))])
    DotPlot(data1,'RNA',c('ins1a-like-1'),cols=c('lightgrey','darkorange3'),scale.by = 'size',dot.min = 0.01,col.min = 0,scale = F)&RotatedAxis()&
    DotPlot(data1,'RNA',c('ashA','pou4-like1'),cols=c('lightgrey','darkgreen'),scale.by = 'radius',dot.min = 0.01,col.min = 0,scale = F)&RotatedAxis()&theme(legend.title = element_text(size = 8),legend.text = element_text(size=8),panel.grid = element_line(linewidth = 0.2,colour = 'grey90'))#
    embeddings <- Embeddings(object = data1, reduction = 'pca')[,1:50]
    data.dims=matrix('0',50L,length(levels(data1)))
    for (i in 1:length(levels(data1)))  
    {  cells <- WhichCells(object = data1, idents = levels(data1)[i])
    data.dims[,i] <- colMeans(embeddings[cells,])
    }
    colnames(x = data.dims) <- levels(x = data1)
    data.dist <- dist(x = t(x = data.dims),method = 'euclidean')
    nj.tree <- ape::nj(data.dist)
    plot.phylo(nj.tree, type = "unrooted", label.offset = 0.5,show.node.label = F,use.edge.length = T,lab4ut = 'axial',tip.color = clust.cp, cex = 0.8,no.margin = T)
    
    # 4d----
    data1@active.assay = 'RNA'
    all.markers=rbind(markers$neural.1,markers$neural.2)
list = NULL
for (i in 1:length(levels(data1@active.ident)))
{
  x = all.markers[all.markers$cluster == levels(data1@active.ident)[i], ][1:min(3, length(which(as.numeric(all.markers$cluster) == i))), 7]
  
  list = c(list, x)
}
DotPlot(all.neurons,'RNA', tbb,scale.by='radius' ,dot.min = 0.1, col.min = 0,cols = c('lightgrey','darkred'))+RotatedAxis() + labs(title='All neurons',subtitle = 'DEGs p.val <=0.001')+theme(legend.position = 'bottom',legend.text = (element_text(size=8)),panel.grid = element_line('grey90',linewidth = 0.1))

goi=c('tba1-like3','tbb-like-1','prfa-like-3')#,'tbb-like-2','tbb-like-7')#'prfa-like-1','prfa-like-2','prfa4-like'
tbb=c('tbb-like-1','tbb-like-2','tbb-like-3','tbb-like-4','tbb-like-5','tbb-like-6','tbb-like-7','tbb-like-8','tbb-like-9','tbb-like-10','tbb-like-11')
all.neurons@active.assay='RNA'
FeaturePlot(all.neurons,goi,order=T, ncol=5)&NoAxes()&NoLegend()&scale_colour_gradientn(colors=gene.cp)
plot_density(all.neurons, goi[1],reduction = 'umap')+
  scale_colour_gradient(low = 'grey80',high = 'red',na.value = 'grey90')+plot_density(all.neurons, goi[2],reduction = 'umap')+
  scale_colour_gradient(low = 'grey80',high = 'red',na.value = 'grey90')+plot_density(all.neurons, goi[3],reduction = 'umap')+
  scale_colour_gradient(low = 'grey80',high = 'aquamarine4',na.value = 'grey90')+plot_layout(ncol=3)&NoAxes()

}

# figS6----
atoh8=Aa_genes$gene_short_name[match(c('AAUR2.25654','AAUR2.19288','AAUR2.45802','AAUR2.21664','AAUR2.42412','AAUR2.10736'),Aa_genes$geneID)] #not expressed:,'AAUR2.45801','AAUR2.7216','AAUR2.15768','AAUR2.17097'
DotPlot(all.neurons,'RNA', features=(atoh8),scale.by='size' ,dot.min = 0.001, idents=levels(all.neurons)[1:13],col.min = 0,cols = c('lightgrey','steelblue'))+RotatedAxis()+theme(legend.position = 'bottom',legend.text = (element_text(size=8)),panel.grid = element_line('grey90',linewidth = 0.1))+coord_flip()

load(file='nem.neur.Robj')
nem.neur.reduced=subset(nem.neur,idents=levels(nem.neur)[c(2:4,6:10,13:21,23,39:46)])
atoh8.nem=genes$gene_short_name[match(c('NV2.11608','NV2.6611','NV2.6610','NV2.13885'),genes$geneID)] #not expressed:,'AAUR2.7216''AAUR2.15768','AAUR2.45801','AAUR2.45802','AAUR2.17097',
DotPlot(nem.neur.reduced,'RNA', features=atoh8.nem,scale.by='size' ,dot.min = 0.01,group.by = 'IDs', col.min = 0,cols = c('lightgrey','red'))+RotatedAxis() +FontSize(8,10)+ labs(title='Nematostella all neurons',subtitle = 'atoh8-like')+theme(legend.position = 'bottom',legend.text = (element_text(size=8)),panel.grid = element_line('grey90',linewidth = 0.1))+coord_flip()
                                     
#Figure S6
DotPlot(Aa.Alldata,'RNA', atoh8, group.by = 'ID.separate',scale.by='size' ,dot.min = 0.001, col.min = 0,cols = c('lightgrey','darkred'))+RotatedAxis() + labs(title='Atoh8 expression across all cell type states',subtitle = 'DEGs p.val <=0.001')+ theme(legend.position = 'bottom',legend.text = (element_text(size=8)),panel.grid = element_line('grey90',linewidth = 0.5))

#Figure S7 ----
# "the bulk of the upregulated smooth muscle shared with outer" -- plot all smooth muscle DEG on 'reduced' dataset.
load(file='Ac_manuscript/all.idents.markers.RData')

DotPlot(data1,'RNA', all.markers.allidents$gene[all.markers.allidents$cluster=='outer.smooth.muscle'], group.by = 'reduced',scale.by='size' ,dot.min = 0.001, col.min = 0,cols = c('lightgrey','darkred'))+RotatedAxis() + labs(title='Atoh8 expression across all cell type states',subtitle = 'DEGs p.val <=0.001')+ theme(legend.position = 'bottom',legend.text = (element_text(size=8)),panel.grid = element_line('grey90',linewidth = 0.5))


# Figure 8 ----
Aa.Alldata<-SetIdent(Aa.Alldata,value='IDs')
coi.st=WhichCells(Aa.Alldata,idents ='muscle.st')

Aa.Alldata<-SetIdent(Aa.Alldata,value='ID.separate')
coi.sm=WhichCells(Aa.Alldata,idents ='outer.smooth.muscle')

data1=Aa.Alldata
Idents(data1) <- data1$IDs
data1$reduced=as.character(data1@active.ident)
data1$reduced[coi.sm]=as.character(Aa.Alldata$ID.separate[coi.sm])
data1<-SetIdent(data1,value='reduced')

#further reduce:
levels(data1$reduced)
new.names=c('muscle.st','muscle.sm','non-muscle','non-muscle','non-muscle','non-muscle','non-muscle','non-muscle','non-muscle')
data1$reduced.muscle=data1$reduced
levels(data1$reduced.muscle)=new.names

#generate gene lists:
calm=Aa_genes$gene_short_name[grep('calm',Aa_genes$gene_short_name)]
actin=Aa_genes$gene_short_name[grep('act',Aa_genes$gene_short_name)]
myosin=Aa_genes$gene_short_name[grep('myosin',Aa_genes$best_OG_desc)]
melc=Aa_genes$gene_short_name[grep('MYL6B',Aa_genes$Preferred_name)]
melc=c(melc,'melc12-like')
myosin=c(myosin,'stmyhc1','mrlc2','mrlc1-like','mypt1-like-1','mypt2-like-2','mypt2-like-4','mypt2-like-5')
x=DotPlot(data1,features=c(ach,glu,myosin,melc,actin,calm),group.by = 'reduced')
x=x$data
x=x[!is.nan(x$avg.exp.scaled),]
#subset the muscle:
x2a=x[x$id == 'muscle.st' | x$id == 'outer.smooth.muscle',]
#filter for expression min. cells
x2=x2a[x2a$pct.exp>=10,]
#filter genes for up-regulation in muscle
goi=x2$features.plot[x2$avg.exp.scaled>0]
goi=unique(as.character(goi))
#filter gene sets:
actin.f=actin[actin %in% goi]
actin.f=actin.f[c(2,1,4,5)] #>5 c(4,3,8,9)
ach.f=ach[ach %in% goi]
glu.f=glu[glu %in% goi]
calm.f=calm[calm %in% goi]

melc.f=melc[melc %in% goi]
melc.f=sort(melc.f)

myosin.f=myosin[myosin %in% goi]
myosin.f=myosin.f[sort(Aa_genes$best_OG_desc[match(myosin.f, Aa_genes$gene_short_name)],index.return=T)$ix]
myosin.ann=Aa_genes[match(myosin.f, Aa_genes$gene_short_name),]
myosin.f
#remove 'other' myosin categories
myosin.f=myosin.f[c(4:10,15,16:20,24:25,13,21,14,26,11,12,27,23)]# filter 5 c(7:15,23,35,34,21,30,24:29,18,19,37,22,36,33)
#order the gene lists for plotting:
goi=c(actin.f,myosin.f[1:13],melc.f,myosin.f[14:19],'smtl1-like',myosin.f[20:22],calm.f)

#Plot with names
p=DotPlot(data1,features=rev(goi),scale.by='size' ,dot.min = 0.01,group.by = 'reduced.muscle', col.min = -0.5,cols = c('lightgrey','red'))+RotatedAxis() +FontSize(16,16)+ labs(title='Ortho-Group description')+theme(legend.position = 'bottom',legend.text = (element_text(size=12)),panel.grid = element_line('black',linewidth = 0.1))+coord_flip()

#update names to gene models:
gene.names=p$data$features.plot
for (i in 1:length(levels(gene.names)))
  levels(gene.names)[i]=Aa_genes$geneID[match(as.character(levels(gene.names)[i]), Aa_genes$gene_short_name)]
p2=p
p2$data$features.plot=gene.names
p+p2 #check this worked correctly:

# plot Figure8:
p2

# Generate CellBrowser Files ----
Aa.Alldata@active.assay='RNA'
AllData <- DietSeurat(Aa.Alldata,data=T,assays='RNA',dimreducs = c('umap.int'))
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
  # jpeg(file=paste('TEMP_CellBrowser/',my.list[i],'.jpg',sep=''))
  print(DimPlot(Aa.Alldata,cells.highlight = colnames(data1.subsets[[i]]),cols.highlight = 'red',reduction='umap.h')&NoAxes()&NoLegend()&labs(title=my.list[i]))
  # dev.off()
}

jpeg(file=paste('TEMP_CellBrowser/SeriesOverview.JPEG'))
print(DimPlot(AllData,cols=LibCP,group.by = 'orig.ident',pt.size = 0.5)&NoAxes()&labs(title = 'All Data'))
dev.off()

write.csv(Aa_genes[,3:4],'TEMP_CellBrowser/gene.table.csv')
