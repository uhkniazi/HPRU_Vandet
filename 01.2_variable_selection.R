# Name: 01_variable_selection.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 08/06/2016
# Desc: data import and variable selection

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

url = 'https://raw.githubusercontent.com/uhkniazi/CGraphClust/master/CGraphClust.R'
download(url, 'CGraphClust.R')

# load the required packages
source('CGraphClust.R')
# delete the file after source
unlink('CGraphClust.R')

p.old = par()
### import the data
dfData = read.csv('Data_external/data.csv', header=T, stringsAsFactors=F)
dfData = na.omit(dfData)

# create short names for proteins
rn = dfData$Majority.protein.IDs
rn = gsub('([\\w\\d]+);[\\w\\d]*', replacement = '\\1', rn, ignore.case = T, perl = T)
rownames(dfData) = rn

## load the annotation file and create grouping factor
dfAnno = read.csv('Data_external/annotation.csv', header=T)
fDos = dfAnno$Dosanjh
table(fDos)
fGroups = rep(NA, length.out = length(fDos))
fGroups[fDos == '1'] = 'D1.2'
fGroups[fDos == '2'] = 'D1.2'
fGroups[!(fDos %in% c('1', '2'))] = 'D3.4'
fGroups = factor(fGroups, levels = c('D3.4', 'D1.2'))
table(fGroups)
dfAnno$fGroups = fGroups

## select matching samples between the two data frame
f = colnames(dfData) %in% as.character(dfAnno$Patient.Study.ID)
dfData = dfData[,f]

i = match(colnames(dfData), as.character(dfAnno$Patient.Study.ID))
dfAnno = dfAnno[i,]

## format data for variable selection
mData = as.matrix(dfData)
mData = t(mData)
dfData = data.frame(mData)
fGroups = dfAnno$fGroups

### clustering using reactome data 
mCounts = mData
library(org.Hs.eg.db)

# convert enterez ids to uniprot as Reactome database file uses UNIPROT ids
dfMap = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'ENTREZID', 'UNIPROT')
dfMap = na.omit(dfMap)

### load the uniprot2reactome mapping obtained from
# http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x

## map reactome ids to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('V1', 'V2')]
dfGraph = na.omit(dfGraph)

# select genes that have a reactome term attached
n = unique(dfGraph$V1)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
# order the count matrix based on grouping factor
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

# create a correlation matrix to decide cor cutoff
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# stabalize the data and check correlation again i.e. calculated bayes adjusted posterior means
mCounts.st = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(mCounts.st) = fGroups

# create a correlation matrix
mCor = cor(mCounts.st)
# check distribution 
hist(mCor, prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation to cluster positively and negatively correlated genes
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6, bSuppressPlots = F)

## general graph structure
## we would like to see how does the graph look like, are the clusters connected or in subgraphs
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))

## community structure
## overview of how the commuinties look like
# plot the main communities in 2 different ways
ig = getFinalGraph(oGr)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
# remove communities smaller than 5 members or choose a size of your liking
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

## if we want to look at the expression profiles of the top genes
# plot a heatmap of these top genes
library(NMF)
m1 = mCounts[,names(V(getFinalGraph(oGr)))]
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

m1 = mCounts[,names(V(getFinalGraph(oGr)))]
m1 = apply(m1, 2, f_ivStabilizeData, fGroups)
rownames(m1) = fGroups
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

cvTopGenes.cl = names(V(getFinalGraph(oGr)))

################ variable selection steps
dfData = data.frame(mData)
fGroups = dfAnno$fGroups
dfData = dfData[,colnames(dfData) %in% cvTopGenes.cl]

set.seed(123)

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)

## random forest step
oVar.r = CVariableSelection.RandomForest(dfData[-test,], groups = fGroups[-test], boot.num=100, big.warn = F)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
cvTopGenes = rownames(dfRF)[1:10]

## subset selection
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.s = CVariableSelection.ReduceModel(dfData[-test,], fGroups[-test], boot.num=100)
plot.var.selection(oVar.s)

## 10 fold cv
## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
for (i in 1:4){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.s, i)
  dfData.train = as.data.frame(dfData[-test, cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  dfData.test = data.frame(dfData[test, cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fGroups[test],
                             train.groups = fGroups[-test], level.predict = 'D1.2', boot.num = 100)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}


