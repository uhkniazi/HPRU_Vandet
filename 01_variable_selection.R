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

################ variable selection steps
set.seed(123)

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)

## random forest step
oVar.r = CVariableSelection.RandomForest(dfData[-test,], groups = fGroups[-test], boot.num=100, big.warn = F)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
cvTopGenes = rownames(dfRF)[1:30]

## subset selection
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.s = CVariableSelection.ReduceModel(dfData[-test,], fGroups[-test], boot.num=100)
plot.var.selection(oVar.s)

## 10 fold cv
## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
for (i in 1:8){
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


