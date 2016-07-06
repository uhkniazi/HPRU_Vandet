# Name: 02_de_analysis.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 04/07/2016
# Desc: blk1_12 dataset 

p.old = par()

### import the data
dfData = read.csv('Data_external/BLK1_12_Umar.csv', header=T, stringsAsFactors=F)

# create short names for proteins
rn = dfData$Protein
rn = gsub('([\\w\\d]+);[\\w\\d]*', replacement = '\\1', rn, ignore.case = T, perl = T)
rownames(dfData) = rn

## load the annotation file and create grouping factor
dfAnno = read.csv('Data_external/BLK1_12_Umar_annotation.csv', header=T)
fDos = dfAnno$Dosanjh
table(fDos)
fGroups = rep(NA, length.out = length(fDos))
fGroups[fDos == '1'] = 'D1.2'
fGroups[fDos == '2'] = 'D1.2'
fGroups[!(fDos %in% c('1', '2'))] = 'D4'
fGroups = factor(fGroups, levels = c('D4', 'D1.2'))
table(fGroups)
dfAnno$fGroups = fGroups

## remove people with lymph node involvement
i = dfAnno$Lymph.node.involvement == 'Yes'
table(i)
dfAnno = dfAnno[!i,]

## select matching samples between the two data frame
f = colnames(dfData) %in% as.character(dfAnno$IDEA.ID)
table(f)
dfData = dfData[,f]

i = match(colnames(dfData), as.character(dfAnno$IDEA.ID))
dfAnno = dfAnno[i,]

## create factors for ethnicity and age
iAge = dfAnno$Age
fEthnicity = gsub('(\\w+):.+', '\\1', dfAnno$Ehtnicity.code.1, perl = T)
i = which(fEthnicity %in% c('Unable/unwilling to respond', 'Mixed'))
fEthnicity[i] = 'Other'
fEthnicity = factor(fEthnicity)
summary(iAge)
summary(fEthnicity)
dfAnno$fEthnicity = fEthnicity

## format data for DE analysis
mData = as.matrix(dfData)
mData = t(mData)
fGroups = dfAnno$fGroups
table(fGroups)

## anova to find DE
p.anova = apply(mData, 2, function(x){
  fit = lm(x ~ fGroups)
  return(anova(fit)$Pr[1])
})

hist(p.anova)
p.anova.adj = p.adjust(p.anova, method='BH')
hist(p.anova.adj)

table(p.anova.adj < 0.1)
table(p.anova.adj < 0.01)

i = which(p.anova.adj < 0.1)

## fit second model with age and ethnicity
mData = mData[,i]

p.anova = apply(mData, 2, function(x){
  fit = lm(x ~ fGroups + fEthnicity + iAge)
  return(anova(fit)$Pr[1])
})

hist(p.anova)
p.anova.adj = p.adjust(p.anova, method='BH')
hist(p.anova.adj)

table(p.anova.adj < 0.1)
table(p.anova.adj < 0.01)

i = which(p.anova.adj < 0.01)

cvTopProteins = names(i)
dfData.sub = data.frame(t(na.omit(t(mData[,i]))))
dim(dfData.sub)
head(dfData.sub)

################ variable selection steps
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfData = dfData.sub

set.seed(123)

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)
table(fGroups[test]); table(fGroups[-test])
## random forest step
oVar.r = CVariableSelection.RandomForest(dfData[-test,], groups = fGroups[-test], boot.num=100, big.warn = F)
plot.var.selection(oVar.r)
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
cvTopGenes = rownames(dfRF)[1:30]

## subset selection
dfData = dfData[,colnames(dfData) %in% cvTopGenes]
oVar.s = CVariableSelection.ReduceModel(dfData[-test,], fGroups[-test], boot.num=100)
plot.var.selection(oVar.s)
dfPrint = NULL
par(mfrow=c(1,2))
## 10 fold cv
## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
for (i in 1:6){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.s, i)
  dfData.train = as.data.frame(dfData[-test, cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  dfData.test = data.frame(dfData[test, cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fGroups[test],
                             train.groups = fGroups[-test], level.predict = 'D1.2', boot.num = 500)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  x = getAUCVector(oCV)
  vc = (paste('Variable Count', i))
  gn = paste(cvTopGenes.sub, collapse = ' ')
  sig = (signif(quantile(x, probs = c(0.025, 0.975)), 2))
  dfPrint = rbind(dfPrint, c(vc, gn, sig))
}

dfPrint

######### save list of DE genes and p-values
dfExport = data.frame(Proteins=names(p.anova), p.value=p.anova, p.adjusted=p.anova.adj)
write.csv(dfExport, file='Temp/proteins.csv')
