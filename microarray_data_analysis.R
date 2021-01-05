

#Task 2 Importing Microarray Data to R

#Install the package affy, which is needed for the task
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("Rcpp")
library(Rcpp)
library(Biobase)
library(affy)
#Set the directory where the CEL files are saved.
setwd("E:/course/Advanced bioinformatic tools/E5/New folder/raw")

#Import the 24 data files.
data<- ReadAffy()
#Check the structure of the data.
str(data)
#Check the probe names and the number of probes. 
data.probenames<-probeNames(data)
length(data.probenames)
[1] 496468
#Check the signal intensities of different samples for PM and MM
data.mm<-mm(data)
data.pm<-pm(data)
#Check the number of PM(perfect match) and MM(mismatch)
nrow(data.pm)
[1] 496468
nrow(data.mm)
[1] 496468
data.mm[1:5,]

summary(data.pm)

#Task 3 Technical Quality Control 


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")
library(preprocessCore)
#Obtain the data matrix for pm 
data.pm<-pm(data)
#Compute logarithm of the intensity values with the base of 2.
data.pm.log2<-log(data.pm, 2)
#Normalize the intensity
data.pm.log2.normalize<-normalize.quantiles(data.pm.log2)
data.pm.log2[1:5,]
data.pm.log2.normalize[1:5,]
#Draw and compare probe intensity boxplots of each array before and after quantile normalization. 
boxplot(data.pm.log2, main="before normalization")
boxplot(data.pm.log2.normalize, main="after normalization")

boxplot(data, main="t")
boxplot(data.pm, main="pm")
data.pm.normalize<-normalize.quantiles(data.pm)

boxplot(data.pm.normalize, main="pm after n")

rma<-rma(data)

#Task 4 Data Preprocessing 
#Perform gcRMA procedure for all the arrays using the un-normalized raw data
#Apply function gcrma to the Affybatch data file 
source("https://bioconductor.org/biocLite.R")
biocLite("gcrma")
library(gcrma)
data.gcRMA<-gcrma(data)
str(data.gcRMA)
#extract expression matrix 
expression<-exprs(data.gcRMA)
str(data.gcRMA)
str(expression)
class(expression)
nrow(expression)

#Compute the variance for each probe set over all arrays. 
#Using mouse4302.db available in Bioconductor, annotate each probe set with corresponding gene symbol.
#Filter the probe sets with low variation
#Select one probe set for each gene symbol

#Install the packages needed for the 4 tasks. The annotation package mouse4302.db is for annotating the 
#probe sets. The package dplyr is for manipulating matrix in this case. The package genefilter is for probe 
#sets or gene filtration.  
library(dplyr)
BiocManager::install("genefilter")
source("https://bioconductor.org/biocLite.R")
biocLite("mouse4302.db")
library(mouse4302.db)
library(genefilter)
#Inspect the content of the annotation package 
ls("package:mouse4302.db")
mouse4302()
featureData(data.gcRMA)
#The following pocesses were implemented as described on 
#http://biolearnr.blogspot.com/2017/05/bfx-clinic-getting-up-to-date.html for importing 
#the annotation information into the ExpressionSet. 
#Some probe sets can have multiple gene annotations, which will results in a matrix with 
#more number of rows than the number of the probe set. However, we want to have one row for each 
#probe set that is mentioned in the ExpressionSet.  Thus, a function was created for collapsing 
#purpose (details can be seen from the link above).  
collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
  }

#Create the annotation matrix 
annots<-AnnotationDbi::select(x=mouse4302.db, keys=rownames(data.gcRMA), columns=c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"), keytype="PROBEID")%>%
group_by(PROBEID) %>%
summarise_each(funs(collapser)) %>%
ungroup
#Inspect the annotation data frame
class(annots)
annots[1:5,]
nrow(annots)
#Import the annotation information into the ExpressionSet
fd <- new("AnnotatedDataFrame",
  data = data.frame(annots[, -1], stringsAsFactors = FALSE)
  )
rownames(fd) <- annots$PROBEID
featureData(data.gcRMA) <- fd
#Inspect the annotation from the ExpressionSet
feature.data<-featureData(data.gcRMA)

str(feature.data)
feature.data@data[1:5,]
nrow(feature.data@data)
#Filter the probe sets using the function nsFilter. The probe sets with low variation are removed  
data.gcRMA.filter<-nsFilter(data.gcRMA)
#Only 10302 features are left after filtration 
data.gcRMA.filter

#Task5 Biological Quality Control 

#Compute the correlation between all the arrays. Plot the correlations as a heatmap to compare the 
#arrays. 

#Inspect the structure of the filtered ExpressionSet 
str(data.gcRMA.filter)
#Extract the expression level value for each sample 
expression.filter<-exprs(data.gcRMA.filter$eset)
#Inspect the expression matrix
class(expression.filter)
expression.filter[10:20,]
#Construct a correlation matrix
expression.filter.cor <- cor(expression.filter)
install.packages("corrplot")
library(corrplot)
#Plot the heat map 
corrplot(expression.filter.cor, is.corr=FALSE)

#Hierarchical clustering can be used to study the relative similarities and differences between 
#samples. Plot a dendrogram of the samples with appropriate distance measure and linkage method. 

#Transpose the expression matrix for subsequent processes
expression.filter.t<-t(expression.filter)
#Calculate the distance matrix r 
dd<-dist(scale(expression.filter.t), method = "euclidean")
#Implement hierarchical analysis 
hc<-hclust(dd,method = "ward.D2")
plot(hc, main = "Cluster dendrogram")

#Task 6  Identification of Differentially Expressed Genes 

#Compute expression fold changes for each gene between mice fed with low fat and high fat diet. 

#A character vector containing the group name (low fat or high fat diet) is created.
group=c("low","low","low","low","low","high","high","high","high","high","high","low","low","low","low","low","low","high","high","high","high","high","high","low")
#A model matrix is created based on the character vector created before. By defining the model matrix with ~ 0, limma will simply calculate the mean for each variable in each group later on.
modelmatirx=model.matrix(~ 0 + group)
#Extract the gene names and create a verctor for them. 
genelist<-data.gcRMA.filter$eset@featureData$SYMBOL
length(genelist)
#Add the the gene names as row names in the expression matrix 
rownames(expression.filter)<-genelist
expression.filter[1:10,]
# Limma package can be used for microarray data analysis using linear models. 
library(limma)
#This function fits linear model for the array. In this case, the mean of the expression of each gene in different groups is calculated. 
expression.filter.fit = lmFit(expression.filter,modelmatirx)
expression.filter.fit$coefficients[1:10,]
#A contrast matrix between the groups is made. The low fat diet group is specified as the baseline. 
contrast.matrix = makeContrasts(grouphigh-grouplow,levels=c("grouphigh","grouplow"))
# The expression between the two groups are compared  
expression.filter.fit.con = contrasts.fit(expression.filter.fit,contrast.matrix)
expression.filter.fit.con$coefficients[1:10,]

#Perform t-test for each gene between mice fed with low fat and high fat diet. 

# The moderated t-test is performed by using the eBayes
expression.filter.fit.con.eb = eBayes(expression.filter.fit.con) # The moderated t-test is performed by using the eBayes


# Using a volcano plot, determine appropriate thresholds of significance for the t-test and fold 
#change. Hint: Ideally, we want to continue with only few hundred of the most differentially expressed 
#genes. 

#Inspect the fold change and p value
expression.filter.fit.con$coefficients[1:5,]
expression.filter.fit.con.eb$p.value[1:5,]
#Combine the fold change and value.
data.combined = cbind(expression.filter.fit.con$coefficients,expression.filter.fit.con.eb$p.value)
class(data.combined)
BiocManager::install("a4Base")
library(a4Base)
#x is a vector for fold changes, y is a vector for p values, and z is a vector for gene names 
x<-data.combined[,1]
y<-data.combined[,2]
z<-rownames(data.combined)
volcanoPlot(x,y,z)


#Based on statistical and biological significance thresholds, obtain a list of genes that are 
#differentially expressed between the mice on low-fat diet and the mice on high-fat diet. 

#Filter out the comparison with p value equal to or larger than 0.05 
DEG1 = data.combined[data.combined[,2] < 0.05, ]
#Filter out the comparison with fold change larger than 0.8 or smaller than -0.8. 
DEG2= DEG1 [(DEG1 [, 1] > 0.8)|( DEG1 [, 1] < -0.8), ]
#291 genes were obtained
nrow(DEG2)
[1] 291

#Task 7 Functional Annotation of the Differentially Expressed Genes


#Using topGO package available in Bioconductor, compute enriched GO 
#Biological Process (BP) terms. Hint: You will have to create topGOdata object, use it as an argument in 
#runTest function, and read the result using GenTable function. 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("topGO")
library(topGO)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
genelist<-DEG2[,2]
GOterms <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="symbol")
class(GOterms)
selection <- function(allScore){ return(allScore < 0.05)}
GOdata <- new("topGOdata",ontology="BP",allGenes=genelist,annot=annFUN.GO2genes,GO2genes=GOterms ,geneSel=selection, nodeSize=10)
results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)