# Commands to download and install the repackage 
# These packages are already installed on the class computer 
# Copy this code and uncomment the three commands below to 
# install the required files on another machine 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager") 
# BiocManager::install("affy") 
# BiocManager::install("limma") 

# load the required libraries for the tutorials... 

library(affy) 
library(limma) 

mydata <- ReadAffy()
mydata

png("densityplot1.png")
# Quality control plots 
hist(mydata) 
dev.off()

png("boxplot1.png")
# And a boxplot with different colour per sample group 
colours <- c(rep("yellow",3),rep("red",2),rep("blue",2), 
             "red", rep("green",3)) 
boxplot(mydata, col=colours, las=2) 
dev.off()

# Normalise the data using RMA 
eset <- rma(mydata) 
eset 
# To obtain a matrix of the expression values, use exprs()  
values <- exprs(eset)

png("boxplot2.png")
# Boxplot to observe the results of normalisation 
# Notice differences with the boxplot from the raw data 
boxplot(values, col=colours, las=2) 
dev.off()

# MA plot of the samples 1, 2 and 4 
png("MVAplot.png")
mva.pairs(values[, c(1,2,4)]) 
# The same plot for the non-normalised raw data 
# Note that the mva.pairs call below only plots a few of the  
# samples – you may wish to plot them all but this is slow 
mva.pairs(pm(mydata)[, c(1,2,4)])
dev.off()

# Load the target file into an AnnotatedDataFrame object 
adf <- read.AnnotatedDataFrame(
  "FGT_T3_targets.txt", 
  header=TRUE, 
  row.names=1, 
  as.is= TRUE
) 
# To facilitate interpretation, let’s replace the columns  
# header,currently 
# displaying the filename, to show the name of each sample  
# (if you have a targets file) 
colnames(values) <- rownames(pData(adf)) 
# Performs hierarchical clustering with average linkage based on 
# Pearson’s Correlation Coefficient 
png("hierarchicalcluster.png")
hc <- hclust(as.dist(1-cor(values, method="pearson")), method="average") 
plot(hc)
dev.off()

# here we request a specific package from a specific archive... 
# install.packages("scatterplot3d",repo="http://cran.ma.imperial.ac.uk") 

# Then load the library 
library(scatterplot3d) 
# Perform PCA 
pca <- prcomp(t(values), scale=T) 
# Plot the PCA results 
png("3dpca.png")
s3d <- scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(1)) 
s3d.coords <- s3d$xyz.convert(pca$x[,1:3]) 
text(s3d.coords$x, s3d.coords$y, labels = colnames(values), 
     pos = 3,offset = 0.5)
dev.off()

# obtaining a matrix of expression values 
exprsvals <- exprs(eset) 
# RMA outputs log2 data while MAS5 outputs linear data 
# To convert from log… 
exprsvals10 <-2^exprsvals 
# check conversion 
exprsvals[1:100,] 
# converted 
exprsvals10[1:100,] 

# check order of sample names 
mysamples <- sampleNames(eset) 
# display the list 
mysamples 
# it is useful to obtain a vector of ProbeIDs here 
probesets <- probeNames(mydata) 
# display the first 100 ProbeSets 
probesets[1:100]

# Calculate the means 
# Note mean of the log is not the same as the log of the mean!! 
ES.mean <- apply(exprsvals10[,c("GSM272753.CEL", 
                                "GSM272836.CEL","GSM272837.CEL")],1,mean) 
iPS_OK.mean <- apply(exprsvals10[,c("GSM272839.CEL", 
                                    "GSM272846.CEL","GSM272890.CEL")],1,mean) 
iPS_4F.mean <- apply(exprsvals10[,c("GSM279200.CEL", 
                                    "GSM279201.CEL","GSM279202.CEL")],1,mean) 
NSC.mean <- 
  apply(exprsvals10[,c("GSM272847.CEL","GSM272848.CEL")],1,mean) 
# calculate some fold changes 
ES_iPS_OK <-ES.mean /iPS_OK.mean 
ES_iPS_4F <-ES.mean /iPS_4F.mean 
ES_NSC <-ES.mean /NSC.mean 
# build a summary table to hold all the data 
all.data= cbind(ES.mean,iPS_OK.mean,iPS_4F.mean, NSC.mean, 
                ES_iPS_OK, 
                ES_iPS_4F, ES_NSC) 
# check the column names 
colnames(all.data) 
# write the table of means as an output 
write.table(all.data,file="group_means.txt", quote=F, 
            sep="\t",col.names=NA) 

