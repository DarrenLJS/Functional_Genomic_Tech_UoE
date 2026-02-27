#load the required libraries for the tutorials... 
library(affy) 
library(limma) 
library(mouse4302.db) 
library(annotate)

#check the tutorial data has loaded into the current 
#working directory. 
dir() 
#You now need to load the .Rdata file into 
#your workspace to recover the saved files 

load("enrichment.Rdata") 

#Note: In theory you could use 
#load("Mm.h.all.v7.1.entrez.rds") 
#But load of this type of file is not supported by all 
#versions of R- so we use readRDS instead 

Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds")  

#Check that you have the required objects 
ls() 

#Show the full contents of the annotation package 
ls("package:mouse4302.db") 

#Show the annotation keys in this database 
keytypes(mouse4302.db)  

#Mapping from 
#http://genomicsclass.github.io/book/pages/mapping_features.ht 
#ml 
sampleNames(eset) 

#Here we select from the annotation a number of keys with the primary key being PROBEID 
res <- select(mouse4302.db, keys = rownames(eset), columns = 
                c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID") 
#View the top of the table 
head(res) 
#find the index of each row of the expression set in the 
#annotation object res 
idx <- match(rownames(eset), res$PROBEID) 
#Use the index to set the phenotypic data in the ExpressionSet 
fData(eset) <- res[idx, ] 
head(fData(eset), 10) 
#Find all rows that don’t have an EntrezID and remove then 
eset_t<-eset[is.na(fData(eset)$ENTREZID)==0,] 

#convert to indexes 
H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID) 


#Pick the most suitable enrichment analysis tool to find 
#enrichment signatures in the data and run this tool So:-if 
#you want to run mroast 
#results <-mroast(eset_t, 
#                 index=H.indices, 
#                 design=design, 
#                 contrast=contrastmatrix[,1], 
#                 adjust.method = "BH") 
#if you want to run camera 
#results <-camera(eset_t, 
#                 index=H.indices, 
#                 design=design, 
#                 contrast=contrastmatrix[,1]) 
#if you want to run romer (takes 5mins or so) 
#results <-romer(eset_t, 
#                index=H.indices, 
#                design=design, 
#                contrast=contrastmatrix[,1] )

run_enrichment <- function(method, eset, indices, design, contrast) {
  
  fun <- switch(method,
                mroast = mroast,
                camera = camera,
                romer  = romer)
  
  fun(eset,
      index = indices,
      design = design,
      contrast = contrast)
}

results <- run_enrichment("camera",
                          eset_t,
                          H.indices,
                          design,
                          contrastmatrix[,1])

#View the results 
results 
#Use help for other parameters. Note we might decide to use 
#exactly the same model as our differential gene analysis for 
#the enrichment analysis- in this case we can extract it from 
#the fit 
#sv <- squeezeVar(fit$sigma^2,df=fit$df.residual)

write.table(results,"enrichment.txt",sep="\t") 
#You can then examine the results in “enrichment.txt”.  It is 
#a text file.  It can be downloaded to view in a spreadsheet 
#such as Excel.