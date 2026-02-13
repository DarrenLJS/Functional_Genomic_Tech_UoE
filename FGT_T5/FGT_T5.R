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