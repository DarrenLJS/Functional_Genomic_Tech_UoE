library(affy)
library(pheatmap)

test <- read.table("afile.txt", sep = "\t", header = TRUE, row.names = 1, as.is = TRUE)
mat1 <- as.matrix(test)
exprset1 <- new("ExpressionSet", exprs = mat1)

mymatrix <- exprs(exprset1)
mymatrix[1:10]
mysnames <- sampleNames(exprset1)
mysnames


