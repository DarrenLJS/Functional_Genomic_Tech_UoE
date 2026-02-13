library(affy) 
library(limma) 
library(mouse4302.db) 
library(annotate)

load("fold_filtering.Rdata")

#Check original sample order 
sampleNames(eset) 
#Rename the samples 
sampleNames(eset) <- c(
  "ESC.1","ESC.2","ESC.3","iPS2.2","iPS2.3",
  "NSC.1","NSC.2","iPS2.1","iPS4.1","iPS4.2","iPS.3"
) 
#Check the samples have renamed 
sampleNames(eset) 

eset@annotation 

library(mouse4302.db)# load chip-specific annotation 
#packages in the annotation package 
ls("package:mouse4302.db") 

#build an annotation table 
ID <- featureNames(eset) 
Symbol <- getSYMBOL(ID, "mouse4302.db") 
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME")) 
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, 
                  stringsAsFactors=F) 
tmp[tmp=="NA"] <- NA #fix padding with NA characters  
#assign as feature data of the current Eset 
fData(eset) <- tmp 

#Build the design matrix 
design <- model.matrix(~-1+factor(c(1,1,1,2,2,3,3,2,4,4,4))) 
colnames(design) <- c("ESC","iPS2","NSC","iPS4") 
#Check it makes sense 
sampleNames(eset) 
#output the design matrix 
design 

#This instructs Limma which comparisons to make 
contrastmatrix <- makeContrasts(ESC-iPS2,ESC-NSC,ESC-iPS4, 
                                levels=design) 
contrastmatrix 

#issue these commands to fit the model 
#and make the contrasts 
fit <- lmFit(eset, design) 

fit2 <- contrasts.fit(fit, contrastmatrix) 

#this last part essentially moderates the t-statistic using  
#the borrowed variance approach described in class 
fit2 <- eBayes(fit2) 

topTable(fit2,coef=1,adjust="fdr") 
myresults <- topTable(fit2,coef=1, adjust="fdr", 
                     number=nrow(eset)) 
write.table(myresults,"myresults.txt") 

clas <- classifyTestsF(fit2) 
vennDiagram(clas) 

# Volcano plot
library(ggplot2)

# Add a significance column
myresults$significant <- ifelse(myresults$adj.P.Val < 0.05 & abs(myresults$logFC) > 2, 
                           "Significant", "Not Significant")

ggplot(myresults, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("Not Significant" = "gray", 
                                "Significant" = "red")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10(P-value)") +
  theme_minimal()

ggplot(myresults, aes(x = AveExpr, y = logFC, color = adj.P.Val < 0.05)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "MA Plot",
       x = "Average Expression",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal()

