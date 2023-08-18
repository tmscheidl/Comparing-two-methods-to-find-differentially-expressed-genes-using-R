library(affy)
library(readr)
library(affyio)
library(genefilter)
library(limma)
library(farms)
library(dplyr)
library(ggrepel)
library(arsenal)
library(ggplot2)
library("ggVennDiagram")



data <- ReadAffy() 
rma.data <- rma(data)
file<- exprs(rma.data)

#C = list.files(pattern = "C-")
#D = list.files(pattern = "D-")

#data_group <- file   
#group_by(d1.CEL, d2.CEL,d3.CEL,d4.CEL,d5.CEL,d6.CEL,d7.CEL,d8.CEL,d9.CEL) %>%
#group_by(c1.CEL,c2.CEL,c3.CEL,c4.CEL,c5.CEL,c6.CEL) %>%
  
#as.data.frame()
#data_group
#group_split(data_group)

#m<- matrix(file)
#m <- mapply(m, FUN=as.numeric)

C_group <-file
C_group<- subset (C_group, select = -c(d1.CEL:d9.CEL))

D_group <-file
D_group<- subset (D_group, select = -c(c1.CEL:c6.CEL))

filter<- nsFilter(rma.data)
esetnew<-filter[[1]]
stats<-filter[[2]]

#a1 <- varFilter(file)
#a2 <- featureFilter(m)

Grp<- c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)
eset<- rma(data)
dm <- model.matrix (~0+factor(Grp))
colnames(dm) <- c("C_group", "D_group")
fit <-lmFit(eset,dm)
fit2 <- eBayes(fit)
DEG <- topTable(fit2, coef = 1)

best50 <- topTable(fit2, number = 50)

rowttests_method<- rowttests(eset, as.factor(Grp))
best50rm <- topTable(rowttests_method, number = 50)


a <- rowttests_method$p.value
b <- rowttests_method$dm
plot(b, -log10(a), pch=".", xlim=c(-1.5,1.5),
     xlab="Fold change", ylab="p.value")

tt_rows <- rownames(rowttests_method) 
limma_rows <- rownames(fit2)

tt_rows %in% limma_rows
match (tt_rows,limma_rows)

z <- intersect(limma_rows,tt_rows)

ggVennDiagram(
  x = list(limma_rows, tt_rows)
)
