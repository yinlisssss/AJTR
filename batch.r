library(sva)
library(limma)
rt <- read.table("bat.txt",sep="\t",header=T,check.names=F) #the preprocessed input file of the four genes from GSE39279 and GSE75008
rt<-as.matrix(rt)
rownames(rt)<-rt[,1]
exp<-rt[,2:ncol(rt)]
dimnames<-list(rownames(exp),colnames(exp))
data<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType<-c(rep(1,56),rep(2,122)) # the first batch contains 56 samples, the second batch containes 122 samples
modType<-c(rep("normal",40),rep("tumor",16),rep("normal",0),rep("tumor",122)) #sample type, normal or tumor
mod <- model.matrix(~as.factor(modType)) 
outTab<-ComBat(data, batchType, mod, par.prior=TRUE) #batch normalization
outTab<-rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)

