library("glmnet")
library("survival")
myexpr <- read.csv("easy_input_exp.csv",header = T,row.names = 1,check.names = F)
mysurv <- read.csv("easy_input_suv.csv",header = T,row.names = 1)
cvfit = cv.glmnet(t(myexpr), Surv(mysurv$futime,mysurv$fustat), nfold=10,
                  family = "cox"
) 
plot(cvfit)
fit <- glmnet(t(myexpr), Surv(mysurv$futime,mysurv$fustat), 
              family = "cox")   #LASSO regression
plot(fit, label = TRUE)

cvfit$lambda.min #check lambda

coef.min = coef(cvfit, s = "lambda.min") 
coef.min
active.min = which(coef.min != 0)
geneids <- rownames(myexpr)[active.min]
geneids

index.min = coef.min[active.min]
index.min
combine <- cbind(geneids, index.min)
write.csv(combine,"gene_index.csv")
signature <- as.matrix(t(myexpr[geneids,])) %*% as.matrix(index.min) 
summary(signature)
colnames(signature)[1] <- "lasso"
row.names = row.names(myexpr)
write.table(signature,"lasso_output.txt",row.names = T, quote = F)
