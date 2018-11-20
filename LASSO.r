library("glmnet")
library("survival")
myexpr <- read.csv("easy_input_exp.csv",header = T,row.names = 1,check.names = F)
mysurv <- read.csv("easy_input_suv.csv",header = T,row.names = 1)
mysurv$futime <- mysurv$futime/365 #convert days to years
cvfit = cv.glmnet(t(myexpr), Surv(mysurv$futime,mysurv$fustat), nfold=10, 
                  family = "cox"   
)  #LASSO regression with 10 folds cross validation
plot(cvfit)
fit <- glmnet(t(myexpr), Surv(mysurv$futime,mysurv$fustat), 
              family = "cox")   
plot(fit, label = TRUE)
cvfit$lambda.min #check best lambda
coef.min = coef(cvfit, s = "lambda.min") 
coef.min
active.min = which(coef.min != 0)
eneids <- rownames(myexpr)[active.min]
geneids #check significant genes
index.min = coef.min[active.min]
index.min #check index
combine <- cbind(geneids, index.min)
write.csv(combine,"gene_index.csv")
signature <- as.matrix(t(myexpr[geneids,])) %*% as.matrix(index.min) 
summary(signature)
colnames(signature)[1] <- "lasso"
row.names = row.names(myexpr)
write.table(signature,"lasso_output.txt",row.names = T, quote = F)
