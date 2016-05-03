## LM_ConfoundRemove_Short.R

options(stringsAsFactors = FALSE)
setwd("C:/Users/jillh/Dropbox/GitHub/RNAseq/")
load("C:/Users/jillh/Dropbox/DHGLab/commonmind/FinalProcData_CM.RData")

datMeta$Dx <- factor(datMeta$Dx, levels=c("Control","BP","SCZ"))
datMeta$Ethnicity <- factor(datMeta$Ethnicity, levels=c("Caucasian","African-American","Hispanic","Asian"))

## Get the covariate data
mDx <- model.matrix(~datMeta$Dx)
mDx <- mDx[,-1]
age <- as.numeric(as.factor(datMeta[,"Age_of_Death"]))
sex <- as.numeric(as.factor(datMeta[,"Gender"]))-1
RIN <- as.numeric(as.character(datMeta[,"RIN"]))
race <- model.matrix(~datMeta$Ethnicity)
race <- race[,-1]
PMI <- as.numeric(as.factor(datMeta$PMI_hrs))
seqStatPC1 <- as.numeric(datMeta$`SeqPC1 (98.7%) - Depth`)
seqStatPC2 <- as.numeric(datMeta$`SeqPC2 (<2%) - GC/Length`)
regvars <- as.data.frame(cbind(mDx,age,sex,RIN,race,seqStatPC1,seqStatPC2))
colnames(regvars)[1:2]=c("Bipolar","Schizo")

## Make object for ttables
ttable = list("Bipolar"=matrix(NA,nrow=nrow(datExpr),ncol=5),
              "Schizo"=matrix(NA,nrow=nrow(datExpr),ncol=5))
rownames(ttable$Bipolar) <- rownames(ttable$Schizo) <- rownames(datExpr)
colnames(ttable$Bipolar) <- colnames(ttable$Schizo) <- c("Estimate","Std.Error","t.value","Pr(>|t|)","p.adj")

## Run the regression
datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

lmmod <- apply(as.matrix(datExpr),1,function(y) lm(y~age+sex+RIN+race+PMI+seqStatPC1+seqStatPC2,data=regvars))

for (i in 1:nrow(datExpr)) {
  if (i%%1000 == 0) {print(i)}
  
  coef <- coef(lmmod[[i]])
  datExpr.reg[i,] <- coef[1] + residuals(lmmod[[i]]) 
  
  lmmod2 <- lm(datExpr.reg[i,]~Bipolar+Schizo,data=regvars)
  lmmod2$df.residual <- lmmod2$df.residual - 9   ## 9 the number of coefficients previously removed
  summary = summary(lmmod2)
  
  ttable$Bipolar[i,1:4] <- coef(summary)[2,]
  ttable$Schizo[i,1:4] <- coef(summary)[3,]
  
  ## The full data minus all of the covariates but condition 
  ## Also equivalent to <- datExpr - coef*var (unwanted variates) model
  ## So we have a datExpr object with only the contributions of Dx
  
  ## Then adjust the degrees of freedom for the residual before we do our t tests
  ## Ignore the warning
}


## Get adjusted p values
for(i in c(1:2)){
  ttable[[i]][,5] <- p.adjust(ttable[[i]][,4],method="fdr")
}

write.csv(ttable$Bipolar,file="DiffExpr_Bipolar_CM.csv")
write.csv(ttable$Schizo,file="DiffExpr_Schizo_CM.csv")

quantile(datExpr[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(datExpr.reg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))

save(datMeta,datSeq,datExpr.reg,file="datExpr.reg_CM.RData")