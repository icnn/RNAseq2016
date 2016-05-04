## PC_Corr.R
## Script to Compare expression principal componenets to datMeta covariates
## Thanks and credit to Neel Parikshak and Michael Gandal for original work done with this script

options(stringsAsFactors = FALSE)
load(file="rlogNorm_OutlierRemovedData_CM.RData")

thisdat <- t(scale(t(datExpr),scale=F)) 
## Centers the mean of all genes - this means the PCA gives us the eigenvectors of the geneXgene covariance matrix, allowing us to assess the proportion of variance each component contributes to the data
PC <- prcomp(thisdat);
topPC <- PC$rotation[,1:5];
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topvar <- varexp[1:5]
colnames(topPC) <- paste(colnames(topPC)," (",signif(100*topvar[1:5],2),"%)",sep="")

## Do a PCA of the sequencing statistics

datSeq_compare <- data.matrix(datSeq[,c(3:28)])
datSeq_compare = datSeq_compare[,-c(2,3,4,8)]

datSeqNorm <- t(scale(datSeq_compare,scale=F))
PC.datSeq <- prcomp(datSeqNorm);
varexp <- (PC.datSeq$sdev)^2 / sum(PC.datSeq$sdev^2)
print(varexp[1:5])
topPC.datSeq <- PC.datSeq$rotation[,1:2]; ## Explains nearly 100% of variance in datSeq
colnames(topPC.datSeq) <- c("SeqPC1","SeqPC2")
### PC1 - amount of reads (coverage)
### PC2 - how well these reads aligned to genome/GC content

pairsdat <- data.frame(Diagnosis=as.factor(datMeta[,"Dx"]),
                       age=as.factor(datMeta[,"Age_of_Death"]),sex=as.factor(datMeta[,"Gender"]),
                       race=as.factor(datMeta[,"Ethnicity"]),PMI=as.factor(datMeta[,"PMI_hrs"]),RIN=as.factor(datMeta[,"RIN"]))

library(gplots)

blackcol <- paste(col2hex("black"),"30",sep="")
redcol <- paste(col2hex("red"),"30",sep="")
greencol <- paste(col2hex("green"),"30",sep="")

colCondition <- datMeta[,"Dx"]; colCondition <- as.character(colCondition) 
colCondition[colCondition=="BP"] <- blackcol; 
colCondition[colCondition=="Control"] <- redcol; 
colCondition[colCondition=="SCZ"] <- greencol;

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs",method="spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

Total_Reads = datSeq$Total_Reads

pdf("./SeqStats_Comparison_OutlierRm_seqPCs.pdf",height=20,width=24)
pairs(cbind(datSeq_compare,topPC.datSeq),col=colCondition,pch=19,upper.panel = panel.cor,main="Sequencing Statistics and Top datSeq PC Comparison (|Spearman's rho| correlation values)")
pairs(cbind(topPC,topPC.datSeq,Total_Reads),col=colCondition,pch=19,upper.panel = panel.cor,main="Top datExpr PC and Top datSeq PC Comparison (|Spearman's rho| correlation values)")
pairs(cbind(pairsdat,topPC,topPC.datSeq),col=colCondition,pch=19,upper.panel = panel.cor,main="datMeta, Top PC datSeq and Top PC datExpr Comparison (|Spearman's rho| correlation values)")
pairs(cbind(pairsdat,topPC,topPC.datSeq,datSeq_compare),col=colCondition,pch=19,upper.panel = panel.cor,main="All Comparison (|Spearman's rho| correlation values)")
dev.off()

datMeta <- cbind(datMeta,topPC.datSeq)
## Use seqStats in your linear model

save(datExpr,datMeta,datSeq,file="FinalProcData_CM.RData")

