## FindConfounds.R
## Script to make covariate Plots - want to look at relationship of specific factors in datMeta with Diagnosis
## Thanks and credit to Neel Parikshak and Michael Gandal for original work done with this script

## Load your datMeta

options(stringsAsFactors = FALSE)
load("FinalProcData_CM.RData")

pdf(file="datMeta_Covariates.pdf")

par(mfrow=c(2,3))

plot(datMeta$Dx, ylab="Number", main="Subjects")

for(i in c(2,6:15,19)){  ## which factors in datMeta are you looking at?
  
  ## Specify which factors of interest are characters (categorical variables)
  if( i == 2 || i == 6 || i == 7 || i== 9 || i == 13){
    print(paste(i,"Character Graph",sep=" "))
    A = anova(lm(as.numeric(datMeta[,i]) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   
    plot(datMeta[,i] ~ datMeta$Dx, main=paste(colnames(datMeta)[i]," p=", signif(p,2)), ylab="", xlab="")
  }
  else{
    print(paste(i,"Number Graph",sep=" "))
    A = anova(lm(as.numeric(datMeta[,i]) ~ datMeta$Dx)); p = A$"Pr(>F)"[1];   
    plot(as.numeric(as.character(datMeta[,i])) ~ datMeta$Dx, main=paste(colnames(datMeta)[i]," p=", signif(p,2)), ylab="", xlab="")
  }
}

dev.off()
