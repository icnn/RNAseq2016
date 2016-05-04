# RNAseq2016
Repository for code shared at the RNAseq Symposium, 5/4/2016.

**FindConfounds.R** -> Look for confounds in datMeta of the CommonMind data.  
**PC_Corr.R** -> See the correlation between datExpr principal componenets and covariates.  
**LM_ConfoundRemove.R** -> Conservatively remove confounds from datExpr and do a t-test with Diagnosis coefficient.  
**rlogNorm_OutlierRemovedData_CMshort.RData** -> datMeta, datSeq, and datExpr to use with PC_Corr.R. Normalized and Outlier Removed Common Mind Consortium Data.  
**FinalProcData_CMshort.RData** -> datMeta, datSeq, and datExpr to use with FindConfounds.R and LM_ConfoundRemove.R. Normalized and Outlier Removed Common Mind Consortium Data with sequencing statistic principal components.  
**datMeta_CoVariates_CM.pdf** -> Example output from FindConfounds.R.  
**SeqStats_Comparison_OutlierRm_seqPCs.pdf** -> Example output from PC_Corr.R  
**DiffExpr_Bipolar_CM.csv** -> Example output from LM_ConfoundRemove.R.  
