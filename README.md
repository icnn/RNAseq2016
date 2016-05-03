# RNAseq2016
Repository for code shared at the RNAseq Symposium, 5/4/2016.

**FindConfounds.R** -> Look for confounds in datMeta of the CommonMind data.  
**LM_ConfoundRemove.R** -> Conservatively remove confounds from datExpr.  
**FinalProcData_CM.RData** -> datMeta, datSeq, and datExpr to use with FindConfounds.R and LM_ConfoundRemove.R. Normalized and Outlier Removed Common Mind Consortium Data.  
**datMeta_CoVariates_CM.pdf** -> Example output from FindConfounds.R.  
**DiffExpr_Bipolar_CM.csv** -> Example output from LM_ConfoundRemove.R.
