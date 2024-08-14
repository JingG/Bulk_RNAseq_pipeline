# Bulk_RNAseq_pipeline

Bulk RNAseq analysis pipeline using bash scripts and R.

HPC Server
  1. fastqc - QC
  2. fastp - trimming of adapters
  3. star - Alignment
  4. featureCounts - quantification
  5. multiQC - summary of statistics from steps above

R - local machine or Server

  6. stats QC
  7. DEseq2
