# Bulk_RNAseq_pipeline

Bulk RNAseq analysis pipeline using bash scripts and R.

HPC Server
  1. fastqc - QC
  2. fastp - trimming of adapters
  3. star - Alignment
  4. featureCounts - quantification
  5. multiQC - summary of statistics from steps above
  6. bash script 6.combine_counts.sh - combine results from featureCounts into one matrix

R - local machine or Server

  1. preprocessing - stats QC, gene ID annotation 
  2. DEseq2
  3. functional pathways
