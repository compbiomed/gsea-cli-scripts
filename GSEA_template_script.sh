#!/bin/bash
# Run GSEA using command-line interface

##### Automatically defined variables

# Project prefix, e.g., "yyyy-mm-dd_Investigator"; defaults to current working directory
project_prefix="$(basename $(pwd))"
# MSigDB version
msigdb_version="6.0"
# SCC project account; defaults to current active group
project_account="$(groups | cut -f1 -d' ')"
# Paths and filenames
gsea_path="/restricted/projectnb/cbmhive/GSEA"
qsub_file="$gsea_path/run_GSEA.qsub"
project_path="/restricted/projectnb/$project_account/${project_prefix/-*/}/$project_prefix"
param_file="$gsea_path/GSEAParameters_v$msigdb_version.txt"

##### Edit the following for each pairwise comparison

# Full path to .rnk file
rnk_file="$project_path/${project_prefix}_*.rnk"
# Textual description of comparison, e.g., "knockout vs wildtype"
comparison=""
# Name for the queue job/output log
job_name="GSEA_$(basename ${rnk_file/.rnk/})"
# Run GSEA using specified variables
qsub -P $project_account -N $job_name -o $job_name.log $qsub_file $rnk_file $param_file "$comparison"
