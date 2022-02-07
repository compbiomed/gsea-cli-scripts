#!/bin/bash
# Run GSEA using command-line interface

# Global variables #############################################################

# Project prefix, e.g., "yyyy-mm-dd_Investigator"
# Defaults to current working directory
project_prefix="$(basename $(pwd))"
# SCC project account; defaults to current active group
project_account="$(groups | cut -f1 -d' ')"

# GSEA version
# gsea_version="4.1.0" # New script-based version, not yet working right on SCC
gsea_version="2.2.1"   # Old .jar-based version
# MSigDB version
msigdb_version="7.5.1"

# By default, look in current working directory for .rnk files; change as needed
rnk_path="."

# Set environment variable to path to GSEA files; this will be passed to the
# queue job with the -V qsub flag so that it can find the GSEA executables
export GSEA_PATH="/rprojectnb/cbmhive/GSEA"
# Define filenames
qsub_file="${GSEA_PATH}/run_GSEA.qsub"
param_file="${GSEA_PATH}/GSEAParameters_v${msigdb_version}.txt"

# Regular expressions for phrases that should be removed from .rnk filenames
# in order to extract a textual description of the pairwise comparison
comparison_regexes=(
  "^${project_prefix}_"
  "_[Mm]oderated_t$" "_Student_t$" "_t$" "_Wald_statistic$" "_log2_fold_change$"
)

# Define the main effects in the model if interaction effect was used for GSEA
# For example:
# effects=("genotype" "treatment")
effects=()

# Submit all *.rnk files in specified directory ################################

# Resolve path name
rnk_path="$(realpath ${rnk_path})"

for rnk_file in "${rnk_path}"/*.rnk
do
  # Extract prefix from .rnk filename
  rnk_file_prefix="$(echo "$(basename ${rnk_file})" | sed -r 's/\.rnk$//')"

  # Extract textual description of comparison (e.g., "knockout vs wildtype")
  # from .rnk filename prefix
  comparison="${rnk_file_prefix}"
  # Remove any specified phrases
  for regex in ${comparison_regexes[@]}
  do
    comparison="$(echo "${comparison}" | sed "s/${regex}//")"
  done
  # Replace any underscores with spaces
  comparison="$(echo "${comparison}" | sed "s/_/ /g")"
  # An interaction effect is labeled differently
  # Note that this code still works even if the 'effects' array is not defined
  if [[ "${comparison}" = "${effects[@]}" ]]
  then
    comparison="interaction of ${effects[0]} and ${effects[1]}"
  fi

  # Run GSEA using specified variables
  # Note: this job is run in parallel with 8 processors
  #       since it will use more than 1 processor of power otherwise,
  #       and the process reaper will kill it.
  # (This was originally set to 4 processors, but the process reaper killed two
  #  runs each using > 5 CPUs.) 2016-06-06 ACG
  job_name="GSEA_${rnk_file_prefix}"
  qsub -P ${project_account} -N ${job_name} -o ${job_name}.log -pe omp 8 -V \
       "${qsub_file}" \
       ${gsea_version} "${rnk_file}" "${param_file}" "${comparison}"
done
