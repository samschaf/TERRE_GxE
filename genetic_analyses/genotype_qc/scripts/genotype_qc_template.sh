#!/bin/bash

#########################################################
#Pre-Imputation QC template
#########################################################
# Inputs:
  # b -> base_data String: prefix for plink 1.9 binary base dataset
  # o -> output_data String: directory
  # c -> chkpt: checkpoints file to continue a previous run at a specified endtime
# Outputs:
  #output_name.chkpt: list of filenames to recontinue checkpoints
while getopts "b:o:" opt; do
  case ${opt} in
    b ) base_data="$OPTARG"
    ;;
    o ) output_data="$OPTARG"
    ;;
    \? ) echo "Usage: processing_plink.sh -b <binary-dataset> -o <output-dataset>"
    ;;
  esac
done

# 1. Pre-cleaning
pre_cleaning(){
  cur_geno=$1
  plink --bfile "${cur_geno}" \
    --make-bed \
    --geno 0.05 \
    --out "${cur_geno}_snpcall"
  cur_geno="${cur_geno}_snpcall"
}
