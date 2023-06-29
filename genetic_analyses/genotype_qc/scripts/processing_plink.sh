#!/bin/bash
################################################################################
# PROCESSING PIPELINE FOR mQTL AND RELATED ANALYSIS
# Used plink 1.90 beta
################################################################################
while getopts "d:o:" opt; do
  case ${opt} in
    d ) file_dir="$OPTARG"
    ;;
    o ) output_name="$OPTARG"
    ;;
    \? ) echo "Usage: processing_plink.sh -d <binary-directory> -o <output_basename>"
    ;;
  esac
done


# merge files (split across chromosome)

ls -d $file_dir/* | 
  grep -E "bed|bim|fam" | 
  sed 's/\..*//g' | 
  sort |
  uniq > "$output_name".manifest.txt

plink --merge-list "$output_name".manifest.txt --out "${output_name}_merged"

cur_step="${output_name}_merged"


# remove samples and snps with missing data

plink --bfile "${cur_step}" \
  --geno 0.01 \
  --mind 0.05 \
  --make_bed \
  --out "${cur_step}_missingness"
cur_step="${cur_step}_missingness"

# maf filter
plink --bfile "${cur_step}" \
  --maf 0.05 \
  --make-bed \
  --out "${cur_step}_maf"
cur_step="${cur_step}_maf"

# Hardy-Weinberg Disequilibrium
plink2 --bfile "${cur_step}" \
  --hwe 0.001 \
  --make-bed \
  --out "${cur_step}_hwe"
cur_step="${cur_step}_hwe"

# Sex validation
plink --bfile "${cur_step}" \
    --check-sex \
    --out "${cur_output}_sex_check"

