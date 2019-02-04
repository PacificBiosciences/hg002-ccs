#!/bin/bash

#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 2

source /mnt/software/Modules/current/init/bash
module load bedtools/2.27.1 gatk/4.0.6.0

echo "Extracting HG19 PAR regions from ${DIP_VCF} and haploid regions from ${HAP_VCF} and merging them to form ${OUTVCF}."

# The Y chromosome in this assembly contains two pseudoautosomal regions (PARs) that were taken from the corresponding regions in the X chromosome and are exact duplicates:
# chrY:10001-2649520 and chrY:59034050-59363566
# chrX:60001-2699520 and chrX:154931044-155260560

read -r -d '' HG19_DIP << EOM
X       60001   2699520
X       154931044       155260560
EOM

read -r -d '' HG19_HAP << EOM
X       1       60000
X       2699521 154931044
X       155260561       155270560
EOM

bedtools intersect -header -a "${HAP_VCF}" -b <(echo "${HG19_HAP}") > "${HAP_VCF%.*}.tocombine.vcf"
bedtools intersect -header -a "${DIP_VCF}" -b <(echo "${HG19_DIP}") > "${DIP_VCF%.*}.tocombine.vcf"

$GATK MergeVcfs \
    --INPUT "${HAP_VCF%.*}.tocombine.vcf" \
    --INPUT "${DIP_VCF%.*}.tocombine.vcf" \
    --OUTPUT "${OUTVCF}"