#!/bin/bash

#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 4

source /mnt/software/Modules/current/init/bash
module load gatk/4.0.6.0
GATK="gatk --java-options -Xmx8G"

# required $GVCF $CHR $PLOIDY
# $GVCF -> *.g.vcf.gz produced by GATKHC
# $CHR -> chr
# $PLOIDY -> expected ploidy
# $REF -> reference fasta, requires associated sequence dict (ref.dict)

echo "Calling variants for ${CHR} reads from ${GVCF} against ${REF}, and ploidy expectation ${PLOIDY}."

$GATK GenotypeGVCFs \
 --reference "${REF}" \
 --variant "${GVCF}" \
 --output "${GVCF%.g.vcf.gz}.vcf.gz" \
 --sample-ploidy "${PLOIDY}" \
 --intervals "${CHR}" \
 --annotation-group StandardAnnotation \
 --annotation-group AS_StandardAnnotation \
 --annotation-group StandardHCAnnotation \
 --standard-min-confidence-threshold-for-calling 2.0