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

# required $BAMLIST $REF $SAMPLE $CHR $PLOIDY $PCRINDELMODEL $MAPQ
# $BAMLIST -> BAM or list of BAM file paths (named bam.list)
# $REF -> reference fasta, requires associated sequence dict (ref.dict)
# $SAMPLE -> sample name from BAM
# $CHR -> chr
# $PCRINDELMODEL -> ['NONE', 'HOSTILE', 'AGGRESSIVE', 'CONSERVATIVE']
# $MAPQ -> minimum MAPQ
# $PLOIDY -> expected ploidy

echo "Calling variants from interval ${CHR} of ${BAMLIST} against ${REF}.  Using pcr-indel-model ${PCRINDELMODEL}, with an expected ploidy of ${PLOIDY}, on alignments with MAPQ >= ${MAPQ}."

$GATK HaplotypeCaller \
 --reference "${REF}" \
 --input "${BAMLIST}" \
 --output "${SAMPLE}.${CHR}.${PLOIDY}ploid.${PCRINDELMODEL}.HaplotypeCaller.g.vcf.gz" \
 --sample-ploidy "${PLOIDY}" \
 --pcr-indel-model "${PCRINDELMODEL}" \
 --intervals "${CHR}" \
 --read-filter MappingQualityReadFilter \
 --read-filter NotSecondaryAlignmentReadFilter \
 --read-filter NotSupplementaryAlignmentReadFilter \
 --minimum-mapping-quality "${MAPQ}" \
 --emit-ref-confidence GVCF \
 --annotation-group StandardAnnotation \
 --annotation-group AS_StandardAnnotation \
 --annotation-group StandardHCAnnotation