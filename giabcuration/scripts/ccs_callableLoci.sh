#!/usr/bin.sh
## Script used to generate callableLoci for use in calculating expected exanded region size. 
GATK=/scratch/nolson/sw/GenomeAnalysisTK.jar
REF=/scratch/nolson/giab/ref_genomes/hs37d5.fa
BAM=/scratch/nolson/giab/HG002-data/HG002.PacBio.15kbCCS.Q20.hs37d5.pbmm2.MAPQ60.HP10xtrioRTG.bam
OUTROOT=data/benchmark_extend/callableLoci/HG002_PacBio_15kbCCS

## GATK version info
java -jar ${GATK} --version

## Running CallableLoci for full genome
java -Xmx64G -jar ${GATK} \
	-T CallableLoci \
	-R ${REF} -I ${BAM} \
	--summary ${OUTROOT}.summary -o ${OUTROOT}.bed \
	-log ${OUTROOT}.log --generate_md5 \
	-minDepth 20 -mmq 20 -maxDepth 52 --allow_potentially_misencoded_quality_scores


## Calc MD5 sums for reproducibility
md5sum ${REF}
md5sum ${BAM}
md5sum ${OUTROOT}.summary
md5sum ${OUTROOT}.bed
