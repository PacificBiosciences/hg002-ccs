
## Generate bed with expanded regions
STRATBEDDIR=/scratch/nolson/giab/benchmarking-tools/resources/stratification-bed-files/LowComplexity
## CCS callable subtracting AllRepeats_gt95percidentity_slop5.bed.gz and SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz from callableLoci bed
ALLREPEATS=${STRATBEDDIR}/AllRepeats_gt95percidentity_slop5.bed.gz
IMPHOMO=${STRATBEDDIR}/SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz 
CALLABLE=data/benchmark_extend/callableLoci/HG002_PacBio_15kbCCS.bed
OUTROOT=data/benchmark_extend
CCSCALLABLE=${OUTROOT}/HG002_PacBio_15kbCCS_callable.bed
BENCHMARK=/scratch/nolson/giab/HG002-data/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed

bedtools --version


## Subtract stratification beds from CCS callable
bedtools subtract -a ${CALLABLE} -b ${ALLREPEATS} | \
    bedtools subtract -a - -b ${IMPHOMO} >  ${CCSCALLABLE} #Check for STDIN

## Substracting benchmark from callable
bedtools subtract -a ${CCSCALLABLE} -b ${BENCHMARK} > ${OUTROOT}/ccs_extended.bed 

## Substract structural variants and superdups merged, maybe selfchain


## MD5s for reproducibility
md5sum ${ALLREPEATS}
md5sum ${IMPHOMO}
md5sum ${CALLABLE}
md5sum ${CCSCALLABLE}
md5sum ${BENCHMARK}
md5sum ${OUTROOT}/ccs_extended.bed