# Small Variant Benchmarking Results
Analysis of small variant benchmarkign results and GIAB benchmark set expansion estimates for the HG002 PacBio CCS 15Kb library.
Variant calling was performed using GATK and DeepVariant.
Callsets were benchmarked against the GIAB HG002 V3.3.2 benchmark set using the vcfeval-happy pipeline though the precisionFDA app. 

These results are part of a larger study:   
[__Highly-accurate long-read sequencing improves variant detection and assembly of a human genome__](https://www.biorxiv.org/content/10.1101/519025v2)

_Abstract_  
The major DNA sequencing technologies in use today produce either highly-accurate short reads or noisy long reads. We developed a protocol based on single-molecule, circular consensus sequencing (CCS) to generate highly-accurate (99.8%) long reads averaging 13.5 kb and applied it to sequence the well-characterized human HG002/NA24385. We optimized existing tools to comprehensively detect variants, achieving precision and recall above 99.91% for SNVs, 95.98% for indels, and 95.99% for structural variants. We estimate that 2,434 discordances are correctable mistakes in the high-quality Genome in a Bottle benchmark. Nearly all (99.64%) variants are phased into haplotypes, which further improves variant detection. De novo assembly produces a highly contiguous and accurate genome with contig N50 above 15 Mb and concordance of 99.998%. CCS reads match short reads for small variant detection, while enabling structural variant detection and de novo assembly at similar contiguity and markedly higher concordance than noisy long reads.
