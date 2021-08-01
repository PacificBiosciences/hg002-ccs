#!/usr/bin/env python3
"""Measure concordance of read alignments to a reference genome"""

from __future__ import division
import argparse
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
signal.signal(signal.SIGINT, signal.SIG_DFL)
collections = None
math = None
pysam = None

"""Usage"""
argparser = argparse.ArgumentParser("""Measure concordance of read alignments to a reference genome""",
    epilog="""
CONCORDANCE
    Measure the concordance of alignments to the reference genome.
    Consider only non-variant positions in high-confidence regions.
""",
    formatter_class=argparse.RawDescriptionHelpFormatter
    )
argparser.add_argument("reffasta", metavar="ref.fasta")
argparser.add_argument("inbam", metavar="in.bam")
argparser.add_argument("outcsv", metavar="out.csv", help="Per aligment concordance statistics, CSV")
argparser.add_argument("--hcregions", metavar="hcregions.bed.gz", help="High-confidence regions, BED.gz (.tbi required)")
argparser.add_argument("--hcvariants", metavar="hcvariants.vcf.gz", help="High-confidence variants, VCF.gz (.tbi required)")
argparser.add_argument("--chrom", metavar="S", help="Limit analysis to the specified chromosome")
argparser.add_argument("--outbam", metavar="out.bam", help="Input BAM with concordance annotations for relevant records")


# Enum for CIGAR operations (pysam convention)
MATCH = 0
INS = 1
DEL = 2
REF_SKIP = 3
SOFT_CLIP = 4
HARD_CLIP = 5
PAD = 6
EQUAL = 7
DIFF = 8

# Groups of related CIGAR operations
QUERYOPS = (MATCH, EQUAL, DIFF, INS) # operations that include query bases
REFOPS = (MATCH, EQUAL, DIFF, DEL)   # operations that include reference bases
CLIPPINGOPS = (HARD_CLIP, SOFT_CLIP) # operations that clip query bases


class AlignmentStats:
    def __init__(self, alignment):
        self.alignment = alignment
        self.matchBp = 0
        self.mismatchBp = 0
        self.nonHpInsertionBp = 0
        self.nonHpDeletionBp = 0
        self.hpInsertionBp = 0
        self.hpDeletionBp = 0
        self.ignoredBp = 0

    @property
    def insertionBp(self):
        return self.nonHpInsertionBp + self.hpInsertionBp

    @property
    def deletionBp(self):
        return self.nonHpDeletionBp + self.hpDeletionBp

    @property
    def errorsBp(self):
        return self.mismatchBp + self.deletionBp + self.insertionBp

    @property
    def hcReadLength(self): # length of read that fell within high-confidence, non-variant positions
        return self.mismatchBp + self.matchBp + self.deletionBp

    @property
    def concordance(self): # M/(M+X+D+I) = read/reference symmetric concordance
        return 1. * self.matchBp / (self.matchBp + self.mismatchBp + self.deletionBp + self.insertionBp)

    @property
    def concordanceQv(self): # concordance in Phred-scale, capped by hcReadLength
        maxQv = 10*math.log10(self.hcReadLength+1) # QV capped by read length with pseudocount of 1
        qv = -10*math.log10(1-self.concordance) if self.errorsBp else maxQv
        return min(maxQv, qv)


def measureConcordance(al, reffasta, hcRegions, hcVariants):
    """Measure the concordance of the alignment to the reference."""

    # Walk through each CIGAR block in the alignment
    stats = AlignmentStats(al)

    rpos = al.reference_start
    qpos = 0
    qseq = al.query_alignment_sequence
    hcRegionsIx = 0
    for op, oplen in al.cigartuples:
        assert op != MATCH, """BAM must be aligned with --eqx, "M" CIGAR operation not supported"""

        if op in REFOPS:
            # record reference-consuming operations
            for i in range(oplen):
                # only consider operations within the high-confidence regions
                while hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromEnd <= rpos:
                    hcRegionsIx += 1
                if hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromStart <= rpos:
                    if op == EQUAL:
                        if rpos in hcVariants:
                            stats.ignoredBp += 1
                        else:
                            stats.matchBp += 1
                    elif op == DIFF:
                        if rpos in hcVariants:
                            stats.ignoredBp += 1
                        else:
                           stats.mismatchBp += 1
                    elif op == DEL:
                        if rpos in hcVariants or (rpos-1) in hcVariants or (rpos-i) in hcVariants or (rpos-i-1) in hcVariants:
                            stats.ignoredBp += 1
                        else:
                            curRefBp = reffasta.fetch(al.reference_name, rpos, rpos+1).upper()
                            prevRefBp = reffasta.fetch(al.reference_name, rpos-1, rpos).upper()
                            nextRefBp = reffasta.fetch(al.reference_name, rpos+1, rpos+2).upper()
                            if curRefBp in (prevRefBp, nextRefBp):
                                stats.hpDeletionBp += 1
                            else:
                                stats.nonHpDeletionBp += 1

                rpos += 1
                qpos += 1 if op in QUERYOPS else 0
        elif op == INS:
            # Only consider insertions within the high-confidence regions.
            while hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromEnd <= rpos:
                hcRegionsIx += 1
            if hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromStart <= rpos:
                if rpos in hcVariants or (rpos-1) in hcVariants:
                    stats.ignoredBp += oplen
                else:
                    curRefBp = reffasta.fetch(al.reference_name, rpos, rpos+1).upper()
                    nextRefBp = reffasta.fetch(al.reference_name, rpos+1, rpos+2).upper()
                    if qseq[qpos:qpos+oplen] in (oplen*curRefBp, oplen*nextRefBp):
                        stats.hpInsertionBp += oplen
                    else:
                        stats.nonHpInsertionBp += oplen

            qpos += oplen

    return stats


class Bed:
    def __init__(self, line):
        cols = line.rstrip("\n").split()
        self.chrom = cols[0]
        self.chromStart = int(cols[1])
        self.chromEnd = int(cols[2])

class Vcf:
    def __init__(self, chrom, chromStart, chromEnd, ref, alt):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.ref = ref
        self.alt = alt

    @staticmethod
    def LineToVcfs(line):
        cols = line.rstrip("\n").split()

        chrom = cols[0]
        chromStart = int(cols[1]) - 1
        ref = cols[3]
        chromEnd = chromStart + len(ref)
        result = []
        for alt in cols[4].split(","):
            result.append(Vcf(chrom, chromStart, chromEnd, ref, alt))

        return result

class DevNullFile:
    def __init__(self):
        pass

    def close(self):
        pass

    def write(self, s):
        pass


def main():
    args = argparser.parse_args()
    global collections; import collections
    global math; import math
    global pysam; import pysam

    # Open files
    inbam = pysam.AlignmentFile(args.inbam, check_sq=False)
    reffasta = pysam.FastaFile(args.reffasta)
    outcsv = open(args.outcsv, "w")
    outcsv.write("#read,readLengthBp,effectiveCoverage,subreadPasses,predictedConcordance,alignmentType," +
                 "alignmentMapq,hcReadLengthBp,concordance,concordanceQv,mismatchBp,nonHpInsertionBp," +
                 "nonHpDeletionBp,hpInsertionBp,hpDeletionBp\n")
    hcregionsbed = pysam.TabixFile(args.hcregions) if args.hcregions else None
    hcrcontigs = set(hcregionsbed.contigs if hcregionsbed else [])
    hcvariantsvcf = pysam.TabixFile(args.hcvariants) if args.hcvariants else None
    hcvcontigs = set(hcvariantsvcf.contigs if hcvariantsvcf else [])
    outbam = pysam.AlignmentFile(args.outbam, "wb", template=inbam) if args.outbam else None

    activeChrom,activeChromLength,activeHcRegions,activeHcVariants = None, None, collections.deque([]), collections.defaultdict(list)
    # Process alignments in sort order.
    for al in (inbam if not args.chrom else inbam.fetch(args.chrom)):
        if not al.is_unmapped:
            # For a new chromosome, update the list of HC regions and HC variants.
            if al.reference_name != activeChrom:
                activeChrom = al.reference_name
                activeChromLength = reffasta.get_reference_length(activeChrom)

                if hcregionsbed:
                    # List of non-overlapping, high-confidence regions, sorted by chromStart.
                    if al.reference_name in hcrcontigs:
                        activeHcRegions = collections.deque([Bed(x) for x in hcregionsbed.fetch(al.reference_name)])
                    else:
                        activeHcRegions = collections.deque([])
                else: # treat the whole chromosome as a high-confidence region
                    activeHcRegions = collections.deque([Bed("%s\t0\t%d" % (activeChrom, activeChromLength))])

                activeHcVariants = collections.defaultdict(list)
                if hcvariantsvcf and al.reference_name in hcvcontigs:
                    # Map from chromStart to known variants.
                    for line in hcvariantsvcf.fetch(al.reference_name):
                        for vcf in Vcf.LineToVcfs(line):
                            activeHcVariants[vcf.chromStart].append(vcf)

            # Pop high confidence regions that end before the current alignment.
            while len(activeHcRegions) and activeHcRegions[0].chromEnd <= al.reference_start:
                activeHcRegions.popleft()

            # Only measure concordance for non-secondary alignments that overlap a high confidence region.
            if (not al.is_secondary and
                len(activeHcRegions) and
                max(activeHcRegions[0].chromStart, al.reference_start) < min(activeHcRegions[0].chromEnd, al.reference_end)):

                alStats = measureConcordance(al, reffasta, activeHcRegions, activeHcVariants)
                effectiveCoverage = str(round(al.get_tag("ec"), 2)) if al.has_tag("ec") else ""
                subreadPasses = str(al.get_tag("np")) if al.has_tag("np") else ""
                predictedConcordance = "%0.6f" % (al.get_tag("rq")) if al.has_tag("rq") else ""
                if alStats.hcReadLength > 0:
                    row = "%s,%d,%s,%s,%s,%s,%d,%d,%0.6f,%0.2f,%d,%d,%d,%d,%d\n" % \
                              (al.query_name, al.query_length, effectiveCoverage, subreadPasses, predictedConcordance,
                              "Supplementary" if al.is_supplementary else "Primary", al.mapping_quality,
                              alStats.hcReadLength, alStats.concordance, alStats.concordanceQv, alStats.mismatchBp,
                              alStats.nonHpInsertionBp, alStats.nonHpDeletionBp, alStats.hpInsertionBp,
                              alStats.hpDeletionBp)
                    outcsv.write(row)

        if outbam:
            outbam.write(al)

    # Close input files
    inbam.close()
    reffasta.close()

    # Close output file
    if outbam:
        outbam.close()
    outcsv.close()


if __name__ == "__main__":
    main()
