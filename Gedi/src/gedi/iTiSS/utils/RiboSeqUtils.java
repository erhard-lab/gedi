package gedi.iTiSS.utils;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTree;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.genomic.Genomic;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.array.SparseMemoryFloatArray;

import java.util.*;

public class RiboSeqUtils {
    public static double getBiggestSingleAaExpression(ReferenceGenomicRegion region, CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codonsCit, int riboCondition, boolean useOnlyStartCodons, boolean useAllStartCodons, Genomic genomic) {
        Map<GenomicRegion, Float> regionToExpressionMap = new HashMap<>();
        for (int i = 0; i < region.getRegion().getTotalLength() - 2; i++) {
            GenomicRegion codon = region.getRegion().map(new ArrayGenomicRegion(i, i+3));
            List<ImmutableReferenceGenomicRegion<SparseMemoryFloatArray>> lst = codonsCit.ei(region.getReference(), codon).filter(f -> {
                return f.getRegion().getStart() == codon.getStart() && f.getRegion().getEnd() == codon.getEnd();
            }).filter(f -> {
                if (!useOnlyStartCodons) {
                    return true;
                }
                return SequenceUtils2.isStartCodon(genomic.getSequence(f).toString(), useAllStartCodons);
            }).list();
            if (lst.size() > 1) {
                System.err.println("Error: Multiple codons with the same position in CIT.");
            }
            if (lst.size() == 1) {
                if (regionToExpressionMap.containsKey(codon)) {
                    System.err.println("Error: Multiple codons with the same position read twice.");
                }
                regionToExpressionMap.put(codon, lst.get(0).getData().getFloat(riboCondition));
            }
        }
        if (regionToExpressionMap.values().size() == 0) {
            return 0;
        }
        return Collections.max(regionToExpressionMap.values());
    }

    public static double getBiggestSingleAaExpression(ReferenceGenomicRegion region, CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codonsCit, int riboCondition, Genomic genomic, String... codonsToCheck) {
        Map<GenomicRegion, Float> regionToExpressionMap = new HashMap<>();
        for (int i = 0; i < region.getRegion().getTotalLength() - 2; i++) {
            GenomicRegion codon = region.getRegion().map(new ArrayGenomicRegion(i, i+3));
            List<ImmutableReferenceGenomicRegion<SparseMemoryFloatArray>> lst = codonsCit.ei(region.getReference(), codon).filter(f -> {
                return f.getRegion().getStart() == codon.getStart() && f.getRegion().getEnd() == codon.getEnd();
            }).filter(f -> {
                return SequenceUtils2.codonEqualsAny(genomic.getSequence(f).toString(), codonsToCheck);
            }).list();
            if (lst.size() > 1) {
                System.err.println("Error: Multiple codons with the same position in CIT.");
            }
            if (lst.size() == 1) {
                if (regionToExpressionMap.containsKey(codon)) {
                    System.err.println("Error: Multiple codons with the same position read twice.");
                }
                regionToExpressionMap.put(codon, lst.get(0).getData().getFloat(riboCondition));
            }
        }
        if (regionToExpressionMap.values().size() == 0) {
            return 0;
        }
        return Collections.max(regionToExpressionMap.values());
    }

    public static double getAaExpression(ReferenceGenomicRegion region, CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codonsCit, int riboCondition) {
        if (region.getRegion().getTotalLength() % 3 != 0) {
            throw new IllegalArgumentException("Region need to be dividable by 3 (codons).");
        }
        double totalExpression = 0;
        for (int i = 0; i < region.getRegion().getTotalLength(); i += 3) {
            GenomicRegion codon = region.getRegion().map(new ArrayGenomicRegion(i, i+3));
            List<ImmutableReferenceGenomicRegion<SparseMemoryFloatArray>> lst = codonsCit.ei(region.getReference(), codon).filter(f -> {
                return f.getRegion().getStart() == codon.getStart() && f.getRegion().getEnd() == codon.getEnd();
            }).list();
            if (lst.size() > 1) {
                System.err.println("Error: Multiple codons with the same position in CIT.");
            }
            if (lst.size() == 1) {
                totalExpression += lst.get(0).getData().getFloat(riboCondition);
            }
        }
        return totalExpression;
    }
}
