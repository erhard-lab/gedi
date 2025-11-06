package gedi.iTiSS.utils;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;

import java.lang.ref.Reference;

public class GenomicRegionUtils {
    public static boolean isContainedInIntron(GenomicRegion region, int pos) {
        return pos >= region.getStart() && pos < region.getEnd() && !region.contains(pos);
    }

    public static boolean notInIntron(GenomicRegion region, int pos) {
        return region.contains(pos) || pos < region.getStart() || pos >= region.getEnd();
    }

    public static int getIntronIndex(GenomicRegion region, int pos) {
        if (!isContainedInIntron(region, pos)) {
            return -1;
        }
        return getPartIndes(region.invert(), pos);
    }

    public static int getPartIndes(GenomicRegion region, int pos) {
        if (!region.contains(pos)) {
            return -1;
        }
        for (int i = 0; i < region.getNumParts(); i++) {
            if (region.getPart(i).asRegion().contains(pos)) {
                return i;
            }
        }
        return -1;
    }

    public static GenomicRegion removeIntron(GenomicRegion region, int intronIndex) {
        if (intronIndex < 0 || intronIndex > region.invert().getNumParts()) {
            throw new IndexOutOfBoundsException("Region intron parts: " + region.invert().getNumParts() + ", index: " + intronIndex);
        }

        int[] newRegionIndices = new int[(region.getNumParts()-1)*2];
        int newRegionIndex = 0;
        for (int i = 0; i < region.getNumParts(); i++) {
            if (i == intronIndex) {
                newRegionIndices[newRegionIndex*2] = region.getPart(i).getStart();
                newRegionIndices[newRegionIndex*2+1] = region.getPart(i+1).getEnd();
                newRegionIndex++;
                i++;
                continue;
            }
            newRegionIndices[newRegionIndex*2] = region.getPart(i).getStart();
            newRegionIndices[newRegionIndex*2+1] = region.getPart(i).getEnd();
            newRegionIndex++;
        }

        return new ArrayGenomicRegion(newRegionIndices);
    }

    public static GenomicRegion createRegion(ReferenceSequence ref, int tss, int tts) {
        if (ref.isPlus()) {
            return new ArrayGenomicRegion(tss, tts + 1);
        } else {
            return new ArrayGenomicRegion(tts, tss + 1);
        }
    }

    public static String toStringReversed(GenomicRegion region) {
        StringBuilder sb = new StringBuilder();
        for (int i = region.getNumParts()-1; i >= 0; i--) {
            sb.append(region.getPart(i).getEnd()).append("-").append(region.getPart(i).getStart());
            if (i != 0) {
                sb.append("|");
            }
        }
        return sb.toString();
    }

    public static String toLocationStringReversed(ReferenceGenomicRegion rgr) {
        return rgr.getReference().toPlusMinusString() + ":" + toStringReversed(rgr.getRegion());
    }

    public static GenomicRegion move(GenomicRegion regToMove, int positionsToMove) {
        int[] newcoords = new int[regToMove.getNumParts()*2];
        for (int i = 0; i < regToMove.getNumParts(); i++) {
            newcoords[i*2] = regToMove.getStart(i) + positionsToMove;
            newcoords[i*2+1] = regToMove.getEnd(i) + positionsToMove;
        }
        return new ArrayGenomicRegion(newcoords);
    }

    public static <T> ReferenceGenomicRegion<T> moveUpstream(ReferenceGenomicRegion<T> rgr, int positionsToMove) {
        if (rgr.getReference().isPlus() || rgr.getReference().isIndependent()) {
            return new ImmutableReferenceGenomicRegion<>(rgr.getReference(), move(rgr.getRegion(), -positionsToMove), rgr.getData());
        } else {
            return new ImmutableReferenceGenomicRegion<>(rgr.getReference(), move(rgr.getRegion(), positionsToMove), rgr.getData());
        }
    }

    public static <T> ReferenceGenomicRegion<T> moveDownstream(ReferenceGenomicRegion<T> rgr, int positionsToMove) {
        if (rgr.getReference().isPlus() || rgr.getReference().isIndependent()) {
            return new ImmutableReferenceGenomicRegion<>(rgr.getReference(), move(rgr.getRegion(), positionsToMove), rgr.getData());
        } else {
            return new ImmutableReferenceGenomicRegion<>(rgr.getReference(), move(rgr.getRegion(), -positionsToMove), rgr.getData());
        }
    }
}
