package gedi.iTiSS.utils;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.*;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.StringUtils;
import gedi.util.functions.EI;

import java.util.*;

public class GenomicUtils {
    public static MemoryIntervalTreeStorage<String> getOverlappingGenes(Genomic genomic) {
        MemoryIntervalTreeStorage<String> genes = genomic.getGenes();
        Map<ReferenceSequence, List<ImmutableReferenceGenomicRegion<String>>> overlappingGenes = new HashMap<>();
        genes.ei().forEachRemaining(g -> {
            // refGenes is filled with all the genes for the current reference in this loop
            List<ImmutableReferenceGenomicRegion<String>> refGenes = overlappingGenes.computeIfAbsent(g.getReference(), absent -> new ArrayList<>());
            for (ImmutableReferenceGenomicRegion<String> gene : refGenes) {
                // while filling refGenes, every new gene is checked whether it intersects with one of the genes in refGenes
                if (gene.getRegion().intersects(g.getRegion())) {
                    // if this is the case, the two intersecting genes are merged by unionising their GenomicRegion
                    refGenes.remove(gene);
                    refGenes.add(g.toMutable().setRegion(g.getRegion().union(gene.getRegion())).setData(gene.getData() + "_" + g.getData()).toImmutable());
                    return;
                }
            }
            refGenes.add(g);
        });
        // We create a MemoryIntervalTreeStorage, which is filled with the unionized genes
        MemoryIntervalTreeStorage<String> out = new MemoryIntervalTreeStorage<>(String.class);
        out.fill(EI.wrap(overlappingGenes.values()).unfold(EI::wrap));
        return out;
    }

    public static MemoryIntervalTreeStorage<String> getUniqueTss(Genomic genomic, int mergingWindow, boolean onlyProteinCoding) {
        MemoryIntervalTreeStorage<String> genes = getOverlappingGenes(genomic);
        MemoryIntervalTreeStorage<Transcript> transcripts = genomic.getTranscripts();
        MemoryIntervalTreeStorage<String> out = new MemoryIntervalTreeStorage<>(String.class);

        out.fill(genes.ei().unfold(g -> {
            // we unfold every gene into the TiSS for each of its transcripts
            ArrayList<ImmutableReferenceGenomicRegion<String>> extendedTssLst = transcripts.ei(g).filter(f -> !onlyProteinCoding || f.getData().isCoding()).map(m -> {
                // all TiSS are extended by mergingWindow
                int tss = GenomicRegionPosition.FivePrime.position(m);
                return new ImmutableReferenceGenomicRegion<>(m.getReference(), new ArrayGenomicRegion(tss - mergingWindow, tss + mergingWindow + 1), g.getData());
            }).list();
            if (extendedTssLst.size() == 0) {
                return EI.wrap(extendedTssLst);
            }
            extendedTssLst.sort(Comparator.comparingInt(GenomicRegionPosition.FivePrime::position));
            ImmutableReferenceGenomicRegion<String> last = extendedTssLst.get(0);
            ArrayList<ImmutableReferenceGenomicRegion<String>> mergedTssLst = new ArrayList<>();
            for (int i = 1; i < extendedTssLst.size(); i++) {
                // We go through the sorted list of TiSS and merge intersecting ones similar as we did in getOverlappingGenes()
                if (last.getRegion().intersects(extendedTssLst.get(i).getRegion())) {
                    last = new ImmutableReferenceGenomicRegion<>(last.getReference(), new ArrayGenomicRegion(last.getRegion().getStart(), extendedTssLst.get(i).getRegion().getEnd()), last.getData());
                } else {
                    mergedTssLst.add(last);
                    last = extendedTssLst.get(i);
                }
            }
            mergedTssLst.add(last);
            return EI.wrap(mergedTssLst);
        }).removeNulls());
        return out;
    }

    public static String getGeneNameFromGeneId(String geneId, Genomic genomic) {
        return genomic.getGeneTable("symbol").apply(geneId);
    }

    public static ReferenceGenomicRegion<NameAnnotation> getRgrFromGeneName(String geneName, Genomic genomic) {
        ReferenceGenomicRegion<?> rgr = genomic.getNameIndex().get(geneName);
        if (rgr == null) {
            return null;
        }
        return new ImmutableReferenceGenomicRegion<>(rgr.getReference(), rgr.getRegion(), new NameAnnotation(rgr.getData().toString()));
    }

    public static ReferenceGenomicRegion<String> getRgrFromGeneId(String geneId, Genomic genomic) {
        return genomic.getGeneMapping().apply(geneId);
    }

    public static List<ImmutableReferenceGenomicRegion<Transcript>> getTranscripts(ReferenceGenomicRegion<String> geneRgr, Genomic genomic) {
        return genomic.getTranscripts().ei(geneRgr).filter(t -> t.getData().getGeneId().equals(geneRgr.getData())).list();
    }

    public static List<ImmutableReferenceGenomicRegion<Transcript>> getTranscripts(ReferenceGenomicRegion<String> geneRgr, MemoryIntervalTreeStorage<Transcript> transcripts) {
        return transcripts.ei(geneRgr).filter(t -> t.getData().getGeneId().equals(geneRgr.getData())).list();
    }

    public static boolean isProteinCoding(ReferenceGenomicRegion<String> gene, Genomic genomic) {
        if (gene == null) {
            return false;
        }
        for (ImmutableReferenceGenomicRegion<Transcript> transcript : getTranscripts(gene, genomic)) {
            if (transcript.getData().isCoding()) {
                return true;
            }
        }
        return false;
    }

    public static boolean isProteinCoding(ReferenceGenomicRegion<String> gene, MemoryIntervalTreeStorage<Transcript> transcripts) {
        if (gene == null) {
            return false;
        }
        for (ImmutableReferenceGenomicRegion<Transcript> transcript : getTranscripts(gene, transcripts)) {
            if (transcript.getData().isCoding()) {
                return true;
            }
        }
        return false;
    }

    public static String toChromosomeAddedName(ReferenceSequence ref) {
        try {
            int i = Integer.parseInt(ref.getName());
            return "chr" + ref.getName();
        } catch (NumberFormatException e) {
            if (ref.getName().equals("X")) {
                return "chrX";
            }
            if (ref.getName().equals("Y")) {
                return "chrY";
            }
            return ref.getName();
        }
    }

    public static boolean isStandardChromosome(ReferenceSequence ref) {
        try {
            int i = Integer.parseInt(ref.getName());
            return true;
        } catch (NumberFormatException e) {
            if (ref.getName().equals("X")) {
                return true;
            }
            if (ref.getName().equals("Y")) {
                return true;
            }
            return false;
        }
    }

    // not very performant. Use only for small genomes or implement using suffix array BWT thingy
    public static MemoryIntervalTreeStorage<String> getSequencePositions(Genomic genomic, final String querySeq) {
        MemoryIntervalTreeStorage<String> out = new MemoryIntervalTreeStorage<>(String.class);
        out.fill(genomic.iterateReferenceSequences().unfold(ref -> {
            int refLength =  genomic.getLength(ref.toPlusMinusString());
            CharSequence sequence = genomic.getSequence(ref, new ArrayGenomicRegion(0, refLength));
            List<ImmutableReferenceGenomicRegion<String>> lst = new ArrayList<>();
            for (int i = 0; i < sequence.length() - (querySeq.length() - 1); i++) {
                if (sequence.subSequence(i, i + querySeq.length()).toString().equals(querySeq)) {
                    GenomicRegion reg;
                    if (ref.isPlus()) {
                        reg = new ArrayGenomicRegion(i, i + querySeq.length());
                    } else {
                        reg = new ArrayGenomicRegion(refLength-i-querySeq.length(),refLength-i);
                    }
                    lst.add(new ImmutableReferenceGenomicRegion<>(ref, reg, querySeq));
                }
            }
            return EI.wrap(lst);
        }));
        return out;
    }
}
