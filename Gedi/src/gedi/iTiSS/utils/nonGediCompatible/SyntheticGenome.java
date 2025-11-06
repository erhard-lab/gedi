package gedi.iTiSS.utils.nonGediCompatible;

import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Annotation;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.functions.EI;

import java.util.ArrayList;
import java.util.List;

public class SyntheticGenome {

    public static Genomic createSyntheticGenomeOnlyTranscripts(List<String> genomicName, List<Integer> genomicLength) {
        Genomic genomic = new Genomic();
        List<ReferenceGenomicRegion<Transcript>> lst = new ArrayList<>();
        for (int i = 0; i < genomicName.size(); i++) {
            createPlusAndMinusTranscripts(genomicName.get(i), genomicLength.get(i), genomic, lst);
        }
        GenomicRegionStorage<Transcript> transcripts = new MemoryIntervalTreeStorage<>(Transcript.class);
        transcripts.fill(EI.wrap(lst));
        genomic.add(new Annotation<Transcript>(Genomic.AnnotationType.Transcripts.name()).set(transcripts));
        return genomic;
    }

    private static void createPlusAndMinusTranscripts(String name, int length, Genomic genomic, List<ReferenceGenomicRegion<Transcript>> lst) {
        if (name.toLowerCase().startsWith("chr")) {
            name = name.substring(3);
        }
        SeqProvider seqProvider = new SeqProvider(name, length);
        genomic.add(seqProvider);
        ReferenceGenomicRegion<Transcript> plus = new ImmutableReferenceGenomicRegion<Transcript>(Chromosome.obtain(name + "+"),
                new ArrayGenomicRegion(0,1), new Transcript("gene","trans",-1,-1));
        ReferenceGenomicRegion<Transcript> minus = new ImmutableReferenceGenomicRegion<Transcript>(Chromosome.obtain(name + "-"),
                new ArrayGenomicRegion(0,1), new Transcript("gene","trans",-1,-1));
        lst.add(plus); lst.add(minus);
    }
}
