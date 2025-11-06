package gedi.iTiSS.utils.loader;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.util.functions.TriConsumer;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.tsv.GenomicTsvFileReader;

import java.util.HashSet;
import java.util.Set;

public class TsrDataFileReader extends GenomicTsvFileReader<TsrData> {

    public TsrDataFileReader(String path) {
        super(path, true, "\t", new TsrDataElementParser(), null, TsrData.class);
    }

    public static class TsrDataElementParser implements TriConsumer<HeaderLine, String[], MutableReferenceGenomicRegion<TsrData>> {

        @Override
        public void accept(HeaderLine a, String[] fields, MutableReferenceGenomicRegion<TsrData> box) {
            Set<Integer> score = new HashSet<>();
            for (int i = 3; i < fields.length-2; i++) {
                int called = Integer.parseInt(fields[i]);
                if (called == 1) {
                    score.add(i-3);
                } else if (called == 0) {
                    // Not called
                } else {
                    throw new IllegalStateException();
                }
            }
            GenomicRegion reg = GenomicRegion.parse(fields[2]);
            int maxTissPos = Integer.parseInt(fields[1]);
            ReferenceSequence ref = Chromosome.obtain(fields[0]);
            box.set(ref, reg, new TsrData(score.size(), maxTissPos, score));
        }
    }
}
