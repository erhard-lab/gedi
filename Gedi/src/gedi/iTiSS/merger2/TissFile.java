package gedi.iTiSS.merger2;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class TissFile {
    private MemoryIntervalTreeStorage<TissFileEntry> entries;
    public final static String READ_COUNT_COLUMN_NAME = "Read-count";
    public final static String FOLD_CHANGE_COLUMN_NAME = "Fold-change";
    public final static String P_VALUE_COLUMN_NAME = "pValue";
    public final static String Z_SCORE_COLUMN_NAME = "zScore";

    private TissFile(MemoryIntervalTreeStorage<TissFileEntry> entries) {
        this.entries = entries;
    }

    public TsrFile reduceToTsr(final int gap) {
        List<TsrFileEntry> tsrs = new ArrayList<>();
        for (ReferenceSequence ref : entries.getReferenceSequences()) {
            TsrFileEntry lastEntry = null;
            for (ImmutableReferenceGenomicRegion<TissFileEntry> currentEntry : entries.ei(ref).loop()) {
                TissFileEntry tissFileEntry = currentEntry.getData();
                if (lastEntry == null) {
                    lastEntry = new TsrFileEntry(ref);
                } else if (tissFileEntry.getPosition() - lastEntry.getOriginalEnd() >= gap){
                    tsrs.add(lastEntry);
                    lastEntry = new TsrFileEntry(ref);
                }
                lastEntry.addTissFileEntry(tissFileEntry);
            }
            if (lastEntry != null) {
                tsrs.add(lastEntry);
            }
        }
        MemoryIntervalTreeStorage<TsrFileEntry> mem = new MemoryIntervalTreeStorage<>(TsrFileEntry.class);
        mem.fill(EI.wrap(tsrs).map(t -> new ImmutableReferenceGenomicRegion<>(t.getReference(), t.getOriginalRegion(), t)));
        return new TsrFile(mem);
    }

    public static TissFile loadFromFile(String path) {
        MemoryIntervalTreeStorage<TissFileEntry> entries = new MemoryIntervalTreeStorage<>(TissFileEntry.class);
        try {
            String[] header = StringUtils.split(EI.lines(path).next(), "\t");
            int readcountcolumn = ArrayUtils.find(header, READ_COUNT_COLUMN_NAME);
            int foldchangecolumn = ArrayUtils.find(header, FOLD_CHANGE_COLUMN_NAME);
            int zscorecolumn = ArrayUtils.find(header, Z_SCORE_COLUMN_NAME);
            int pvaluecolumn = ArrayUtils.find(header, P_VALUE_COLUMN_NAME);
            entries.fill(EI.lines(path).skip(1).map(l -> {
                String[] split = StringUtils.split(l, "\t");
                double readcount = 0;
                double value = 0;
                if (readcountcolumn > 1) {
                    readcount = Double.parseDouble(split[readcountcolumn]);
                }
                if (foldchangecolumn > 1) {
                    value = Double.parseDouble(split[foldchangecolumn]);
                } else if (zscorecolumn > 1) {
                    value = Double.parseDouble(split[zscorecolumn]);
                } else if (pvaluecolumn > 1) {
                    value = Double.parseDouble(split[pvaluecolumn]);
                }
                int pos = Integer.parseInt(split[1]);
                ReferenceSequence ref = Chromosome.obtain(split[0]);
                TissFileEntry entry = new TissFileEntry(ref, pos, readcount, value, Paths.get(path).getFileName().toString().hashCode());
                return new ImmutableReferenceGenomicRegion<>(ref, new ArrayGenomicRegion(pos, pos + 1), entry);
            }));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new TissFile(entries);
    }

    public ExtendedIterator<ImmutableReferenceGenomicRegion<TissFileEntry>> ei() {
        return entries.ei();
    }

    public ExtendedIterator<ImmutableReferenceGenomicRegion<TissFileEntry>> ei(ReferenceSequence ref) {
        return entries.ei(ref);
    }

    public ExtendedIterator<ImmutableReferenceGenomicRegion<TissFileEntry>> ei(ReferenceSequence ref, GenomicRegion region) {
        return entries.ei(ref, region);
    }
}
