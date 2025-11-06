package gedi.iTiSS.merger2;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.StringUtils;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

import java.io.IOException;
import java.text.ParseException;
import java.util.*;

public class TsrFile {
    private MemoryIntervalTreeStorage<TsrFileEntry> entries;
    private String[] originNames;
    private Map<Integer, String> hashCodeToOriginNames;

    protected TsrFile(MemoryIntervalTreeStorage<TsrFileEntry> entries) {
        this.entries = entries;
    }

    public TsrFile filter(int minScore) {
        MemoryIntervalTreeStorage<TsrFileEntry> filtered = new MemoryIntervalTreeStorage<>(TsrFileEntry.class);
        filtered.fill(entries.ei().filter(e -> e.getData().getScore() >= minScore));
        return new TsrFile(filtered);
    }

    public Iterable<TsrFileEntry> loop() {
        return entries.ei().map(e -> e.getData()).loop();
    }

    public ExtendedIterator<TsrFileEntry> ei() {
        return entries.ei().map(e -> e.getData());
    }

    public ExtendedIterator<TsrFileEntry> ei(ReferenceSequence ref, GenomicRegion reg) {
        return entries.ei(ref, reg).map(e -> e.getData());
    }

    public ExtendedIterator<TsrFileEntry> ei(ReferenceGenomicRegion rgr) {
        return ei(rgr.getReference(), rgr.getRegion());
    }

    public ExtendedIterator<TsrFileEntry> ei(ReferenceSequence ref) {
        return entries.ei(ref).map(e -> e.getData());
    }

    public TsrFile merge(TsrFile other, int extension) {
        Map<ReferenceSequence, IntervalTree<GenomicRegion, TsrFileEntry>> tmpmem = new HashMap<>();
        entries.getReferenceSequences().forEach(ref -> tmpmem.put(ref, new IntervalTree<>(ref)));
        entries.ei().forEachRemaining(entry -> tmpmem.get(entry.getReference()).put(entry.getRegion(), entry.getData()));
        for (TsrFileEntry entry : other.loop()) {
            if (!tmpmem.containsKey(entry.getReference())) {
                tmpmem.put(entry.getReference(), new IntervalTree<>(entry.getReference()));
            }
            List<TsrFileEntry> overlaps = EI.wrap(tmpmem.get(entry.getReference()).iterateIntervalsIntersecting(entry.getOriginalStart()-extension, entry.getOriginalStop()+extension, i->true)).map(Map.Entry::getValue).list();
//            List<TsrFileEntry> overlaps = tmpmem.get(entry.getReference()).ei(entry.getOriginalStart()-extension, entry.getOriginalEnd()+extension).map(e -> e.getValue()).list();
            if (overlaps.size() == 0) {
                tmpmem.get(entry.getReference()).put(entry.getOriginalRegion(), entry);
            } else {
                for (TsrFileEntry e : overlaps) {
                    tmpmem.get(e.getReference()).remove(e.getOriginalRegion());
                }
                TsrFileEntry mergedEntry = entry.merge(overlaps);
                tmpmem.get(mergedEntry.getReference()).put(mergedEntry.getOriginalRegion(), mergedEntry);
            }
        }
        MemoryIntervalTreeStorage<TsrFileEntry> mem = new MemoryIntervalTreeStorage<>(TsrFileEntry.class);
        Set<ImmutableReferenceGenomicRegion<TsrFileEntry>> set = new HashSet<>();
        EI.wrap(tmpmem.keySet()).forEachRemaining(ref -> {
            set.addAll(tmpmem.get(ref).ei().set());
        });
        mem.fill(EI.wrap(set));
        return new TsrFile(mem);
    }

    public void writeToFile(String path, String[] originNames) throws IOException {
        if (!path.endsWith(".tsr")) {
            path += ".tsr";
        }
        LineWriter writer = new LineOrientedFile(path).write();
        writer.writeLine("# " + StringUtils.concat(",", originNames));
        writer.writeLine("Reference\tRegion\tReadcounts\tValues\tOriginIds\tPositions\tScore");
        for (TsrFileEntry entry : entries.ei().map(e -> e.getData()).loop()) {
            writer.writeLine(entry.toString());
        }
        writer.close();
    }

    public static TsrFile loadFromFile(String path) throws IOException {
        LineIterator lines = EI.lines(path);
        String comment = lines.next();
        String[] originNames = null;
        if (comment.startsWith("# ")) {
            originNames = StringUtils.split(StringUtils.splitField(comment, " ", 1), ",");
            lines.next();
        }
        MemoryIntervalTreeStorage<TsrFileEntry> entries = new MemoryIntervalTreeStorage<>(TsrFileEntry.class);
        entries.fill(lines.map(l -> {
            try {
                TsrFileEntry entry = TsrFileEntry.parse(l);
                return new ImmutableReferenceGenomicRegion<>(entry.getReference(), entry.getOriginalRegion(), entry);
            } catch (IOException e) {
                e.printStackTrace();
                return null;
            }
        }).removeNulls());
        TsrFile file = new TsrFile(entries);
        file.originNames = originNames;
        return file;
    }

    public void resetFileAssociations(int newFileId) {
        originNames = null;
        entries.ei().forEachRemaining(e -> e.getData().resetFileAssociations(newFileId));
    }

    public String[] getOriginNames() {
        return originNames;
    }

    private void createHashCodeToOriginNamesMap() {
        hashCodeToOriginNames = new HashMap<>();
        for (String originName : getOriginNames()) {
            hashCodeToOriginNames.put(originName.hashCode(), originName);
        }
    }

    public Set<String> getOriginFiles(TsrFileEntry tsrEntry) {
        if (hashCodeToOriginNames == null) {
            createHashCodeToOriginNamesMap();
        }
        Set<String> originFiles = new HashSet<>();
        for (int originId : tsrEntry.getOriginIds()) {
            originFiles.add(hashCodeToOriginNames.get(originId));
        }
        return originFiles;
    }
}
