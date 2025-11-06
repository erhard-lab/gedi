package gedi.iTiSS.utils.datastructures;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;

import java.io.IOException;
import java.util.*;
import java.util.function.Function;

public class TissFileOld<T> {
    private Map<ReferenceSequence, Map<Integer, TissFileEntryOld<T>>> entries;

    public TissFileOld() {
        this.entries = new HashMap<>();
    }

    public void add(TissFileEntryOld<T> toAdd) {
        entries.computeIfAbsent(toAdd.getReference(), absent -> new HashMap<>()).put(toAdd.getPosition(), toAdd);
    }

    public void addAll(Collection<TissFileEntryOld<T>> toAdd) {
        toAdd.forEach(this::add);
    }

    public void addAll(TissFileOld<T> toAdd) {
        toAdd.iterateEntries().forEachRemaining(this::add);
    }

    public boolean contains(TissFileEntryOld<T> entry) {
        Map<Integer, TissFileEntryOld<T>> map = entries.get(entry.getReference());
        if (map == null) {
            return false;
        }
        TissFileEntryOld<T> thisEntry = map.get(entry.getPosition());
        if (thisEntry == null) {
            return false;
        }
        return thisEntry.equals(entry);
    }

    public Iterator<TissFileEntryOld<T>> iterateEntries() {
        return EI.wrap(entries.values()).unfold(e -> e.values().iterator());
    }

    public Iterator<TissFileEntryOld<T>> iterateEntries(ReferenceSequence reference) {
        if (!entries.containsKey(reference)) {
            return EI.empty();
        }
        return entries.get(reference).values().iterator();
    }

    public Set<ReferenceSequence> referenceSequences() {
        return entries.keySet();
    }

    public TissFileEntryOld<T> getNearest(ReferenceSequence reference, int position) {
        if (!entries.containsKey(reference)) {
            return null;
        }
        int nearest = Integer.MIN_VALUE;
        int nearestDist = Integer.MAX_VALUE;
        for (int i : entries.get(reference).keySet()) {
            int currestDist = Math.abs(position-i);
            if (currestDist < nearestDist) {
                nearestDist = currestDist;
                nearest = i;
            }
        }

        return entries.get(reference).get(nearest);
    }

    public List<TissFileEntryOld<T>> getWithinRangeSorted(ReferenceSequence reference, int position, int range, boolean lookLeft, boolean lookRight) {
        if (!entries.containsKey(reference)) {
            return new ArrayList<>();
        }
        List<TissFileEntryOld<T>> out = new ArrayList<>();
        for (int i : entries.get(reference).keySet()) {
            if (lookLeft && lookRight) {
                if (Math.abs(i - position) <= range) {
                    out.add(entries.get(reference).get(i));
                }
            } else {
                if (lookLeft) {
                    if (position - i > 0 && position - i <= range) {
                        out.add(entries.get(reference).get(i));
                    }
                } else {
                    if (i - position > 0 && i - position <= range) {
                        out.add(entries.get(reference).get(i));
                    }
                }
            }
        }

        out.sort(Comparator.comparingInt(TissFileEntryOld::getPosition));

        return out;
    }

    public static <T> TissFileOld<T> read(String filePath, Function<String, T> dataParse) throws IOException {
        String[] firstLine = StringUtils.split(EI.lines(filePath).next(), "\t");
        if (firstLine.length < 2) {
            throw new IOException("Not a valid TiSS file from iTiSS");
        }
        if (!firstLine[0].equals("Reference (Strand)") || !firstLine[1].equals("Position")) {
            if (firstLine[1].contains("-")) {
                GenomicRegion reg = GenomicRegion.parse(firstLine[1]);
                if (reg.getTotalLength() != 1) {
                    throw new IOException("TiSS regions only supported with a length of 1");
                }
            } else {
                if (!StringUtils.isInt(firstLine[1])) {
                    throw new IOException("Second column is neither a region nor a position");
                }
            }
        }

        TissFileOld<T> spf = new TissFileOld<>();

        EI.lines(filePath).skip(1).forEachRemaining(l -> {
            String[] split = StringUtils.split(l,"\t");
            String data = StringUtils.concat("\t", ArrayUtils.slice(split, 2));
            Integer pos = split[1].contains("-") ? GenomicRegion.parse(split[1]).getStart() : Integer.parseInt(split[1]);
            TissFileEntryOld<T> entry = new TissFileEntryOld<>(Chromosome.obtain(split[0]), pos, dataParse.apply(data));
            spf.add(entry);
        });

        return spf;
    }
}
