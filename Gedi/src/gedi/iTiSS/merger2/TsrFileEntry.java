package gedi.iTiSS.merger2;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;

import java.io.IOException;
import java.text.ParseException;
import java.util.*;

public class TsrFileEntry {
    private List<TissFileEntry> tissEntries;
    private GenomicRegion originalRegion;
    private ReferenceSequence reference;

    public TsrFileEntry(ReferenceSequence reference) {
        this.reference = reference;
        tissEntries = new ArrayList<>();
    }

    public TsrFileEntry(ReferenceSequence reference, GenomicRegion region, List<TissFileEntry> entries) {
        this.tissEntries = entries;
        this.reference = reference;
        this.originalRegion = region;
    }

    public void addTissFileEntry(TissFileEntry entry) {
        if (!entry.getReference().equals(reference)) {
            throw new IllegalArgumentException("New entry needs to be of the same reference as the TRSFileEntry. This: " + reference + ", new Entry: " + entry.getReference());
        }
        tissEntries.add(entry);
        if (originalRegion == null) {
            originalRegion = new ArrayGenomicRegion(entry.getPosition(), entry.getPosition()+1);
            return;
        }
        if (entry.getPosition() >= originalRegion.getEnd()) {
            originalRegion = originalRegion.extendBack(entry.getPosition() - originalRegion.getStop());
        } else if (entry.getPosition() < originalRegion.getStart()) {
            originalRegion = originalRegion.extendFront(originalRegion.getStart() - entry.getPosition());
        }
    }

    public String toLocationString() {
        return getReference().toPlusMinusString() + ":" + getOriginalRegion().toString2();
    }

    public ReferenceGenomicRegion<TsrFileEntry> toReferenceGenomicRegion() {
        return new ImmutableReferenceGenomicRegion<>(getReference(), getOriginalRegion(), this);
    }

    public TsrFileEntry merge(List<TsrFileEntry> others) {
        int start = EI.wrap(others).mapToInt(TsrFileEntry::getOriginalStart).min();
        int end = EI.wrap(others).mapToInt(TsrFileEntry::getOriginalEnd).max();
        start = originalRegion.getStart() < start ? originalRegion.getStart() : start;
        end = originalRegion.getEnd() > end ? originalRegion.getEnd() : end;
        GenomicRegion newRegion = new ArrayGenomicRegion(start, end);
        List<TissFileEntry> newEntries = new ArrayList<>();
        others.forEach(o -> newEntries.addAll(o.tissEntries));
        newEntries.addAll(tissEntries);
        return new TsrFileEntry(reference, newRegion, newEntries);
    }

    // TODO: should be calculated during creation. Careful with resetFileAssociations method!
    public int getScore() {
        Set<Integer> fileCodes = new HashSet<>();
        tissEntries.forEach(t -> fileCodes.add(t.getOriginId()));
        return fileCodes.size();
    }

    public int getOriginalEnd() {
        return originalRegion.getEnd();
    }

    public int getOriginalStop() {
        return originalRegion.getStop();
    }

    public int getOriginalStart() {
        return originalRegion.getStart();
    }

    public ReferenceSequence getReference() {
        return reference;
    }

    public GenomicRegion getOriginalRegion() {
        return originalRegion;
    }

    public double[] getReadCounts() {
        double[] readcounts = new double[tissEntries.size()];
        for (int i = 0; i < tissEntries.size(); i++) {
            readcounts[i] = tissEntries.get(i).getReadcount();
        }
        return readcounts;
    }

    public int computeNumberOfDisjunctPeaks() {
        int count = 1;
        List<TissFileEntry> tmp = new ArrayList<>(tissEntries);
        tmp.sort(Comparator.comparingInt(TissFileEntry::getPosition));
        int lastPosition = tmp.get(0).getPosition();
        for (int i = 1; i < tmp.size(); i++) {
            if (tmp.get(i).getPosition() - lastPosition > 1) {
                count++;
            }
            lastPosition = tmp.get(i).getPosition();
        }
        return count;
    }

    // Per non-disjunct TiSS-clusters the max TiSS is returned
    public double[] getReadCountsFromDisjunctPeaks() {
        List<TissFileEntry> biggestDisjunctTiss = new ArrayList<>();
        List<TissFileEntry> lastTiss = new ArrayList<>();
        List<TissFileEntry> tmp = new ArrayList<>(tissEntries);
        tmp.sort(Comparator.comparingInt(TissFileEntry::getPosition));
        int lastPosition = tmp.get(0).getPosition();
        lastTiss.add(tmp.get(0));
        for (int i = 1; i < tmp.size(); i++) {
            if (tmp.get(i).getPosition() - lastPosition > 1) {
                biggestDisjunctTiss.add(lastTiss.stream().max(Comparator.comparingDouble(TissFileEntry::getReadcount)).get());
                lastTiss = new ArrayList<>();
            }
            lastPosition = tmp.get(i).getPosition();
            lastTiss.add(tmp.get(i));
        }
        biggestDisjunctTiss.add(lastTiss.stream().max(Comparator.comparingDouble(TissFileEntry::getReadcount)).get());
        return biggestDisjunctTiss.stream().mapToDouble(TissFileEntry::getReadcount).toArray();
    }

    public double[] getReadCountsForId(int id) {
        List<Double> lst = new ArrayList<>();
        double[] rc = getReadCounts();
        int[] ids = getOriginIds();
        for (int i = 0; i < ids.length; i++) {
            if (ids[i] == id) {
                lst.add(rc[i]);
            }
        }
        if (lst.size() == 0) {
            return new double[0];
        }
        return lst.stream().mapToDouble(d -> d).toArray();
    }

    public double getMaxValueForId(int id) {
        List<Double> lst = new ArrayList<>();
        double[] vals = getValues();
        int[] ids = getOriginIds();
        for (int i = 0; i < ids.length; i++) {
            if (ids[i] == id) {
                lst.add(vals[i]);
            }
        }
        if (lst.size() == 0) {
            return 0;
        }
        return lst.stream().mapToDouble(d -> d).max().getAsDouble();
    }

    public double getMaxReadCountForId(int id) {
        List<Double> lst = new ArrayList<>();
        double[] rc = getReadCounts();
        int[] ids = getOriginIds();
        for (int i = 0; i < ids.length; i++) {
            if (ids[i] == id) {
                lst.add(rc[i]);
            }
        }
        if (lst.size() == 0) {
            return 0;
        }
        return lst.stream().mapToDouble(d -> d).max().getAsDouble();
    }

    public double getMaxReadCount() {
        return EI.wrap(tissEntries).mapToDouble(TissFileEntry::getReadcount).max();
    }

    public int getMaxReadCountTissPos() {
        TissFileEntry maxTiss = tissEntries.get(0);
        for (int i = 1; i < tissEntries.size(); i++) {
            maxTiss = tissEntries.get(i).getReadcount() > maxTiss.getReadcount() ? tissEntries.get(i) : maxTiss;
        }
        return maxTiss.getPosition();
    }

    public int getMaxReadCountTissPosForId(int id) {
        TissFileEntry maxTiss = null;
        for (int i = 0; i < tissEntries.size(); i++) {
            if (tissEntries.get(i).getOriginId() != id) {
                continue;
            }
            if (maxTiss == null) {
                maxTiss = tissEntries.get(i);
                continue;
            }
            maxTiss = tissEntries.get(i).getReadcount() > maxTiss.getReadcount() ? tissEntries.get(i) : maxTiss;
        }
        return maxTiss == null ? -1 : maxTiss.getPosition();
    }

    public boolean isValidatedByFile(String fileName) {
        int hash = fileName.hashCode();
        return Arrays.stream(getOriginIds()).anyMatch(x -> x == hash);
    }

    public double[] getValues() {
        double[] values = new double[tissEntries.size()];
        for (int i = 0; i < tissEntries.size(); i++) {
            values[i] = tissEntries.get(i).getValue();
        }
        return values;
    }

    public double getMaxValue() {
        double[] values = new double[tissEntries.size()];
        for (int i = 0; i < tissEntries.size(); i++) {
            values[i] = tissEntries.get(i).getValue();
        }
        if (values.length == 0) {
            return 0;
        }
        return Arrays.stream(values).max().getAsDouble();
    }

    public int[] getOriginIds() {
        int[] originIds = new int[tissEntries.size()];
        for (int i = 0; i < tissEntries.size(); i++) {
            originIds[i] = tissEntries.get(i).getOriginId();
        }
        return originIds;
    }

    public int[] getPositions() {
        int[] readcounts = new int[tissEntries.size()];
        for (int i = 0; i < tissEntries.size(); i++) {
            readcounts[i] = tissEntries.get(i).getPosition();
        }
        return readcounts;
    }

    public void resetFileAssociations(int newFileId) {
        tissEntries.forEach(t -> t.setOriginId(newFileId));
    }

    @Override
    public String toString() {
        return reference + "\t" + originalRegion + "\t" +
                StringUtils.concat(",", getReadCounts()) + "\t" +
                StringUtils.concat(",", getValues()) + "\t" +
                StringUtils.concat(",", getOriginIds()) + "\t" +
                StringUtils.concat(",", getPositions()) + "\t" +
                getScore();
    }

    @Override
    public int hashCode() {
        return (getReference().hashCode() * 59 + getOriginalRegion().hashCode2()) * 59 + tissEntries.size();
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof TsrFileEntry)) {
            return false;
        }
        TsrFileEntry entry = (TsrFileEntry) obj;
        if (!entry.getReference().equals(getReference())) {
            return false;
        }
        if (!entry.getOriginalRegion().equals2(getOriginalRegion())) {
            return false;
        }
        for (TissFileEntry tiss : tissEntries) {
            if (!entry.tissEntries.contains(tiss)) {
                return false;
            }
        }
        if (tissEntries.size() != entry.tissEntries.size()) {
            return false;
        }
        return true;
    }

    public static TsrFileEntry parse(String line) throws IOException {
        String[] split = StringUtils.split(line, "\t");
        ReferenceSequence ref = Chromosome.obtain(split[0]);
        GenomicRegion reg = GenomicRegion.parse(split[1]);
        double[] readcounts = EI.wrap(StringUtils.split(split[2],",")).mapToDouble(Double::parseDouble).toDoubleArray();
        double[] values = EI.wrap(StringUtils.split(split[3],",")).mapToDouble(Double::parseDouble).toDoubleArray();
        int[] originIds = EI.wrap(StringUtils.split(split[4],",")).mapToInt(Integer::parseInt).toIntArray();
        int[] positions = EI.wrap(StringUtils.split(split[5],",")).mapToInt(Integer::parseInt).toIntArray();
        if (readcounts.length != values.length || values.length != originIds.length || originIds.length != positions.length) {
            throw new IOException(line);
        }
        List<TissFileEntry> entries = new ArrayList<>();
        for (int i = 0; i < readcounts.length; i++) {
            entries.add(new TissFileEntry(ref, positions[i], readcounts[i], values[i], originIds[i]));
        }
        return new TsrFileEntry(ref, reg, entries);
    }
}
