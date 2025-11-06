package gedi.iTiSS.merger;

import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.util.functions.EI;

import java.util.*;

public class Tsr {
    private Set<Integer> calledBy;
    private GenomicRegion region;
    private List<Tiss> maxTiss;

    public Tsr(Tiss tiss) {
        calledBy = new HashSet<>();
        calledBy.addAll(tiss.getCalledBy());
        region = new ArrayGenomicRegion(tiss.getTissPos(), tiss.getTissPos()+1);
        maxTiss = new ArrayList<>();
        maxTiss.add(tiss);
    }

    public Tsr(GenomicRegion region) {
        calledBy = new HashSet<>();
        maxTiss = new ArrayList<>();
        this.region = region;
    }

    public Tsr(int tiss) {
        this(new Tiss(new HashSet<>(), tiss));
    }

    public void add(Tiss tiss) {
        if (tiss.getTissPos() > region.getStop()) {
            region = region.extendBack(tiss.getTissPos()-region.getStop());
        } else if (tiss.getTissPos() < region.getStart()) {
            region = region.extendFront(region.getStart()-tiss.getTissPos());
        }
        calledBy.addAll(tiss.getCalledBy());
        addMaxTiss(tiss);
    }

    public void addAll(Collection<Tiss> tissLst) {
        for (Tiss tss : tissLst) {
            add(tss);
        }
    }

    public void add(int tiss) {
        add(new Tiss(new HashSet<>(), tiss));
    }

    private void addMaxTiss(Tiss tiss) {
        if (maxTiss.contains(tiss)) {
            maxTiss.get(maxTiss.indexOf(tiss)).getCalledBy().addAll(tiss.getCalledBy());
        } else {
            maxTiss.add(tiss);
        }
    }

    public void extendFron(int delta) {
        region = region.extendFront(delta);
    }

    public void extendBack(int delta) {
        region = region.extendBack(delta);
    }

    public int getEnd() {
        return region.getEnd();
    }

    public int getStart() {
        return region.getStart();
    }

    public int getRightmostTissPosition() {
        int max = -1;
        for (Tiss tiss : getMaxTiss()) {
            if (tiss.getTissPos() > max) {
                max = tiss.getTissPos();
            }
        }
        if (max == -1) {
            System.err.println("Careful, no TiSS inside the TSR. This should never happen.");
        }
        return max;
    }

    public int getLeftmostTissPosition() {
        int min = Integer.MAX_VALUE;
        for (Tiss tiss : getMaxTiss()) {
            if (tiss.getTissPos() < min) {
                min = tiss.getTissPos();
            }
        }
        if (min == Integer.MAX_VALUE) {
            System.err.println("Careful, no TiSS inside the TSR. This should never happen.");
        }
        return min;
    }

    public Set<Integer> getCalledBy() {
        return calledBy;
    }

    public List<Tiss> getMaxTiss() {
        return maxTiss;
    }

    public boolean contains(Tsr tsr) {
        return region.intersects(new ArrayGenomicRegion(tsr.getStart(), tsr.getEnd()));
    }

    public boolean contains(int pos) {
        return region.contains(pos);
    }

    public boolean containsMaxTiSS(Tsr tsr) {
        if (!contains(tsr)) {
            return false;
        }
        for (Tiss maxTss : tsr.getMaxTiss()) {
            if (region.contains(maxTss.getTissPos())) {
                return true;
            }
        }
        return false;
    }

    public Set<Tiss> getContainingMaxTiSS(Tsr tsr) {
        Set<Tiss> contained = new HashSet<>();
        for (Tiss maxTss : tsr.getMaxTiss()) {
            if (region.contains(maxTss.getTissPos())) {
                contained.add(maxTss);
            }
        }
        return contained;
    }

    public Tiss getSingleMaxTiSS() {
        return Collections.max(getMaxTiss(), Comparator.comparingInt(a -> a.getCalledBy().size()));
    }

    public int getSingleMaxTiSSScore() {
        return Collections.max(getMaxTiss(), Comparator.comparingInt(a -> a.getCalledBy().size())).getCalledBy().size();
    }

    public List<Tsr> breakUp(int gap) {
        List<Tsr> outList = new ArrayList<>();
        if (maxTiss.size() == 0) {
            return outList;
        }
        List<Tiss> maxTissClone = new ArrayList<>(maxTiss);
        maxTissClone.sort(Tiss::compare);
        Tiss lastTiss = maxTissClone.get(0);
        Tsr tsr = new Tsr(lastTiss);
        for (int i = 1; i < maxTissClone.size(); i++) {
            if (maxTissClone.get(i).getTissPos() - lastTiss.getTissPos() <= gap) {
                tsr.add(maxTissClone.get(i));
            } else {
                tsr.extendBack(gap);
                tsr.extendFron(gap);
                outList.add(tsr);
                tsr = new Tsr(maxTissClone.get(i));
            }
            lastTiss = maxTissClone.get(i);
        }
        tsr.extendBack(gap);
        tsr.extendFron(gap);
        outList.add(tsr);
        return outList;
    }

    public void removeTiss(int minTissScore) {
        List<Tiss> tissToAdd = new ArrayList<>();
        for (Tiss tiss : maxTiss) {
            if (tiss.getCalledBy().size() >= minTissScore) {
                tissToAdd.add(tiss);
            }
        }
        calledBy = new HashSet<>();
        maxTiss = new ArrayList<>();
        tissToAdd.forEach(this::add);
    }

    public static int compare(Tsr tsr1, Tsr tsr2) {
        return Integer.compare(tsr1.getStart(), tsr2.getStart());
    }

    @Override
    public int hashCode() {
        return getCalledBy().hashCode() + getMaxTiss().hashCode() + region.hashCode2();
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof Tsr)) {
            return false;
        }
        Tsr other = (Tsr)obj;
        return maxTiss.containsAll(other.getMaxTiss()) && other.getMaxTiss().containsAll(maxTiss) && region.equals2(new ArrayGenomicRegion(other.getStart(), other.getEnd()));
    }

    @Override
    public String toString() {
        StringBuilder outString = new StringBuilder(region.getStart() + "-" + region.getEnd());
        for (Integer i : calledBy) {
            outString.append("\t").append(i);
        }
        return outString.toString();
    }
}
