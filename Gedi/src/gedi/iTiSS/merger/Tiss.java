package gedi.iTiSS.merger;

import java.util.List;
import java.util.Set;

public class Tiss {
    private Set<Integer> calledBy;
    private int tissPos;

    public Tiss(Set<Integer> calledBy, int tissPos) {
        this.calledBy = calledBy;
        this.tissPos = tissPos;
    }

    public Set<Integer> getCalledBy() {
        return calledBy;
    }

    public int getTissPos() {
        return tissPos;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Integer) {
            Integer other = (Integer)obj;
            return tissPos == other;
        }
        if (!(obj instanceof Tiss)) {
            return false;
        }
        Tiss other = (Tiss)obj;
        return tissPos == other.getTissPos();
    }

    @Override
    public int hashCode() {
        return new Integer(tissPos).hashCode();
    }

    public static int compare(Tiss t1, Tiss t2) {
        return Integer.compare(t1.tissPos, t2.tissPos);
    }
}
