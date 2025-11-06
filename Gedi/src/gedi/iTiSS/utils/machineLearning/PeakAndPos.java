package gedi.iTiSS.utils.machineLearning;

import org.apache.commons.math3.ml.clustering.Clusterable;

public class PeakAndPos implements Clusterable, Comparable<PeakAndPos>{
    private int pos;
    private double value;

    public PeakAndPos(int pos, double value) {
        this.pos = pos;
        this.value = value;
    }

    public int getPos() {
        return pos;
    }

    public double getValue() {
        return value;
    }

    @Override
    public double[] getPoint() {
        return new double[] {value};
    }

    @Override
    public int compareTo(PeakAndPos o) {
        return Double.compare(value, o.getValue());
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof PeakAndPos)) {
            return false;
        }
        PeakAndPos other = (PeakAndPos) obj;
        return pos == other.getPos();
    }
}
