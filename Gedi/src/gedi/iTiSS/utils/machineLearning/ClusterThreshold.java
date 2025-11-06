package gedi.iTiSS.utils.machineLearning;

import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTuple;

public class ClusterThreshold extends MutablePair<Double, Double> {
    public ClusterThreshold(Double lowThresh, Double highThresh) {
        super(lowThresh, highThresh);
    }

    public double getLowThresh() {
        return this.Item1;
    }

    public double getHighThresh() {
        return this.Item2;
    }

    public double getMeanThresh() {
        return (this.Item1 + this.Item2) / 2.0;
    }
}
