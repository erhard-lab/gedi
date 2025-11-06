package gedi.iTiSS.utils.machineLearning;

import gedi.util.functions.EI;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class PeakClusterer {
    private List<PeakAndPos> data;
    private List<PeakAndPos> topCluster;
    private List<PeakAndPos> bottomCluster;
    private ClusterThreshold threshold;

    public PeakClusterer(List<PeakAndPos> data) {
        this.data = data;
    }

    public void cluster() {
        KMeansPlusPlusClusterer<PeakAndPos> kmeans = new KMeansPlusPlusClusterer<>(2);
        kmeans.getRandomGenerator().setSeed(1);
        List<CentroidCluster<PeakAndPos>> cl = kmeans.cluster(data);
        List<PeakAndPos> left = cl.get(0).getPoints();
        List<PeakAndPos> right = cl.get(1).getPoints();
        double leftAvg = EI.wrap(left).mapToDouble(PeakAndPos::getValue).sum()/left.size();
        double rightAvg = EI.wrap(right).mapToDouble(PeakAndPos::getValue).sum()/right.size();
        topCluster = leftAvg > rightAvg ? left : right;
        bottomCluster = leftAvg < rightAvg ? left : right;
    }

    private void calculateThreshold() {
        double lowThresh = Collections.max(bottomCluster, Comparator.comparingDouble(PeakAndPos::getValue)).getValue();
        double highThresh = Collections.min(topCluster, Comparator.comparingDouble(PeakAndPos::getValue)).getValue();
        threshold = new ClusterThreshold(lowThresh, highThresh);
    }

    public List<PeakAndPos> getTopCluster() {
        return topCluster;
    }

    public List<PeakAndPos> getBottomCluster() {
        return bottomCluster;
    }

    public ClusterThreshold getThreshold() {
        if (threshold == null) {
            calculateThreshold();
        }
        return threshold;
    }
}
