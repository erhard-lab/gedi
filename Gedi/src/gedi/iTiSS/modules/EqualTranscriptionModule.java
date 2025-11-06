package gedi.iTiSS.modules;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.utils.sortedNodesList.StaticSizeSortedArrayList;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

public class EqualTranscriptionModule extends ModuleBase {
    private int window = 500;
    private int tissMaskWindow = 40;
    private double threshold = 0.001;
    private boolean useML = false;
    private String prefix;
    private double minDownstreamReadMean = 1.0;
    private double upstreamDownstreamEqualityThresh = 0.05;

    /**
     * Uses a kolmogorovSmirnovTest to check if the read start distribution upstream does not follow a normal
     * distribution. Further checks are:
     * - upstream and downstream read start distributions are significantly different
     * - upstream distribution is significantly less likely a normal distribution than the downstream distribution.
     * - downstream mean is greater than the upstream mean
     * @param moduleName
     * @param lane
     * @param window
     * @param tissMaskWindow
     * @param threshold
     * @param useML
     * @param prefix
     * @param minDownstreamReadMean
     * @param upstreamDownstreamEqualityThresh
     */
    public EqualTranscriptionModule(String moduleName, Data lane, int window, int tissMaskWindow, double threshold, boolean useML, String prefix, double minDownstreamReadMean, double upstreamDownstreamEqualityThresh) {
        super(moduleName, lane);
        this.window = window;
        this.tissMaskWindow = tissMaskWindow;
        this.threshold = threshold;
        this.useML = useML;
        this.prefix = prefix;
        this.minDownstreamReadMean = minDownstreamReadMean;
        this.upstreamDownstreamEqualityThresh = upstreamDownstreamEqualityThresh;
    }

    public EqualTranscriptionModule(String moduleName, Data lane) {
        super(moduleName, lane);
    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        NumericArray reads = data[0];

        double[] upstreamWindow = new double[window];
        double[] downstreamWindow = new double[window];
        // TODO use maskWindow size
        double[] tissMaskMiddleWindow = new double[window];

        for (int i = 0; i < window; i++) {
            upstreamWindow[i] = reads.getDouble(i);
        }
        // TODO is this correct?
        for (int i = window; i < window+tissMaskWindow; i++) {
            tissMaskMiddleWindow[i-window] = reads.getDouble(i);
        }
        for (int i = window + tissMaskWindow; i < window*2+tissMaskWindow; i++) {
            downstreamWindow[i-(window+tissMaskWindow)] = reads.getDouble(i);
        }
        if (ref.isMinus()) {
            double[] tmp = upstreamWindow;
            upstreamWindow = downstreamWindow;
            downstreamWindow = tmp;
        }

        double upstreamMean = calculateMean(upstreamWindow);
        double downstreamMean = calculateMean(downstreamWindow);
        StaticSizeSortedArrayList<Double> upstreamWindowList = new StaticSizeSortedArrayList<>(EI.wrap(upstreamWindow).list(), Comparator.comparingDouble(d -> -d));
        StaticSizeSortedArrayList<Double> downstreamWindowList = new StaticSizeSortedArrayList<>(EI.wrap(downstreamWindow).list(), Comparator.comparingDouble(d -> -d));
        StaticSizeSortedArrayList<Double> tissMaskWindowList = new StaticSizeSortedArrayList<>(EI.wrap(tissMaskMiddleWindow).list(), Comparator.comparingDouble(d -> -d));

        Map<Integer, Map<String, Double>> foundPeaks = this.res.computeIfAbsent(ref, k -> new HashMap<>());

        for (int i = 0; i < reads.length() - (window*2+tissMaskWindow)-1; i++) {
            if (i%10000000 == 0) {
                System.err.print("Progress " + i + "/" + reads.length() + ", found peaks: " + foundPeaks.keySet().size() + "\r");
            }

            // This TiSS needs to be the highest peak, otherwise this probably isn't the right TiSS
            if (downstreamWindowList.getValueAtIndex(0) > tissMaskWindowList.getValueAtIndex(0) || upstreamWindowList.getValueAtIndex(0) > tissMaskWindowList.getValueAtIndex(0)) {
                if (ref.isPlus()) {
                    updateLists(reads, upstreamWindowList, downstreamWindowList, tissMaskWindowList, i);
                    double rem = reads.getDouble(i+window+tissMaskWindow);
                    double add = reads.getDouble(i+window*2+tissMaskWindow);
                    downstreamMean -= rem/window;
                    downstreamMean += add/window;
                    rem = reads.getDouble(i);
                    add = reads.getDouble(i+window);
                    upstreamMean -= rem/window;
                    upstreamMean += add/window;
                } else {
                    updateLists(reads, downstreamWindowList, upstreamWindowList, tissMaskWindowList, i);
                    double rem = reads.getDouble(i);
                    double add = reads.getDouble(i+window);
                    downstreamMean -= rem/window;
                    downstreamMean += add/window;
                    rem = reads.getDouble(i+window+tissMaskWindow);
                    add = reads.getDouble(i+window*2+tissMaskWindow);
                    upstreamMean -= rem/window;
                    upstreamMean += add/window;
                }
                continue;
            }

            if (downstreamMean < minDownstreamReadMean || downstreamMean <= upstreamMean) {
                if (ref.isPlus()) {
                    updateLists(reads, upstreamWindowList, downstreamWindowList, tissMaskWindowList, i);
                    double rem = reads.getDouble(i+window+tissMaskWindow);
                    double add = reads.getDouble(i+window*2+tissMaskWindow);
                    downstreamMean -= rem/window;
                    downstreamMean += add/window;
                    rem = reads.getDouble(i);
                    add = reads.getDouble(i+window);
                    upstreamMean -= rem/window;
                    upstreamMean += add/window;
                } else {
                    updateLists(reads, downstreamWindowList, upstreamWindowList, tissMaskWindowList, i);
                    double rem = reads.getDouble(i);
                    double add = reads.getDouble(i+window);
                    downstreamMean -= rem/window;
                    downstreamMean += add/window;
                    rem = reads.getDouble(i+window+tissMaskWindow);
                    add = reads.getDouble(i+window*2+tissMaskWindow);
                    upstreamMean -= rem/window;
                    upstreamMean += add/window;
                }
                continue;
            }

            KolmogorovSmirnovTest test = new KolmogorovSmirnovTest();
            double[] downstreamAry = EI.wrap(downstreamWindowList.toList()).toDoubleArray();
            double[] upstreamAry = EI.wrap(upstreamWindowList.toList()).toDoubleArray();
            double bothEqualPVal = test.kolmogorovSmirnovTest(upstreamAry, downstreamAry);

            if (bothEqualPVal > upstreamDownstreamEqualityThresh) {
                if (ref.isPlus()) {
                    updateLists(reads, upstreamWindowList, downstreamWindowList, tissMaskWindowList, i);
                    double rem = reads.getDouble(i+window+tissMaskWindow);
                    double add = reads.getDouble(i+window*2+tissMaskWindow);
                    downstreamMean -= rem/window;
                    downstreamMean += add/window;
                    rem = reads.getDouble(i);
                    add = reads.getDouble(i+window);
                    upstreamMean -= rem/window;
                    upstreamMean += add/window;
                } else {
                    updateLists(reads, downstreamWindowList, upstreamWindowList, tissMaskWindowList, i);
                    double rem = reads.getDouble(i);
                    double add = reads.getDouble(i+window);
                    downstreamMean -= rem/window;
                    downstreamMean += add/window;
                    rem = reads.getDouble(i+window+tissMaskWindow);
                    add = reads.getDouble(i+window*2+tissMaskWindow);
                    upstreamMean -= rem/window;
                    upstreamMean += add/window;
                }
                continue;
            }

            double downstreamSd = sampleSd(downstreamAry, downstreamMean);
            double upstreamSd = sampleSd(upstreamAry, upstreamMean);
            double downstreamPVal = test.kolmogorovSmirnovTest(new NormalDistribution(downstreamMean, downstreamSd), EI.wrap(downstreamWindowList.toList()).toDoubleArray());
            double upstreamPVal = test.kolmogorovSmirnovTest(new NormalDistribution(upstreamMean, upstreamSd), EI.wrap(upstreamWindowList.toList()).toDoubleArray());

            if (upstreamPVal < downstreamPVal) {
                Map<String, Double> info = new HashMap<>();
                info.put("downP", downstreamPVal);
                info.put("upP", upstreamPVal);
                foundPeaks.put(i + window + tissMaskWindow / 2, info);
            }

            if (ref.isPlus()) {
                updateLists(reads, upstreamWindowList, downstreamWindowList, tissMaskWindowList, i);
                double rem = reads.getDouble(i+window+tissMaskWindow);
                double add = reads.getDouble(i+window*2+tissMaskWindow);
                downstreamMean -= rem/window;
                downstreamMean += add/window;
                rem = reads.getDouble(i);
                add = reads.getDouble(i+window);
                upstreamMean -= rem/window;
                upstreamMean += add/window;
            } else {
                updateLists(reads, downstreamWindowList, upstreamWindowList, tissMaskWindowList, i);
                double rem = reads.getDouble(i);
                double add = reads.getDouble(i+window);
                downstreamMean -= rem/window;
                downstreamMean += add/window;
                rem = reads.getDouble(i+window+tissMaskWindow);
                add = reads.getDouble(i+window*2+tissMaskWindow);
                upstreamMean -= rem/window;
                upstreamMean += add/window;
            }
        }
    }

    private void updateLists(NumericArray reads, StaticSizeSortedArrayList<Double> upstreamWindowList, StaticSizeSortedArrayList<Double> downstreamWindowList, StaticSizeSortedArrayList<Double> tissMaskWindowList, int i) {
        upstreamWindowList.insertSortedAndDelete(reads.getDouble(i+window), reads.getDouble(i));
        downstreamWindowList.insertSortedAndDelete(reads.getDouble(i+window*2+tissMaskWindow), reads.getDouble(i+window+tissMaskWindow));
        tissMaskWindowList.insertSortedAndDelete(reads.getDouble(i+window+tissMaskWindow), reads.getDouble(i+window));
    }

    private double sampleSd(double[] sample, double mean) {
        double sqrtSum = 0;
        for (int i = 0; i < sample.length; i++) {
            sqrtSum += Math.pow(sample[i] - mean, 2);
        }
        return Math.sqrt(sqrtSum/(sample.length-1));
    }

    private double calculateMean(double[] ary) {
        double sum = ArrayUtils.sum(ary);
        return sum/ary.length;
    }
}
