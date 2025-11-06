package gedi.iTiSS.modules;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.utils.ArrayUtils2;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.mutable.MutablePair;

import java.util.ArrayList;
import java.util.List;

/**
 * This module detects peaks in histone datasets. This means that not single nucleotide positions are called
 * but regions in density plots, where significantly more reads map to.
 */
public class HistoneDetectionModule extends ModuleBase {
    private int smoothingRange = 500;
    private double minReads = 10;
    private int minDist = 1000;

    public HistoneDetectionModule(String moduleName, Data lane) {
        super(moduleName, lane);
    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        double[] smoothed = ArrayUtils2.smooth(data[0].toDoubleArray(), smoothingRange);

        boolean isIncreasing = false;
        double lastVal = smoothed[smoothingRange/2];
        List<MutablePair<Integer, Double>> potentialPeaks = new ArrayList<>();
        for (int i = smoothingRange/2; i < smoothed.length; i++) {
            if (smoothed[i] > lastVal) {
                isIncreasing = true;
            } else if (isIncreasing) {
                potentialPeaks.add(new MutablePair<>(i, lastVal));
                isIncreasing = false;
            }
            lastVal = smoothed[i];
        }

        List<MutablePair<Integer, Double>> peaks = new ArrayList<>();
        MutablePair<Integer, Double> lastPeak = potentialPeaks.get(0);
        for (int i = 1; i < potentialPeaks.size(); i++) {
            if (lastPeak.Item2 < minReads) {
                lastPeak = potentialPeaks.get(i);
                continue;
            }
            if (potentialPeaks.get(i).Item1 - lastPeak.Item1 < minDist) {

            }
        }
    }


}
