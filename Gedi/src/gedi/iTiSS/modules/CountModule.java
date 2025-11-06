package gedi.iTiSS.modules;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.functions.NumericArrayFunction;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class CountModule extends ModuleBase {
    private int upstreamCount;
    private int downstreamCount;
    private double countThreshold;

    public CountModule(int upstreamCount, int downstreamCount, double countThreshold, Data lane, String moduleName) {
        super(moduleName, lane);
        this.upstreamCount = upstreamCount;
        this.downstreamCount = downstreamCount;
        this.countThreshold = countThreshold;
    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        NumericArray dat = data[0];
        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());

        int pos = ref.isPlus() ? upstreamCount : downstreamCount;
        double currentSum = dat.slice(0, upstreamCount+downstreamCount).sum();
        for (int i = 0; i < dat.length() - (upstreamCount+downstreamCount+1); i++) {
            currentSum += dat.getDouble(i + downstreamCount+upstreamCount);
            if (currentSum > countThreshold) {
                Map<String, Double> info = new HashMap<>();
                info.put("Number of peaks", currentSum);
                foundPeaksNew.put(pos, info);
            }
            currentSum -= dat.getDouble(i);
            pos++;
        }
    }
}
