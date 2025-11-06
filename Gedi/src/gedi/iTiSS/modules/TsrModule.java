package gedi.iTiSS.modules;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.iTiSS.data.Data;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.decorators.NumericArraySlice;
import gedi.util.datastructure.array.functions.NumericArrayFunction;
import gedi.util.functions.EI;

import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

public class TsrModule extends ModuleBase {
    private int tsrSize;
    private int minReads;
    private int bufferUpstream;
    private int bufferDownstream;

    public TsrModule(int tsrSize, int minReads, int bufferUpstream, int bufferDownstream, Data lane, String moduleName) {
        super(moduleName, lane);
        this.tsrSize = tsrSize;
        this.minReads = minReads;
        this.bufferUpstream = bufferUpstream;
        this.bufferDownstream = bufferDownstream;
    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        NumericArray readCounts = data[0];
        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());

        final GenomicRegion[] gr = {new ArrayGenomicRegion()};
        EI.seq(0,readCounts.length()).
                sort((o1,o2) -> Double.compare(readCounts.getDouble(o2), readCounts.getDouble(o1))).
                forEachRemaining(i -> {
                    if (i-tsrSize < 0) {
                        return;
                    }
                    if (gr[0].contains(i-tsrSize) || gr[0].contains(i+1)) {
                        return;
                    }
                    NumericArraySlice slice = readCounts.slice(i - tsrSize + 1, i + 1);
                    if (slice.sum() >= minReads) {
                        gr[0] = gr[0].union(new ArrayGenomicRegion(i-tsrSize, i+1));
                        int maxPos = argMax(slice) + (i-tsrSize+1);
                        Map<String, Double> info = new HashMap<>();
                        info.put("tsr max", readCounts.getDouble(maxPos));
                        info.put("tsr sum", slice.sum());
                        foundPeaksNew.put(maxPos, info);
                    }
                });
    }

    private int argMax(NumericArray ary) {
        double max = Double.MIN_VALUE;
        int maxPos = -1;
        for (int i = 0; i < ary.length(); i++) {
            if (ary.getDouble(i) > max) {
                max = ary.getDouble(i);
                maxPos = i;
            }
        }
        return maxPos;
    }
}
