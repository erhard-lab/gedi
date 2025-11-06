package gedi.iTiSS.modules;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.data.DataWrapper;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutablePair;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public abstract class ModuleBase {
    protected Data lane;
    protected Map<ReferenceSequence, Map<Integer, Map<String,Double>>> res = new HashMap<>();
    protected Map<ReferenceSequence, List<Double>> allMmData = new HashMap<>();
    protected Map<ReferenceSequence, Double> thresholds = new HashMap<>();
    protected double globalThreshold;
    protected String moduleName;

    public ModuleBase(String moduleName, Data lane){
        this.lane = lane;
        this.moduleName = moduleName;
    }

    public abstract void findTiSS(NumericArray[] data, ReferenceSequence ref);

    public final int laneHashCode() {
        Set<Integer> tmp = new HashSet<>();
        for (int i : lane.getLane()) {
            tmp.add(i);
        }
        return tmp.hashCode();
    }

    public Map<ReferenceSequence, Map<Integer, Map<String, Double>>> getResultsNew() {
        return res;
    }

    public Map<ReferenceSequence, List<Double>> getAllMmData() {
        return allMmData;
    }

    public Map<ReferenceSequence, Double> getThresholds() {
        return thresholds;
    }

    public Data getLane() {
        return lane;
    }

    public String getModuleName() {
        return moduleName;
    }

    public double getGlobalThreshold() {
        return globalThreshold;
    }

    public void calculateMLResults(String prefix, boolean plot) throws IOException {
        return;
    }
}
