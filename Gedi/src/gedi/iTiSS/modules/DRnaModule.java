package gedi.iTiSS.modules;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.merger2.TissFile;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.iTiSS.utils.sortedNodesList.StaticSizeSortedArrayList;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.dynamic.impl.DoubleDynamicObject;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import gedi.util.r.RRunner;

import static gedi.iTiSS.utils.TiSSUtils.log2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class DRnaModule extends ModuleBase {
    private int windowSize;
    private double pseudoCount;
    private double peakCallThreshold;
    private boolean machineLearning;
    private int cleanupThresh;
    private int minReadNum;

    public DRnaModule(int windowSize, double pseudoCount, double peakCallThreshold, int cleanupThresh, int minReadNum, Data lane, boolean useMM, String name) {
        super(name, lane);
        this.windowSize = windowSize;
        this.pseudoCount = pseudoCount;
        this.peakCallThreshold = peakCallThreshold;
        this.machineLearning = useMM;
        this.cleanupThresh = cleanupThresh;
        this.minReadNum = minReadNum;
    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        NumericArray xrn1 = data[0];
        List<PeakAndPos> mmData = machineLearning ? new LinkedList<>() : null;

        System.err.println("Adding pseudocounts");
        addPseudoCount(xrn1, pseudoCount);

        System.err.println("Calculating sum");
        List<PeakAndPos> ptp = new LinkedList<>();
        double[] tmp = xrn1.toDoubleArray(0, windowSize*2+1);
        for (int i = 0; i < tmp.length; i++) {
            ptp.add(new PeakAndPos(i, tmp[i]));
        }
        StaticSizeSortedArrayList<PeakAndPos> window = new StaticSizeSortedArrayList<>(ptp, PeakAndPos::compareTo);
        double windowSum = 0;
        for (int i = 0; i < window.getSize(); i++) {
            windowSum += window.getValueAtIndex(i).getValue();
        }

        System.err.println("Starting TisS identification");
        List<MutablePair<Integer, Map<String, Double>>> tiss = new ArrayList<>();
        for (int i = windowSize; i < xrn1.length() - (windowSize+1); i++) {
            if (mmData != null && i%10000000 == 0) {
                System.err.print("Progress " + i + "/" + xrn1.length() + ", mmData size: " + mmData.size() + "\r");
            }
            double thresholdPeak = (windowSum-xrn1.getDouble(i)) / (window.getSize()-1);
            if (mmData != null) {
                if (xrn1.getDouble(i) > pseudoCount && xrn1.getDouble(i) > thresholdPeak) {
                    mmData.add(new PeakAndPos(i, log2(xrn1.getDouble(i) / thresholdPeak)));
                }
            } else {
                if (thresholdPeak >= pseudoCount && xrn1.getDouble(i) > peakCallThreshold * thresholdPeak && xrn1.getDouble(i) >= minReadNum) {
                    Map<String, Double> infos = new HashMap<>();
                    infos.put("peak height", xrn1.getDouble(i));
                    infos.put("threshold peak height", thresholdPeak);
                    tiss.add(new MutablePair<>(i, infos));
                }
            }

            window.insertSortedAndDelete(new PeakAndPos(i+windowSize+1, xrn1.getDouble(i+windowSize+1)), new PeakAndPos(i-windowSize, xrn1.getDouble(i-windowSize)));
            windowSum += xrn1.getDouble(i+windowSize+1);
            windowSum -= xrn1.getDouble(i-windowSize);
        }

        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());
        List<Double> mmDatOut = this.allMmData.computeIfAbsent(ref, k-> new ArrayList<>());

        if (mmData != null) {
            mmDatOut.addAll(mmData.stream().map(PeakAndPos::getValue).collect(Collectors.toList()));
            mmData.forEach(m -> {
                Map<String, Double> info = new HashMap<>();
                info.put(TissFile.FOLD_CHANGE_COLUMN_NAME, m.getValue());
                info.put(TissFile.READ_COUNT_COLUMN_NAME, xrn1.getDouble(m.getPos()));
                foundPeaksNew.put(m.getPos(), info);
            });
        } else {
            for (MutablePair<Integer, Map<String, Double>> i : tiss) {
                foundPeaksNew.put(i.Item1, i.Item2);
            }
        }

        System.err.println(foundPeaksNew.keySet().size() + " peaks found for ref " + ref.toPlusMinusString());
    }

    private List<PeakAndPos> filterMultiPeaks(List<PeakAndPos> lst) {
        Map<Long, Integer> peakCount = new HashMap<>();

        for (PeakAndPos pap : lst) {
            long value = Double.doubleToLongBits(pap.getValue());
            if (peakCount.containsKey(value)) {
                peakCount.put(value, peakCount.get(value) + 1);
            } else{
                peakCount.put(value, 1);
            }
        }

        return EI.wrap(lst).filter(pap -> peakCount.get(Double.doubleToLongBits(pap.getValue())) < cleanupThresh).list();
    }

    @Override
    public void calculateMLResults(String prefix, boolean plot) throws IOException {
        LineWriter writerup = new LineOrientedFile(prefix + "sparsePeakThresholdData.tsv").write();
        writerup.writeLine("Ref\tPos\t" + TissFile.FOLD_CHANGE_COLUMN_NAME + "\t" + TissFile.READ_COUNT_COLUMN_NAME);

        List<PeakAndPos> mlData = new ArrayList<>();
        for (ReferenceSequence ref : res.keySet()) {
            Map<Integer, Map<String, Double>> tss = res.get(ref);
            for (int tssPos : tss.keySet()) {
                Map<String, Double> info = tss.get(tssPos);
                mlData.add(new PeakAndPos(tssPos, info.get(TissFile.FOLD_CHANGE_COLUMN_NAME)));
                writerup.writeLine(ref.toPlusMinusString() + "\t" + tssPos + "\t" + info.get(TissFile.FOLD_CHANGE_COLUMN_NAME) + "\t" + info.get(TissFile.READ_COUNT_COLUMN_NAME));
            }
        }
        writerup.close();
        int movingAverage = (int)(mlData.size()*0.2);
        double upThresh = TiSSUtils.calculateThreshold(mlData, movingAverage);

        for (ReferenceSequence ref : res.keySet()) {
            Map<Integer, Map<String,Double>> filtered = new HashMap<>();

            Map<Integer, Map<String, Double>> tss2val = res.get(ref);
            if (cleanupThresh > 0) {
                // Triple<Pos, Fold-change, Read-count>
                List<MutableTriple<Integer, Double, Double>> vals = EI.wrap(res.get(ref).keySet()).map(t -> new MutableTriple<>(t, res.get(ref).get(t).get(TissFile.FOLD_CHANGE_COLUMN_NAME), res.get(ref).get(t).get(TissFile.READ_COUNT_COLUMN_NAME))).list();
                vals = TiSSUtils.cleanUpMultiValueDataTriple(vals, cleanupThresh);
                Map<Integer,Map<String,Double>> tss2valFiltered = new HashMap<>();
                EI.wrap(vals).forEachRemaining(u -> {
                    Map<String, Double> tmpmap = tss2valFiltered.computeIfAbsent(u.Item1, absent->new HashMap<>());
                    tmpmap.put(TissFile.FOLD_CHANGE_COLUMN_NAME, u.Item2);
                    tmpmap.put(TissFile.READ_COUNT_COLUMN_NAME, u.Item3);
                });
                tss2val = tss2valFiltered;
            }

            for (int tss : tss2val.keySet()) {
                if (tss2val.get(tss).get(TissFile.FOLD_CHANGE_COLUMN_NAME) > upThresh && tss2val.get(tss).get(TissFile.READ_COUNT_COLUMN_NAME) >= minReadNum) {
                    filtered.put(tss, tss2val.get(tss));
                }
            }
            res.put(ref, filtered);
        }

        if (plot) {
            plotData(prefix, prefix + "sparsePeakThresholdData.tsv", upThresh);
        }
    }

    private void plotData(String prefix, String dataFilePath, double thresh) throws IOException {
        RRunner r = new RRunner(prefix+".plotSparsePeakThresholdData.R");
        r.set("dataFile",dataFilePath);
        r.set("pdfFile",prefix+".sparsePeakThresholds.pdf");
        r.set("prefix", prefix);
        r.set("customThresh", new DoubleDynamicObject(thresh));
        r.set("threshAutoParam", new DoubleDynamicObject(thresh));
        r.set("manualSelectionFile", prefix + "sparsePeakManualSelection.tsv");
        r.set("manualSelectionPdf", prefix + "sparsePeakManuelThresh.pdf");
        r.addSource(getClass().getResourceAsStream("/resources/plotSparseThresh.R"));
        r.run(false);
    }

    private void addPseudoCount(NumericArray ary, double pseudoCount) {
        for (int i = 0; i < ary.length(); i++) {
            ary.setDouble(i, ary.getDouble(i)+pseudoCount);
        }
    }
}
