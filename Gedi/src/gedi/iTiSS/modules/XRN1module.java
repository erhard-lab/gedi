package gedi.iTiSS.modules;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.iTiSS.utils.machineLearning.PeakClusterer;
import gedi.iTiSS.utils.sortedNodesList.StaticSizeSortedArrayList;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutablePair;

import static gedi.iTiSS.utils.TiSSUtils.log2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class XRN1module extends ModuleBase {
    private int windowSize;
    private double pseudoCount;
    private double maxWater2XRN1FoldChange;
    private int thresholdPeakNumber;
    private double peakCallThreshold;
    private boolean waterCheck;
    private int peakProximity;
    private int minReadNum;
    private boolean machineLearning;
    private boolean dumpMmVal;

    public XRN1module(int windowSize, double pseudoCount, double maxWater2XRN1FoldChange, int thresholdPeakNumber,
                      double peakCallThreshold, Data lane, boolean waterCheck, int peakProximity, int minReadNum,
                      boolean useMM, boolean dumpMmVal, String name) {
        super(name, lane);
        this.windowSize = windowSize;
        this.pseudoCount = pseudoCount;
        this.maxWater2XRN1FoldChange = maxWater2XRN1FoldChange;
        this.thresholdPeakNumber = thresholdPeakNumber;
        this.peakCallThreshold = peakCallThreshold;
        this.waterCheck = waterCheck;
        this.peakProximity = peakProximity;
        this.minReadNum = minReadNum;
        this.machineLearning = useMM;
        this.dumpMmVal = dumpMmVal;
    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        NumericArray xrn1 = data[0];
        NumericArray water = waterCheck ? data[1] : null;
        List<PeakAndPos> mmData = machineLearning ? new ArrayList<>() : null;

        addPseudoCount(xrn1, pseudoCount);

        List<PeakAndPos> ptp = new LinkedList<>();
        double[] tmp = xrn1.toDoubleArray(0, windowSize*2+1);
        for (int i = 0; i < tmp.length; i++) {
            ptp.add(new PeakAndPos(i, tmp[i]));
        }
        StaticSizeSortedArrayList<PeakAndPos> window = new StaticSizeSortedArrayList<>(ptp, PeakAndPos::compareTo);

        List<MutablePair<Integer, Map<String, Double>>> tiss = new ArrayList<>();
        for (int i = windowSize; i < xrn1.length() - (windowSize+1); i++) {
            PeakAndPos thresholdPeak = thresholdPeak(window);
            if (mmData != null) {
                if (xrn1.getDouble(i) > pseudoCount && xrn1.getDouble(i)/thresholdPeak.getValue() > 1) {
                    mmData.add(new PeakAndPos(i, log2(xrn1.getDouble(i)/thresholdPeak.getValue())));
                }
            }
//            if (i == 195974 && !ref.isSense()) {
//                for (int j = i-20; j < i+20; j++) {
//                    System.err.println("Pos: " + j + ", Height: " + xrn1.getDouble(j));
//                }
//                System.err.println("=====");
//                System.err.println("Threshold Peak Pos: " + thresholdPeak.getTssPos());
//                System.err.println("Threshold Peak Height: " + thresholdPeak.getValue());
//                System.err.println("Greater Pseudocount: " + (thresholdPeak.getValue() >= pseudoCount ? "Yes" : "No"));
//                System.err.println("Greater Threshold: " + (xrn1.getDouble(i) > peakCallThreshold*thresholdPeak.getValue() ? "Yes" : "No"));
//                System.err.println("=====");
//            }
            if (thresholdPeak.getValue() >= pseudoCount && xrn1.getDouble(i) > peakCallThreshold*thresholdPeak.getValue()) {
                Map<String, Double> infos = new HashMap<>();
                infos.put("peak height", xrn1.getDouble(i));
                infos.put("threshold peak pos", (double)thresholdPeak.getPos());
                infos.put("threshold peak height", thresholdPeak.getValue());
                tiss.add(new MutablePair<>(i, infos));
            }
            window.insertSortedAndDelete(new PeakAndPos(i+windowSize+1, xrn1.getDouble(i+windowSize+1)), new PeakAndPos(i-windowSize, xrn1.getDouble(i-windowSize)));
        }

        System.err.println(tiss.size() + " peaks found in XRN1 with ref " + ref.toPlusMinusString());

        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());
        List<Double> mmDatOut = this.allMmData.computeIfAbsent(ref, k-> new ArrayList<>());

        if (mmData != null) {
            mmDatOut.addAll(mmData.stream().map(PeakAndPos::getValue).collect(Collectors.toList()));
            int movingAverage = (int)(mmData.size()*0.2);
            final double thresh = dumpMmVal ? -9999 : TiSSUtils.calculateThreshold(mmData, movingAverage);
            this.thresholds.put(ref, thresh);
            List<PeakAndPos> peaks = mmData.stream().filter(a -> a.getValue() > thresh).collect(Collectors.toList());
            System.err.println("--== Threshold: " + thresh + " ==--");
            System.err.println("--== Moving avg: " + movingAverage + " ==--");
            tiss = EI.wrap(peaks)
                    .map(a -> {
                        Map<String, Double> tmpMap = new HashMap<>();
                        tmpMap.put("Fold-change", a.getValue());
                        return new MutablePair<>(a.getPos(), tmpMap);
                    }).list();
        }

        for (MutablePair<Integer, Map<String, Double>> i : tiss) {
            if (waterCheck) {
//                if ((xrn1.getDouble(i.Item1) + xrn1Mean) * maxWater2XRN1FoldChange >= water.getDouble(i.Item1)) {
                if ((xrn1.getDouble(i.Item1)) * maxWater2XRN1FoldChange >= water.getDouble(i.Item1)) {
                    i.Item2.put("water height", water.getDouble(i.Item1));
                    foundPeaksNew.put(i.Item1, i.Item2);
                }
            } else {
                i.Item2.put("water height", Double.NaN);
                foundPeaksNew.put(i.Item1, i.Item2);
            }
        }

        System.err.println(foundPeaksNew.keySet().size() + " peaks left after reduction with water for ref " + ref.toPlusMinusString());
    }

    @Override
    public void calculateMLResults(String prefix, boolean plot) throws IOException {
        LineWriter writer = new LineOrientedFile(prefix + ".tsv").write();
        writer.writeLine("Ref\tValue");
        for (ReferenceSequence ref : getAllMmData().keySet()) {
            for (double d : getAllMmData().get(ref)) {
                writer.writeLine(ref.toPlusMinusString() + "\t" + d);
            }
        }
        writer.close();
        List<PeakAndPos> mlData = getAllMmData().values().stream()
                .map(lst -> new ArrayList<>(lst.stream().map(d -> new PeakAndPos(-1, d)).collect(Collectors.toList())))
                .reduce(new ArrayList<>(), (out, next) -> {
                    out.addAll(next);
                    return out;
                });
        int movingAverage = (int)(mlData.size()*0.2);
        System.err.println("Using " + mlData.size() + " values to calculate threshold.");
        final double thresh = TiSSUtils.calculateThreshold(mlData, movingAverage);
        globalThreshold = thresh;
        thresholds.put(Chromosome.obtain("global"), thresh);

        for (ReferenceSequence ref : res.keySet()) {
            Map<Integer, Map<String,Double>> filtered = new HashMap<>();
            for (int tss : res.get(ref).keySet()) {
                if (res.get(ref).get(tss).get("Fold-change") > thresh) {
                    filtered.put(tss, res.get(ref).get(tss));
                }
            }
            res.put(ref, filtered);
        }
    }

    private PeakAndPos thresholdPeak(StaticSizeSortedArrayList<PeakAndPos> ptpList) {
        PeakAndPos first = ptpList.getValueAtIndex(ptpList.getSize()-1);
        PeakAndPos threshold = null;
        int pos = ptpList.getSize()-2;
        int thresholdPos = 0;
        while (threshold == null) {
            if (Math.abs(ptpList.getValueAtIndex(pos).getPos()-first.getPos())>peakProximity) {
                thresholdPos++;
            }
            if (thresholdPos == thresholdPeakNumber) {
                threshold = ptpList.getValueAtIndex(pos);
            }
            pos--;
        }
        return threshold;
    }

    private void addPseudoCount(NumericArray ary, double pseudoCount) {
        for (int i = 0; i < ary.length(); i++) {
            ary.setDouble(i, ary.getDouble(i)+pseudoCount);
        }
    }
}
