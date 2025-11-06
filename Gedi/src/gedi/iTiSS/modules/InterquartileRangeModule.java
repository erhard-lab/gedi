package gedi.iTiSS.modules;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.iTiSS.utils.sortedNodesList.StaticSizeSortedArrayList;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

import static gedi.iTiSS.utils.TiSSUtils.log2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class InterquartileRangeModule extends ModuleBase{

    private double iqrThreshold;
    private int windowSize;
    private double minReadDensity;
    private final double ZERO_THRESHOLD = 0.5;
    private boolean machineLearning;
    private boolean dumpMmVal;

    private int firstQuartileEnd;
    private int thirdQuartileStart;

    public InterquartileRangeModule(double iqrThreshold, int windowSize, double minReadDensity, Data lane,
                                    boolean useMM, boolean dumpMmVal, String name) {
        super(name, lane);
        this.iqrThreshold = iqrThreshold;
        this.windowSize = windowSize;
        this.minReadDensity = minReadDensity;
        this.machineLearning = useMM;
        this.dumpMmVal = dumpMmVal;
        firstQuartileEnd = windowSize/4;
        thirdQuartileStart = (int)((((double)windowSize)/4.)*3.);
    }

    @Override
    public void findTiSS(NumericArray[] inData, ReferenceSequence ref) {
        // Should only be one lane, the total lane
        NumericArray data = inData[0];
        List<PeakAndPos> mmDataL = machineLearning ? new ArrayList<>() : null;
        List<PeakAndPos> mmDataC = machineLearning ? new ArrayList<>() : null;
        List<PeakAndPos> mmDataR = machineLearning ? new ArrayList<>() : null;
        StaticSizeSortedArrayList<Double> windowLeft = new StaticSizeSortedArrayList<>(EI.wrap(data.toDoubleArray(0, windowSize)).toList(), Double::compareTo);
        int leftCenter = (windowSize % 2 == 0) ? windowSize/2 : windowSize/2+1;
        int rightCenter = (windowSize*2-leftCenter)+1;
        List<Double> lst = EI.wrap(data.toDoubleArray(leftCenter, windowSize)).chain(EI.wrap(data.toDoubleArray(windowSize+1, rightCenter))).toList();
        StaticSizeSortedArrayList<Double> windowCenter = new StaticSizeSortedArrayList<>(lst, Double::compareTo);
        StaticSizeSortedArrayList<Double> windowRight = new StaticSizeSortedArrayList<>(EI.wrap(data.toDoubleArray(windowSize+1, windowSize*2+1)).toList(), Double::compareTo);
        int leftWindowZeroCount = zeroCount(data.toDoubleArray(0, windowSize), ZERO_THRESHOLD);
        int rightWindowZeroCount = zeroCount(data.toDoubleArray(windowSize+1, windowSize*2+1), ZERO_THRESHOLD);
        int centerWindowZeroCount = zeroCount(ArrayUtils.concat(data.toDoubleArray(leftCenter, windowSize), data.toDoubleArray(windowSize+1, rightCenter)), ZERO_THRESHOLD);

        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());
        List<Double> mmDatOut = this.allMmData.computeIfAbsent(ref, k-> new ArrayList<>());

        for (int i = windowSize; i < data.length() - (windowSize+1); i++) {
            double valueOfInterest = data.getDouble(i);
            boolean ignoreLeft = notEnoughReads(leftWindowZeroCount);
            boolean ignoreCenter = notEnoughReads(centerWindowZeroCount);
            boolean ignoreRight = notEnoughReads(rightWindowZeroCount);
            if (mmDataL != null) {
                if (!ignoreCenter && !ignoreRight && !ignoreLeft) {
                    if (valueGreaterIQRthreshold(valueOfInterest, windowLeft, 1.0)) {
                        mmDataL.add(new PeakAndPos(i, log2(valueIqr(valueOfInterest, windowLeft) / windowIqr(windowLeft))));
                    }
                    if (valueGreaterIQRthreshold(valueOfInterest, windowCenter, 1.0)) {
                        mmDataC.add(new PeakAndPos(i, log2(valueIqr(valueOfInterest, windowCenter) / windowIqr(windowCenter))));
                    }
                    if (valueGreaterIQRthreshold(valueOfInterest, windowRight, 1.0)) {
                        mmDataR.add(new PeakAndPos(i, log2(valueIqr(valueOfInterest, windowRight) / windowIqr(windowRight))));
                    }
                }
            } else {
                if ((valueGreaterIQRthreshold(valueOfInterest, windowLeft, iqrThreshold) || ignoreLeft) &&
                        (valueGreaterIQRthreshold(valueOfInterest, windowCenter, iqrThreshold) || ignoreCenter) &&
                        (valueGreaterIQRthreshold(valueOfInterest, windowRight, iqrThreshold) || ignoreRight) &&
                        !(ignoreLeft && ignoreCenter && ignoreRight)) {
                    Map<String, Double> infos = new HashMap<>();
                    infos.put("peak height", valueOfInterest);
                    infos.put("read coverage right", (1 - (rightWindowZeroCount / (double) windowSize)));
                    infos.put("read coverage center", (1 - (centerWindowZeroCount / (double) windowSize)));
                    infos.put("read coverage left", (1 - (leftWindowZeroCount / (double) windowSize)));
                    infos.put("IQR right", (windowIqr(windowRight)));
                    infos.put("IQR threshold right", (windowIqr(windowRight) * iqrThreshold));
                    infos.put("IQR value right", (valueIqr(valueOfInterest, windowRight)));
                    infos.put("IQR center", (windowIqr(windowCenter)));
                    infos.put("IQR threshold center", (windowIqr(windowCenter) * iqrThreshold));
                    infos.put("IQR value center", (valueIqr(valueOfInterest, windowCenter)));
                    infos.put("IQR left", (windowIqr(windowLeft)));
                    infos.put("IQR threshold left", (windowIqr(windowLeft) * iqrThreshold));
                    infos.put("IQR value left", (valueIqr(valueOfInterest, windowLeft)));
                    foundPeaksNew.put(i, infos);
                }
            }
            if (data.getDouble(i) <= ZERO_THRESHOLD) {
                leftWindowZeroCount++;
            }
            if (data.getDouble(i - windowSize) <= ZERO_THRESHOLD) {
                leftWindowZeroCount--;
            }
            if (data.getDouble(rightCenter + (i - windowSize)) <= ZERO_THRESHOLD) {
                centerWindowZeroCount++;
            }
            if (data.getDouble(leftCenter + (i - windowSize)) <= ZERO_THRESHOLD) {
                centerWindowZeroCount--;
            }
            if (data.getDouble(i + windowSize + 1) <= ZERO_THRESHOLD) {
                rightWindowZeroCount++;
            }
            if (data.getDouble(i + 1) <= ZERO_THRESHOLD) {
                rightWindowZeroCount--;
            }

            windowLeft.insertSortedAndDelete(data.getDouble(i), data.getDouble(i-windowSize));
            windowCenter.insertSortedAndDelete(data.getDouble(rightCenter+(i-windowSize)), data.getDouble(leftCenter+(i-windowSize)));
            windowRight.insertSortedAndDelete(data.getDouble(i+windowSize+1), data.getDouble(i+1));

        }

        if (mmDataL != null) {
//            PeakClusterer peakClusterer = new PeakClusterer(mmDataC);
//            peakClusterer.cluster();
//            double thresh = peakClusterer.getThreshold().getMeanThresh();
//            int movingAverage = (int)(mmDataC.size()*0.2);
//            double thresh = dumpMmVal ? -999.0 : TiSSUtils.calculateThreshold(mmDataC, movingAverage);
            mmDatOut.addAll(mmDataC.stream().map(PeakAndPos::getValue).collect(Collectors.toList()));
            mmDataC.forEach(m -> {
                Map<String, Double> info = new HashMap<>();
                info.put("Fold-change", m.getValue());
                foundPeaksNew.put(m.getPos(), info);
            });
//            this.thresholds.put(ref, thresh);
//            System.err.println("--== Threshold: " + thresh + " ==--");
//            System.err.println("--== Moving avg: " + movingAverage + " ==--");
//            System.err.println("Threshold! Low: " + peakClusterer.getThreshold().getLowThresh() + ", High: " + peakClusterer.getThreshold().getHighThresh());

//            List<PeakAndPos> peaksL = mmDataL.stream().filter(a -> a.getValue() > thresh).collect(Collectors.toList());
//            List<PeakAndPos> nonPeaksL = mmDataL.stream().filter(a -> a.getValue() < thresh).collect(Collectors.toList());
//            List<PeakAndPos> peaksR = mmDataR.stream().filter(a -> a.getValue() > thresh).collect(Collectors.toList());
//            List<PeakAndPos> nonPeaksR = mmDataR.stream().filter(a -> a.getValue() < thresh).collect(Collectors.toList());
//            List<PeakAndPos> peaksC = mmDataC.stream().filter(a -> a.getValue() > thresh).collect(Collectors.toList());
//            List<PeakAndPos> nonPeaksC = mmDataC.stream().filter(a -> a.getValue() < thresh).collect(Collectors.toList());
//            List<PeakAndPos> peaksC = peakClusterer.getTopCluster();
//            List<PeakAndPos> nonPeaksC = peakClusterer.getBottomCluster();

//            foundPeaksNew.clear();
//            PeakClusterer peakClustererL = new PeakClusterer(mmDataL);
//            peakClustererL.cluster();
//            List<PeakAndPos> peaksL = peakClustererL.getTopCluster();
//            List<PeakAndPos> nonPeaksL = peakClustererL.getBottomCluster();
//            PeakClusterer peakClustererC = new PeakClusterer(mmDataC);
//            peakClustererC.cluster();
//            List<PeakAndPos> peaksC = peakClustererC.getTopCluster();
//            List<PeakAndPos> nonPeaksC = peakClustererC.getBottomCluster();
//            PeakClusterer peakClustererR = new PeakClusterer(mmDataR);
//            peakClustererR.cluster();
//            List<PeakAndPos> peaksR = peakClustererR.getTopCluster();
//            List<PeakAndPos> nonPeaksR = peakClustererR.getBottomCluster();
//            identifyFoundPeaksMM(peaksL, peaksC, peaksR, nonPeaksC, nonPeaksR, foundPeaksNew);
//            identifyFoundPeaksMM(peaksC, peaksL, peaksR, nonPeaksL, nonPeaksR, foundPeaksNew);
//            identifyFoundPeaksMM(peaksR, peaksL, peaksC, nonPeaksL, nonPeaksC, foundPeaksNew);
        }

//        System.err.println("Finished searching for peaks in the IQR module.");
//        System.err.println("Found peaks: " + foundPeaks.size());
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

    private void identifyFoundPeaksMM(List<PeakAndPos> query, List<PeakAndPos> peakA, List<PeakAndPos> peakB,
                                      List<PeakAndPos> nonPeakA, List<PeakAndPos> nonPeakB,
                                      Map<Integer, Map<String, Double>> foundPeaksNew) {
        for (PeakAndPos pap : query) {
            if (!peakA.contains(pap) || !peakB.contains(pap)) {
                continue;
            }
            removePAP(pap, peakA, peakB);
            Map<String, Double> info = new HashMap<>();
            info.put("LogIQRfc", pap.getValue());
            foundPeaksNew.put(pap.getPos(), info);
        }
    }

    private void removePAP(PeakAndPos pap, List<PeakAndPos> lst1, List<PeakAndPos> lst2) {
        lst1.remove(pap);
        lst2.remove(pap);
    }

    private boolean notEnoughReads(int zeroCount) {
        return ((double)zeroCount)/((double)windowSize) >= 1.-minReadDensity;
    }

    private boolean valueGreaterIQRthreshold(double value, StaticSizeSortedArrayList<Double> window, double thresh) {
        double iqr = windowIqr(window);
        return valueIqr(value, window) > thresh*iqr;
    }

    private double windowIqr(StaticSizeSortedArrayList<Double> window) {
        return window.getValueAtIndex(thirdQuartileStart) - window.getValueAtIndex(firstQuartileEnd);
    }

    private double valueIqr(double value, StaticSizeSortedArrayList<Double> window) {
        return value-window.getValueAtIndex(thirdQuartileStart);
    }

    private int zeroCount(double[] array, double zeroThreshold) {
        int zeros = 0;
        for (int i = 0; i<array.length; i++) {
            if (array[i] <= zeroThreshold) {
                zeros++;
            }
        }
        return zeros;
    }

    private class ReadCount {
        double readCount;
        int peakCount;

        public ReadCount(double readCount) {
            this.readCount = readCount;
            this.peakCount = 0;
        }
    }
}
