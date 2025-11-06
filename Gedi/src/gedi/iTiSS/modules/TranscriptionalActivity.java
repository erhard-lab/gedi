package gedi.iTiSS.modules;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.merger2.TissFile;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.iTiSS.utils.sortedNodesList.StaticSizeSortedArrayList;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.dynamic.impl.DoubleDynamicObject;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.testing.FishersExact;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import gedi.util.r.RRunner;

import java.io.IOException;
import java.util.*;

public class TranscriptionalActivity extends ModuleBase {
    private double significanceThresh;
    private int windowSize;
    private double minReadDensity;
    private final double ZERO_THRESHOLD = 0.5;
    private boolean useML;
    private int cleanupThresh;
    private String writerPath;
    private LineWriter mlWriter;
    private int minReadNum;

    public TranscriptionalActivity(double significanceThresh, int windowSize, double minReadDensity, int cleanupThresh, int minReadNum, boolean useML, String prefix, Data lane, String name) {
        super(name, lane);
        this.significanceThresh = significanceThresh;
        this.windowSize = windowSize;
        this.useML = useML;
        this.minReadDensity = minReadDensity;
        this.cleanupThresh = cleanupThresh;
        this.minReadNum = minReadNum;
        if (useML) {
            writerPath = prefix + "densityThresholdData.tsv";
            mlWriter = new LineOrientedFile(writerPath).write();
            try {
                mlWriter.writeLine("Ref\tPos\tValue\tReadCount");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public void findTiSS(NumericArray[] inData, ReferenceSequence ref) {
        //should only be one lane, the total lane
        NumericArray data = inData[0];
        StaticSizeSortedArrayList<Double> windowUpstream = new StaticSizeSortedArrayList<>(EI.wrap(data.toDoubleArray(0, windowSize)).toList(), Double::compareTo);
        StaticSizeSortedArrayList<Double> windowDownstream = new StaticSizeSortedArrayList<>(EI.wrap(data.toDoubleArray(windowSize+1, windowSize*2+1)).toList(), Double::compareTo);
        if (ref.isMinus()) {
            StaticSizeSortedArrayList<Double> tmp = windowUpstream;
            windowUpstream = windowDownstream;
            windowDownstream = tmp;
        }
        double downstreamWindowSum = sum(data, windowSize+1, windowSize*2+1);
        if (ref.isMinus()) {
            downstreamWindowSum = sum(data, 0, windowSize);
        }
        double threshold = downstreamWindowSum/(double)windowSize;
        int zeroReadsInRightWindow = getZeroReads(data, windowSize+1, windowSize*2+1);

        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());
        List<Double> mmDatOut = this.allMmData.computeIfAbsent(ref, k-> new ArrayList<>());
        List<MutableTriple<Integer, Double, Double>> posPVal = new ArrayList<>(data.length()/100);

        for (int i = windowSize; i < data.length() - (windowSize+1); i++) {
            if (useML && i%10000000 == 0) {
                System.err.print("Progress " + i + "/" + data.length() + ", mmData size: " + posPVal.size() + "\r");
            }

            if (zeroReadsInRightWindow < windowDownstream.getSize()) {
                if (zeroReadsInRightWindow / (double) windowSize < 1. - minReadDensity) {
                    threshold = threshold < 1. ? 1 : threshold;
                    // a = left window under threshold
                    // b = left window over threshold
                    // c = right window under threshold
                    // d = right window over threshold
                    int a, b, c, d;
                    a = windowUpstream.getNextLowerIndex(threshold) + 1;
                    b = windowSize - a;
                    c = windowDownstream.getNextLowerIndex(threshold) + 1;
                    d = windowSize - c;

                    FishersExact fishersExact = new FishersExact(a + b + c + d);
                    double p = fishersExact.getTwoTailedP(a, b, c, d);

                    if (useML){
                        if(p <= 0.1 && d > b) {
                            Map<String, Double> info = new HashMap<>();
                            posPVal.add(new MutableTriple<>(i, p, data.getDouble(i)));
//                        info.put("pValue", p);
//                        foundPeaksNew.put(i, info);
                        }
                    } else if (p <= significanceThresh && d > b && data.getDouble(i) >= minReadNum) {
                        Map<String, Double> infos = new HashMap<>();
                        infos.put(TissFile.P_VALUE_COLUMN_NAME, p);
                        infos.put(TissFile.READ_COUNT_COLUMN_NAME, data.getDouble(i));
                        foundPeaksNew.put(i, infos);
                    }
                }
            }

            if (ref.isPlus()) {
                if (data.getDouble(i+1) <= ZERO_THRESHOLD) {
                    zeroReadsInRightWindow--;
                }
                if (data.getDouble(i+windowSize+1) <= ZERO_THRESHOLD) {
                    zeroReadsInRightWindow++;
                }
                downstreamWindowSum -= data.getDouble(i+1);
//                threshold = rightWindowSum/(double) (windowSize-1);
                downstreamWindowSum += data.getDouble(i+windowSize+1);
                threshold = downstreamWindowSum/(double) windowSize;
//                threshold = ((windowSize-1)/(double)windowSize) * threshold + (1./windowSize) * data.getDouble(i+windowSize+1);
                windowUpstream.insertSortedAndDelete(data.getDouble(i), data.getDouble(i - windowSize));
                windowDownstream.insertSortedAndDelete(data.getDouble(i + windowSize + 1), data.getDouble(i + 1));
            } else {
                if (data.getDouble(i - windowSize) <= ZERO_THRESHOLD) {
                    zeroReadsInRightWindow--;
                }
                if (data.getDouble(i) <= ZERO_THRESHOLD) {
                    zeroReadsInRightWindow++;
                }
                downstreamWindowSum -= data.getDouble(i - windowSize);
//                threshold = rightWindowSum/(double) (windowSize-1);
                downstreamWindowSum += data.getDouble(i);
//                threshold = ((windowSize-1)/(double)windowSize) * threshold + (1./windowSize) * data.getDouble(i);
                threshold = downstreamWindowSum/(double) windowSize;
                windowDownstream.insertSortedAndDelete(data.getDouble(i), data.getDouble(i - windowSize));
                windowUpstream.insertSortedAndDelete(data.getDouble(i + windowSize + 1), data.getDouble(i + 1));
            }
        }
        try {
            for (MutableTriple<Integer, Double, Double> tss : posPVal) {
                mlWriter.writeLine(ref.toPlusMinusString() + "\t" + tss.Item1 + "\t" + tss.Item2 + "\t" + tss.Item3);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private double sum(NumericArray ary, int from, int to) {
        double sum = 0;
        for (int i = from; i < to; i++) {
            sum+=ary.getDouble(i);
        }
        return sum;
    }

    private int getZeroReads(NumericArray ary, int from, int to) {
        int zeroReads = 0;
        for (int i = from; i < to; i++) {
            if (ary.getDouble(i) <= ZERO_THRESHOLD) {
                zeroReads++;
            }
        }
        return zeroReads;
    }

    @Override
    public void calculateMLResults(String prefix, boolean plot) throws IOException {
        mlWriter.close();
        Map<ReferenceSequence, List<MutableTriple<Integer, Double, Double>>> refPosPVal = new HashMap<>();
        List<PeakAndPos> mlData = new ArrayList<>();
        EI.lines(writerPath).skip(1).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            mlData.add(new PeakAndPos(Integer.parseInt(split[1]), 1.-Double.parseDouble(split[2])));
            refPosPVal.computeIfAbsent(Chromosome.obtain(split[0]), absent->new ArrayList<>()).add(new MutableTriple<>(Integer.parseInt(split[1]), Double.parseDouble(split[2]), Double.parseDouble(split[3])));
        });

        int movingAverage = (int)(mlData.size()*0.2);
        double upThresh = TiSSUtils.calculateThreshold(mlData, movingAverage);
        upThresh = 1.-upThresh;

        for (ReferenceSequence ref : refPosPVal.keySet()) {
            Map<Integer, Map<String,Double>> filtered = new HashMap<>();
            List<MutableTriple<Integer, Double, Double>> tss2val = refPosPVal.get(ref);
            if (cleanupThresh > 0) {
                tss2val = TiSSUtils.cleanUpMultiValueDataTriple(tss2val, cleanupThresh);
            }

            for (MutableTriple<Integer, Double, Double> tss : tss2val) {
                if (tss.Item2 < upThresh && tss.Item3 >= minReadNum) {
                    Map<String, Double> info = new HashMap<>();
                    info.put(TissFile.P_VALUE_COLUMN_NAME, tss.Item2);
                    info.put(TissFile.READ_COUNT_COLUMN_NAME, tss.Item3);
                    filtered.put(tss.Item1, info);
                }
            }
            res.put(ref, filtered);
        }

        if (plot) {
            plotData(prefix, prefix + "densityThresholdData.tsv", upThresh);
        }
    }

    private void plotData(String prefix, String dataFilePath, double thresh) throws IOException {
        RRunner r = new RRunner(prefix+".plotDensityThresholdData.R");
        r.set("dataFile",dataFilePath);
        r.set("pdfFile",prefix+".densityThresholds.pdf");
        r.set("customThresh", new DoubleDynamicObject(thresh));
        r.set("threshAutoParam", new DoubleDynamicObject(thresh));
        r.set("manualSelectionFile", prefix + "densityManualSelection.tsv");
        r.set("manualSelectionPdf", prefix + "densityManuelThresh.pdf");
        r.addSource(getClass().getResourceAsStream("/resources/plotPValues.R"));
        r.run(false);
    }
}
