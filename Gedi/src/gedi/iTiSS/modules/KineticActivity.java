package gedi.iTiSS.modules;

import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.merger2.TissFile;
import gedi.iTiSS.utils.ArrayUtils2;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.dynamic.impl.DoubleDynamicObject;
import gedi.util.dynamic.impl.IntDynamicObject;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.testing.DirichletLikelihoodRatioTest;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import gedi.util.r.RRunner;

import java.io.IOException;
import java.util.*;

public class KineticActivity extends ModuleBase {
    private int windowSize;
    private double significanceThresh;
    private double pseudoCount;
    private boolean useML;
    private int cleanupThresh;
    private boolean useUpAndDownstream;
    private String writerPath;
    private LineWriter mlWriter;
    private int minReadNum;

    public KineticActivity(int windowSize, double significanceThresh, double pseudoCount, int cleanupThresh, boolean useML, boolean useUpAndDownstream, int minReadNum, String prefix, Data lane, String name) {
        super(name, lane);
        this.windowSize = windowSize;
        this.significanceThresh = significanceThresh;
        this.pseudoCount = pseudoCount;
        this.useML = useML;
        this.cleanupThresh = cleanupThresh;
        this.useUpAndDownstream = useUpAndDownstream;
        this.minReadNum = minReadNum;
        if (useML) {
            writerPath = prefix + "kineticThresholdData.tsv";
            mlWriter = new LineOrientedFile(writerPath).write();
            try {
                mlWriter.writeLine("Ref\tPos\tValue\tReadCount");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }

    @Override
    public void findTiSS(NumericArray[] data, ReferenceSequence ref) {
        // data contains the timeseries in ascending order
        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());
        List<MutableTriple<Integer, Double, Double>> posPVal = new ArrayList<>(data[0].length()/100);

        double[] windowSum = sum(data, 0, windowSize);
        if (ref.isMinus()) {
            windowSum = sum(data, windowSize+1, windowSize*2+1);
        }
        if (useUpAndDownstream) {
            windowSum = sum(data, 0, windowSize*2+1);
            ArrayUtils2.subtract(windowSum, sum(data, windowSize, windowSize+1));
        }
        for (int i = windowSize; i < data[0].length()-(windowSize+1); i++) {
            if (useML && i%10000000 == 0) {
                System.err.print("Progress " + i + "/" + data[0].length() + ", mmData size: " + posPVal.size() + "\r");
            }
            double[] peaks = new double[windowSum.length];
            for (int j = 0; j < windowSum.length; j++) {
                peaks[j] = data[j].getDouble(i);
            }
            double[] windowMeans = new double[windowSum.length];
            for (int j = 0; j < windowSum.length; j++) {
                windowMeans[j] = windowSum[j]/(useUpAndDownstream ? windowSize*2 : windowSize);
            }
            double p = DirichletLikelihoodRatioTest.testMultinomials(pseudoCount, windowMeans, peaks);
            if (useML) {
                if (p < 0.1) {
//                Map<String, Double> info = new HashMap<>();
//                info.put("pValue", p);
//                foundPeaksNew.put(i, info);
                    posPVal.add(new MutableTriple<>(i, p, ArrayUtils.max(peaks)));
                }
            } else if (p <= significanceThresh && ArrayUtils.max(peaks) >= minReadNum) {
                Map<String,Double> infos = new HashMap<>();
                infos.put(TissFile.P_VALUE_COLUMN_NAME, p);
                infos.put(TissFile.READ_COUNT_COLUMN_NAME, ArrayUtils.max(peaks));
                foundPeaksNew.put(i, infos);
            }

            for (int j = 0; j < windowSum.length; j++) {
                double subtract = ref.isPlus() || useUpAndDownstream ? data[j].getDouble(i-windowSize) : data[j].getDouble(i+1);
                double add = ref.isPlus() && !useUpAndDownstream ? data[j].getDouble(i) : data[j].getDouble(i+windowSize+1);
                windowSum[j] -= subtract;
                windowSum[j] += add;
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
        double sum = 0.;
        for (int i = from; i < to; i++) {
            sum += ary.getDouble(i);
        }
        return sum;
    }

    private double[] sum(NumericArray[] ary, int from, int to) {
        double[] sum = new double[ary.length];
        for (int i = 0; i < sum.length; i++) {
            sum[i] = sum(ary[i], from, to);
        }
        return sum;
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
            plotData(prefix, prefix + "kineticThresholdData.tsv", upThresh);
        }
    }

    private void plotData(String prefix, String dataFilePath, double thresh) throws IOException {
        RRunner r = new RRunner(prefix+".plotKineticThresholdData.R");
        r.set("dataFile",dataFilePath);
        r.set("pdfFile",prefix+".kineticThresholds.pdf");
        r.set("customThresh", new DoubleDynamicObject(thresh));
        r.set("threshAutoParam", new DoubleDynamicObject(thresh));
        r.set("manualSelectionFile", prefix + "kineticManualSelection.tsv");
        r.set("manualSelectionPdf", prefix + "kineticManuelThresh.pdf");
        r.addSource(getClass().getResourceAsStream("/resources/plotPValues.R"));
        r.run(false);
    }
}
