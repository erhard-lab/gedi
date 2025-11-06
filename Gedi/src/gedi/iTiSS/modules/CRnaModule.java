package gedi.iTiSS.modules;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.merger2.TissFile;
import gedi.iTiSS.merger2.TsrFile;
import gedi.iTiSS.merger2.TsrFileMerger;
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

import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.IOException;
import java.util.*;

public class CRnaModule extends ModuleBase{

    private double zScoreThresh;
    private int windowSize;
    private boolean useMM;
    private double pseudoCount;
    private int cleanupThresh;
    private int minReadNum;

    public CRnaModule(double zScoreThresh, double pseudoCount, int windowSize, int cleanupThresh, int minReadNum, Data lane,
                      boolean useMM, String name) {
        super(name, lane);
        this.zScoreThresh = zScoreThresh;
        this.windowSize = windowSize;
        this.pseudoCount = pseudoCount;
        this.useMM = useMM;
        this.cleanupThresh = cleanupThresh;
        this.minReadNum = minReadNum;
    }

    @Override
    public void findTiSS(NumericArray[] inData, ReferenceSequence ref) {
        // Should only be one lane, the total lane
        NumericArray data = inData[0];

        List<PeakAndPos> mmDataL = useMM ? new ArrayList<>() : null;
        List<PeakAndPos> mmDataR = useMM ? new ArrayList<>() : null;

        double[] upstream = data.toDoubleArray(0, windowSize);
        double[] downstream = data.toDoubleArray(windowSize+1, windowSize*2+1);
        for (int i = 0; i < windowSize; i++) {
            upstream[i] = upstream[i] + pseudoCount;
            downstream[i] = downstream[i] + pseudoCount;
        }

        StaticSizeSortedArrayList<Double> windowUpstream = new StaticSizeSortedArrayList<>(EI.wrap(upstream).toList(), Double::compareTo);
        StaticSizeSortedArrayList<Double> windowDownstream = new StaticSizeSortedArrayList<>(EI.wrap(downstream).toList(), Double::compareTo);

        if (ref.isMinus()) {
            StaticSizeSortedArrayList<Double> windowTmp = windowUpstream;
            windowUpstream = windowDownstream;
            windowDownstream = windowTmp;
        }

        Map<Integer, Map<String, Double>> foundPeaksNew = this.res.computeIfAbsent(ref, k -> new HashMap<>());

        double sumUpstream = sum(windowUpstream);
        double sumDownstream = sum(windowDownstream);
        double meanUpstream = sumUpstream/(windowUpstream.getSize()-1);
        double meanDownstream = sumDownstream/(windowDownstream.getSize()-1);
        double sdUpsteam = sampleSd(windowUpstream,sumUpstream/windowUpstream.getSize());
        double sdDownsteam = sampleSd(windowDownstream, sumDownstream/windowDownstream.getSize());

        for (int i = windowSize; i < data.length() - (windowSize+1); i++) {
            if (mmDataL != null && i%10000000 == 0) {
                //System.err.print("Progress " + i + "/" + data.length() + ", mmDataL size: " + mmDataL.size() + ", mmDataR size: " + mmDataR.size() + ", upstream: " + sumUpstream + ", downstream: " + sumDownstream + "\r");
                System.err.print(String.format("Progress %d / %d, mmDataL size: %d, mmDataR size: %d, upstream: %.2f, downstream: %.2f \r", i, data.length(), mmDataL.size(), mmDataR.size(), sumUpstream, sumDownstream));
            }
            double valueOfInterest = data.getDouble(i) + pseudoCount;
            if (sdUpsteam <= pseudoCount) {
                sdUpsteam = pseudoCount;
            }
            if (sdDownsteam <= pseudoCount) {
                sdDownsteam = pseudoCount;
            }

            double upstreamZ = (valueOfInterest - meanUpstream) / sdUpsteam;
            double downstreamZ = (valueOfInterest - meanDownstream) / sdDownsteam;
            if (useMM && mmDataL != null) {
                if (upstreamZ > 2 && downstreamZ > 2) {
                    mmDataL.add(new PeakAndPos(i, upstreamZ));
                    mmDataR.add(new PeakAndPos(i, downstreamZ));
                }
            } else {
                if (upstreamZ > zScoreThresh && downstreamZ > zScoreThresh && data.getDouble(i) >= minReadNum) {
                    Map<String, Double> infos = new HashMap<>();
                    infos.put("zScore-US", upstreamZ);
                    infos.put("zScore-DS", downstreamZ);
                    infos.put(TissFile.Z_SCORE_COLUMN_NAME, upstreamZ);
                    infos.put(TissFile.READ_COUNT_COLUMN_NAME, data.getDouble(i));
                    foundPeaksNew.put(i, infos);
                }
            }

            boolean indelSuccess;
            if (ref.isPlus()) {
                indelSuccess = windowUpstream.insertSortedAndDelete(data.getDouble(i) + pseudoCount, data.getDouble(i - windowSize) + pseudoCount);
                indelSuccess &= windowDownstream.insertSortedAndDelete(data.getDouble(i + windowSize + 1) + pseudoCount, data.getDouble(i + 1) + pseudoCount);
                sumUpstream = (sumUpstream - (data.getDouble(i-windowSize)+pseudoCount)) + data.getDouble(i)+pseudoCount;
                sumDownstream = (sumDownstream - (data.getDouble(i+1)+pseudoCount)) + data.getDouble(i+windowSize+1)+pseudoCount;
            } else {
                indelSuccess = windowDownstream.insertSortedAndDelete(data.getDouble(i) + pseudoCount, data.getDouble(i - windowSize) + pseudoCount);
                indelSuccess &= windowUpstream.insertSortedAndDelete(data.getDouble(i + windowSize + 1) + pseudoCount, data.getDouble(i + 1) + pseudoCount);
                sumDownstream = (sumUpstream - (data.getDouble(i-windowSize)+pseudoCount)) + data.getDouble(i)+pseudoCount;
                sumUpstream = (sumDownstream - (data.getDouble(i+1)+pseudoCount)) + data.getDouble(i+windowSize+1)+pseudoCount;
            }

            if (!indelSuccess) {
                System.err.println("Error while trying to delete value. Index: " + i);
            }

            meanUpstream = sumUpstream/(windowUpstream.getSize()-1);
            meanDownstream = sumDownstream/(windowDownstream.getSize()-1);
            sdUpsteam = sampleSd(windowUpstream,sumUpstream/windowUpstream.getSize());
            sdDownsteam = sampleSd(windowDownstream, sumDownstream/windowDownstream.getSize());
        }

        if (mmDataL != null) {
            for (int i = 0; i < mmDataL.size(); i++) {
                Map<String, Double> info = new HashMap<>();
                info.put("zScore-US", mmDataL.get(i).getValue());
                info.put("zScore-DS", mmDataR.get(i).getValue());
                info.put(TissFile.Z_SCORE_COLUMN_NAME, mmDataL.get(i).getValue());
                info.put(TissFile.READ_COUNT_COLUMN_NAME, data.getDouble(mmDataL.get(i).getPos()));
                foundPeaksNew.put(mmDataL.get(i).getPos(), info);
            }
        }
    }

    private double sum(StaticSizeSortedArrayList<Double> lst) {
        double sum = 0;
        for (int i = 0; i < lst.getSize(); i++) {
            sum += lst.getValueAtIndex(i);
        }
        return sum;
    }

    private double sampleSd(StaticSizeSortedArrayList<Double> lst, double mean) {
        double sum = 0;
        for (int i = 0; i < lst.getSize(); i++) {
            sum += Math.pow(lst.getValueAtIndex(i) - mean, 2);
        }
        return Math.sqrt(sum/(lst.getSize()-1));
    }

    @Override
    public void calculateMLResults(String prefix, boolean plot) throws IOException {
        LineWriter writerup = new LineOrientedFile(prefix + "densePeakThresholdData.tsv").write();
        writerup.writeLine("Ref\tPos\tValue1\tValue2\t" + TissFile.Z_SCORE_COLUMN_NAME + "\t" + TissFile.READ_COUNT_COLUMN_NAME);

        List<PeakAndPos> mlDataUp = new ArrayList<>();
        List<PeakAndPos> mlDataDown = new ArrayList<>();
        for (ReferenceSequence ref : res.keySet()) {
            Map<Integer, Map<String, Double>> tss = res.get(ref);
            for (int tssPos : tss.keySet()) {
                Map<String, Double> info = tss.get(tssPos);
                mlDataUp.add(new PeakAndPos(tssPos, info.get("zScore-US")));
                mlDataDown.add(new PeakAndPos(tssPos, info.get("zScore-DS")));
                writerup.writeLine(ref.toPlusMinusString() + "\t" + tssPos + "\t" + info.get("zScore-US") + "\t" + info.get("zScore-DS") + "\t" + info.get(TissFile.Z_SCORE_COLUMN_NAME) + "\t" + info.get(TissFile.READ_COUNT_COLUMN_NAME));
            }
        }
        writerup.close();
        int movingAverageUp = (int)(mlDataUp.size()*0.2);
        int movingAverageDown = (int)(mlDataUp.size()*0.2);
        double upThresh = TiSSUtils.calculateThreshold(mlDataUp, movingAverageUp);
        double downThresh = TiSSUtils.calculateThreshold(mlDataDown, movingAverageDown);

        for (ReferenceSequence ref : res.keySet()) {
            Map<Integer, Map<String,Double>> filtered = new HashMap<>();
            Map<Integer, Map<String, Double>> tss2val = res.get(ref);
            if (cleanupThresh > 0) {
                List<MutableTriple<Integer, Double, Double>> upstream = EI.wrap(res.get(ref).keySet()).map(t -> new MutableTriple<>(t, res.get(ref).get(t).get("zScore-US"), res.get(ref).get(t).get(TissFile.READ_COUNT_COLUMN_NAME))).list();
                List<MutableTriple<Integer, Double, Double>> downstream = EI.wrap(res.get(ref).keySet()).map(t -> new MutableTriple<>(t, res.get(ref).get(t).get("zScore-DS"), res.get(ref).get(t).get(TissFile.READ_COUNT_COLUMN_NAME))).list();
                upstream = TiSSUtils.cleanUpMultiValueDataTriple(upstream, cleanupThresh);
                downstream = TiSSUtils.cleanUpMultiValueDataTriple(downstream, cleanupThresh);
                Map<Integer,Map<String,Double>> tss2valFiltered = new HashMap<>();
                EI.wrap(upstream).forEachRemaining(u -> {
                    Map<String, Double> map = tss2valFiltered.computeIfAbsent(u.Item1, absent -> new HashMap<>());
                    map.put("zScore-US", u.Item2);
                    map.put(TissFile.Z_SCORE_COLUMN_NAME, u.Item2);
                    map.put(TissFile.READ_COUNT_COLUMN_NAME, u.Item3);
                });
                EI.wrap(downstream).forEachRemaining(u -> {
                    Map<String, Double> map = tss2valFiltered.computeIfAbsent(u.Item1, absent -> new HashMap<>());
                    map.put("zScore-DS", u.Item2);
                });
                tss2val = tss2valFiltered;
            }
            for (int tss : tss2val.keySet()) {
                if (tss2val.get(tss).containsKey("zScore-US") && tss2val.get(tss).containsKey("zScore-DS") && tss2val.get(tss).get("zScore-US") > upThresh && tss2val.get(tss).get("zScore-DS") > downThresh && tss2val.get(tss).get(TissFile.READ_COUNT_COLUMN_NAME) >= minReadNum) {
                    filtered.put(tss, tss2val.get(tss));
                }
            }
            res.put(ref, filtered);
        }

        if (plot) {
            plotData(prefix, prefix + "densePeakThresholdData.tsv", upThresh, downThresh);
        }
    }

    private void plotData(String prefix, String dataFilePath, double threshUp, double threshDown) throws IOException {
        RRunner r = new RRunner(prefix+".plotDensePeakThresholdData.R");
        r.set("dataFile",dataFilePath);
        r.set("pdfFile",prefix+".densePeakThresholds.pdf");
        r.set("customThreshUp", new DoubleDynamicObject(threshUp));
        r.set("customThreshDown", new DoubleDynamicObject(threshDown));
        r.set("upThreshAutoParam", new DoubleDynamicObject(threshUp));
        r.set("downThreshAutoParam", new DoubleDynamicObject(threshDown));
        r.set("manualSelectionFile", prefix + "densePeakManualSelection.tsv");
        r.set("manualSelectionPdf", prefix + "densePeakManuelThresh.pdf");
        r.addSource(getClass().getResourceAsStream("/resources/plotDenseThresh.R"));
        r.run(false);
    }
}
