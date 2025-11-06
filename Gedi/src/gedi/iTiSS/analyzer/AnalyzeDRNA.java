package gedi.iTiSS.analyzer;

import gedi.iTiSS.data.Data;
import gedi.iTiSS.modules.XRN1module;
import gedi.util.io.text.LineOrientedFile;

import java.io.IOException;
import java.util.List;

public class AnalyzeDRNA extends AnalyzerBase {
    public AnalyzeDRNA(int windowSize, double pseudoCount, double maxWater2XRN1FoldChange, int thresholdPeakNumber,
                       double peakCallThreshold, boolean waterCheck, int peakProximity, int minReadNum, List<Data> lanes,
                       boolean useMM, boolean dumpMmVal) {
        for (Data lane : lanes) {
            modules.add(new XRN1module(windowSize, pseudoCount, maxWater2XRN1FoldChange, thresholdPeakNumber,
                    peakCallThreshold, lane, waterCheck, peakProximity, minReadNum, useMM, dumpMmVal, "XRN1"));
        }
    }

    public void writeOutXRN1Results(LineOrientedFile outputFile) throws IOException {

        int tissCounter = this.writeOutResults(outputFile, XRN1module.class);

        System.err.println("Number of TiSS in dRNA-seq: " + tissCounter);
    }

    public void writeOutAllMmData(LineOrientedFile outputFile, LineOrientedFile thresholdFile) throws IOException {
        this.writeOutMmData(outputFile, XRN1module.class);
        this.writeOutThresholds(thresholdFile, XRN1module.class);
    }
//
//    private boolean wasPeakFoundInAllLanes(Integer i, ReferenceSequence ref, List<ModuleBase> modulesToCheck) {
//        for (int j = 1; j < modulesToCheck.size(); j++) {
//            if (!modulesToCheck.get(j).getResults().get(ref).contains(i)) {
//                return false;
//            }
//        }
//        return true;
//    }
}
