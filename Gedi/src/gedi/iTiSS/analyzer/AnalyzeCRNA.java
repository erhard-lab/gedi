package gedi.iTiSS.analyzer;

import gedi.iTiSS.data.Data;
import gedi.iTiSS.modules.InterquartileRangeModule;
import gedi.iTiSS.modules.KineticActivity;
import gedi.iTiSS.modules.TranscriptionalActivity;
import gedi.util.io.text.LineOrientedFile;

import java.io.IOException;
import java.util.List;

public class AnalyzeCRNA extends AnalyzerBase {

    public AnalyzeCRNA(int kaWindow, double kaThresh, double kaPseudoCount, double taThresh, int taWindow,
                       double taMinReadDensity, double iqrThreshold, int windowSize, double minReadDensity,
                       List<Data> iqrLanes, List<Data> taLanes, List<Data> timeSeriesLane, boolean useMM, boolean dumpMmVal) {
        for (Data lane : iqrLanes) {
            modules.add(new InterquartileRangeModule(iqrThreshold, windowSize, minReadDensity, lane, useMM, dumpMmVal, "IQR"));
        }
//        for (Data lane : taLanes) {
//            modules.add(new TranscriptionalActivity(taThresh, taWindow, taMinReadDensity, lane, "Transcriptional activity"));
//        }
//        for (Data lane : timeSeriesLane) {
//            modules.add(new KineticActivity(kaWindow, kaThresh, kaPseudoCount, lane, "Kinetic activity"));
//        }
    }

    public void writeOutAllResults(LineOrientedFile iqF, LineOrientedFile taF, LineOrientedFile kaF) throws IOException {

        int iq = this.writeOutResults(iqF, InterquartileRangeModule.class);
//        int ta = this.writeOutResults(taF, TranscriptionalActivity.class);
//        int ka = this.writeOutResults(kaF, KineticActivity.class);

        System.err.println("Number of TiSS in iq: " + iq);
//        System.err.println("Number of TiSS in ta: " + ta);
//        System.err.println("Number of TiSS in ka: " + ka);
    }

    public void writeOutAllMmData(LineOrientedFile outputFile, LineOrientedFile thresholdFile) throws IOException {
        this.writeOutMmData(outputFile, InterquartileRangeModule.class);
        this.writeOutThresholds(thresholdFile, InterquartileRangeModule.class);
    }

//    private boolean wasPeakFoundInAllLanes(Integer i, ReferenceSequence ref, List<ModuleBase> modulesToCheck) {
//        for (int j = 1; j < modulesToCheck.size(); j++) {
//            if (!modulesToCheck.get(j).getResults().get(ref).contains(i)) {
//                return false;
//            }
//        }
//        return true;
//    }
}
