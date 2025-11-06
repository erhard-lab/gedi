package gedi.iTiSS.analyzer;

import gedi.iTiSS.data.Data;
import gedi.iTiSS.modules.TsrModule;
import gedi.iTiSS.modules.XRN1module;
import gedi.util.io.text.LineOrientedFile;

import java.io.IOException;
import java.util.List;

public class AnalyzePROseq extends AnalyzerBase {
    public AnalyzePROseq(int tsrSize, int minReads, int bufferUpstream, int bufferDownstream,
                         int windowSize, double pseudoCount, int thresholdPeakNumber, double peakCallThreshold,
                         int peakProximity, int minReadNum, List<Data> lanes, boolean useMM, boolean dumpMmVals) {
        for (Data data : lanes) {
            modules.add(new TsrModule(tsrSize, minReads, bufferUpstream, bufferDownstream, data, "TSR-module"));
            modules.add(new XRN1module(windowSize, pseudoCount, Double.NaN, thresholdPeakNumber, peakCallThreshold, data, false, peakProximity, minReadNum, useMM, dumpMmVals, "XRN1-module"));
        }
    }

    public void writeOutPROseqResults(LineOrientedFile tsrOut, LineOrientedFile xrnOut) throws IOException {
        int tissCount1 = this.writeOutResults(tsrOut, TsrModule.class);
        int tissCount2 = this.writeOutResults(xrnOut, XRN1module.class);

        System.err.println("Number of TiSS with TSR-model: " + tissCount1);
        System.err.println("Number of TiSS with XRN1-model: " + tissCount2);
    }
}
