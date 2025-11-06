package gedi.iTiSS.analyzer;

import gedi.iTiSS.data.Data;
import gedi.iTiSS.modules.CountModule;
import gedi.util.io.text.LineOrientedFile;

import java.io.IOException;
import java.util.List;

public class AnalyzeCount extends AnalyzerBase {
    public AnalyzeCount(int upstreamCount, int downstreamCount, double countThreshold, List<Data> lanes) {
        for (Data lane : lanes) {
            modules.add(new CountModule(upstreamCount, downstreamCount, countThreshold, lane, "Count"));
        }
    }

    public void writeOutCountResults(LineOrientedFile outputFile) throws IOException {

        int tissCounter = this.writeOutResults(outputFile, CountModule.class);

        System.err.println("Number of TiSS with simple count model: " + tissCounter);
    }
}
