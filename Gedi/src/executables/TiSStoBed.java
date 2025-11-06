package executables;

import gedi.iTiSS.TiSStoBedParameterSet;
import gedi.iTiSS.bedConverter.TsrToBed;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class TiSStoBed {
    public static void main(String[] args) {
        TiSStoBedParameterSet params = new TiSStoBedParameterSet();

        GediProgram tiss = GediProgram.create("TiSStoBed",
                new TsrToBed(params));

        GediProgram.run(tiss, params.paramFile, new CommandLineHandler("TiSStoBed",
                "TiSStoBed converts the final .TSRs.tsv file into bed format", args));
    }
}
