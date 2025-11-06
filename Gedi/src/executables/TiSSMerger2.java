package executables;

import gedi.iTiSS.merger2.TissMerger2ParameterSet;
import gedi.iTiSS.merger2.TsrFileMerger;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class TiSSMerger2 {
    public static void main(String[] args) {
        TissMerger2ParameterSet params = new TissMerger2ParameterSet();

        GediProgram tiss = GediProgram.create("TiSSMerger2",
                new TsrFileMerger(params));

        GediProgram.run(tiss, params.paramFile, new CommandLineHandler("TiSSMerger2",
                "TiSSMerger2 is used to merge the TiSS found by the various TiSS modules", args));
    }
}
