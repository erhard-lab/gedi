package executables;

import gedi.iTiSS.TiSSMergerParameterSet;
import gedi.iTiSS.merger.MergeTiSS;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class TiSSMerger {
    public static void main(String[] args) {
        TiSSMergerParameterSet params = new TiSSMergerParameterSet();

        GediProgram tiss = GediProgram.create("TiSSMerger",
                new MergeTiSS(params));

        GediProgram.run(tiss, params.paramFile, new CommandLineHandler("TiSSMerger",
                "TiSSMerger is used to merge the TiSS found by the various TiSS modules", args));
    }
}
