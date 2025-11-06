package executables;

import gedi.iTiSS.TiSSController;
import gedi.iTiSS.TiSSParameterSet;
import gedi.iTiSS.utils.nonGediCompatible.ConvertBam;
import gedi.iTiSS.utils.nonGediCompatible.GenomicCreate;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class TiSS {
    public static void main(String[] args) {
        TiSSParameterSet params = new TiSSParameterSet();

        GediProgram tiss = GediProgram.create("TiSS",
                new GenomicCreate(params),
                new ConvertBam(params),
                new TiSSController(params));

        GediProgram.run(tiss, params.paramFile, new CommandLineHandler("TiSS",
                "TiSS is used for TiSS-detection in TiSS-profiling datasets", args));
    }
}
