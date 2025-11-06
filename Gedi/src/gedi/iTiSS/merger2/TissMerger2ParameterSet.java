package gedi.iTiSS.merger2;

import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;

import java.io.File;

public class TissMerger2ParameterSet extends GediParameterSet {
    // General stuff
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call TiSS", false, new FileParameterType());
    public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());

    public GediParameter<String> inTissFiles = new GediParameter<String>(this,"inTiss", "input files from iTiSS", true, new StringParameterType(), null, true);
    public GediParameter<String> inTsrFiles = new GediParameter<String>(this,"inTsr", "input files from TissMerger2", true, new StringParameterType(), null, true);
    public GediParameter<Integer> extension = new GediParameter<Integer>(this,"ext", "Amount of positions to consider left AND right", false, new IntParameterType(), 5, true);
    public GediParameter<Integer> gap = new GediParameter<Integer>(this,"gap", "When to consider two TiSS of a single file to be combined into a TSR", false, new IntParameterType(), 5, true);
    public GediParameter<Boolean> keepFileIds = new GediParameter<>(this, "keepFileIds", "Keep file-associations inside TSR files", false, new BooleanParameterType());

    public GediParameter<Integer> minScore = new GediParameter<Integer>(this,"minScore", "Minimum number of datasets to support a TSR. Removed otherwise", false, new IntParameterType(), 1, true);

    public GediParameter<File> outFile = new GediParameter<File>(this, "${prefix}.tsr", "The final output of merged peaks", false, new FileParameterType());
}
