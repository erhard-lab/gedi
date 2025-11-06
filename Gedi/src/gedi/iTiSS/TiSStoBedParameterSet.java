package gedi.iTiSS;

import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;

import java.io.File;

public class TiSStoBedParameterSet extends GediParameterSet {
    // General stuff
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call TiSS", false, new FileParameterType());

    public GediParameter<String> tsrFile = new GediParameter<String>(this,"tsr", "The final .TSRs.tsv file", false, new StringParameterType(), true);
    public GediParameter<Integer> minScore = new GediParameter<Integer>(this,"minScore", "min totalDelta threshold to include", false, new IntParameterType(), 0, true);

    public GediParameter<File> bedOut = new GediParameter<File>(this, "${prefix}.tsr.bed", "The converted TSRs in BED-format", false, new FileParameterType());
}
