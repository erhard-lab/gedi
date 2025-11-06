package gedi.iTiSS;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionStorage;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.*;

import java.io.File;

public class TiSSMergerParameterSet extends GediParameterSet {
    // General stuff
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());
    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call TiSS", false, new FileParameterType());
    // General parameters
    public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "Read data in CIT-format.", true, new StorageParameterType<AlignedReadsData>());
    public GediParameter<String> inputFiles = new GediParameter<String>(this,"in", "input files", true, new StringParameterType());
//    0:1x4x5>2  0=file at index 0, 1=file at index 1, 4=TSR-totalDelta threshold, 5=TiSS-totalDelta threshold (if omitted, only TSR-totalDelta is used), 2=output index used for next dependency
    public GediParameter<String> dependencies = new GediParameter<String>(this,"dep", "0-based indexed. number at certain position depicts the dependency of this file to the one located at the given index", false, new StringParameterType());
    public GediParameter<File> outFile = new GediParameter<File>(this, "${prefix}.tsv", "The final output of merged peaks", false, new FileParameterType());

    public GediParameter<String> pacbioIndices = new GediParameter<String>(this,"pacbio", "A string to identify pacbio files (for ex.: _xx__x -> files at index 1,2 and 5 are pacbio) underscore character (_) for skip", false, new StringParameterType(), "", true);
    public GediParameter<String> minionIndices = new GediParameter<String>(this,"minion", "A string to identify pacbio files (for ex.: _xx__x -> files at index 1,2 and 5 are minion) underscore character (_) for skip", false, new StringParameterType(), "", true);
    public GediParameter<Integer> pacbioDelta = new GediParameter<Integer>(this,"pbDelta", "Amount of positions to consider left AND right from a pacbio read starting point", false, new IntParameterType(), 5, true);
    public GediParameter<Integer> minionDelta = new GediParameter<Integer>(this,"mionDelta", "Amount of positions to consider upstream ONLY from a minion read starting point", false, new IntParameterType(), 20, true);
    public GediParameter<String> priorityDatasets = new GediParameter<String>(this,"prio", "A string to identify priority datasets (for ex.: _xx__x -> files at index 1,2 and 5 are minion) underscore character (_) for skip. Those will be used to collapse successive tiss", false, new StringParameterType(), "", true);
    public GediParameter<Integer> gapMerge = new GediParameter<Integer>(this,"gap", "The gaps between TiSS that should be merged into one", false, new IntParameterType(), 5, true);
}
