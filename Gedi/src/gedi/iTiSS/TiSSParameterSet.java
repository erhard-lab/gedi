package gedi.iTiSS;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.iTiSS.utils.AnalysisModuleType;
import gedi.iTiSS.utils.AnalyzeType;
import gedi.iTiSS.utils.ReadType;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.*;

import java.io.File;

public class TiSSParameterSet extends GediParameterSet {
    // General stuff
    public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
    public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call TiSS", false, new FileParameterType());
    public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());

    // General parameters
    public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "Read data in CIT-format.", true, new StorageParameterType<AlignedReadsData>());
    public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"genomic", "The indexed GEDI genome.", true, new GenomicParameterType());
    public GediParameter<AnalysisModuleType> analyzeModuleType = new GediParameter<AnalysisModuleType>(this, "modType", "The type of analyzis [DENSITY, KINETIC, DENSE_PEAK, SPARSE_PEAK]", true, new EnumParameterType<>(AnalysisModuleType.class));
    public GediParameter<ReadType> readType = new GediParameter<ReadType>(this, "readType", "How to count reads, three-prime not yet supported [FIVE_PRIME, DENSITY, THREE_PRIME]", false, new EnumParameterType<>(ReadType.class), ReadType.FIVE_PRIME, true);
    public GediParameter<Integer> wSize = new GediParameter<Integer>(this,"windowSize", "The size of the moving window", false, new IntParameterType(), 100, true);
    public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Which strandness.", false, new EnumParameterType<>(Strandness.class), Strandness.Sense, true);
    public GediParameter<Double> pseudoCount = new GediParameter<Double>(this,"pseudo", "A pseudo count be added to each position", false, new DoubleParameterType(), 1., true);
    public GediParameter<Integer> cleanupThresh = new GediParameter<Integer>(this,"cleanupThresh", "Threshold at which multi-occurrences of a value will be filtered out", false, new IntParameterType(), 100, true);
    public GediParameter<String> testChromosomes = new GediParameter<String>(this,"testChr", "The chromosomes to use (for testing purposes, individual chromosomes separated by comma, i.e. 1+,1-,...)", false, new StringParameterType(), true);
    public GediParameter<Integer> minReadNum = new GediParameter<Integer>(this,"minReadNum", "Minimum amount of reads to call a TiSS", false, new IntParameterType(), 0, true);

    // GenomicCreate
    public GediParameter<String> chromSizes = new GediParameter<String>(this,"chromSizes", "tsv file containing the sizes of each chromosome", true, new StringParameterType(), true);
    public GediParameter<String> bamFiles = new GediParameter<String>(this,"bams", "Whitespace separated list of bam-files", true, new StringParameterType(), true);

    // dense_peak parameters
    public GediParameter<String> replicates = new GediParameter<String>(this,"rep", "A string to identify samples to combine (for ex.: XX_X -> combines read counts from sample 0, 1 and 3 with 2 being ignored) underscore character (_) for skip", false, new StringParameterType(), true);
    public GediParameter<String> timecourses = new GediParameter<String>(this,"timepoints", "A string to identify conditions (for ex.: 12_3 -> 3 conditions at index 0, 1 and 3, with 2 being ignored) underscore character (_) for skip", false, new StringParameterType(), "", true);
    public GediParameter<Double> iqr = new GediParameter<Double>(this,"zscore", "z-totalDelta threshold to call a TiSS (not used if -autoparam is set)", false, new DoubleParameterType(), 5., true);
    public GediParameter<Double> pValThresh = new GediParameter<Double>(this,"pVal", "p-Value threshold to call a TiSS (not used if -autoparam is set)", false, new DoubleParameterType(), 0.0001, true);
    public GediParameter<Double> minReadDens = new GediParameter<Double>(this,"minReadDens", "The minimum read density to look for", false, new DoubleParameterType(), 0., true);
    public GediParameter<File> outIQR = new GediParameter<File>(this, "${prefix}.densePeak.tsv", "The final output of called peaks for the dense-peak analysis", false, new FileParameterType());
    public GediParameter<File> outTA = new GediParameter<File>(this, "${prefix}.density.tsv", "The final output of called TiSS for density-changes analysis", false, new FileParameterType());
    public GediParameter<File> outKA = new GediParameter<File>(this, "${prefix}.kinetic.tsv", "The final output of called TiSS for kinetic analysis", false, new FileParameterType());

    // Kinetic
    public GediParameter<Boolean> useUpAndDownstream = new GediParameter<Boolean>(this, "upAndDown", "[Kinetic] wrap window around TiSS instead of only upstream", false, new BooleanParameterType(), false, true);

    // parameter estimation
    public GediParameter<Boolean> useAutoparam = new GediParameter<Boolean>(this, "autoparam", "automatically set thresholds based on the data", false, new BooleanParameterType(), false, true);
    public GediParameter<Boolean> plotMM = new GediParameter<Boolean>(this, "plotParams", "plot thresholds", false, new BooleanParameterType(), false, true);

    // sparse_peak parameters
    public GediParameter<Double> peakFCThreshold = new GediParameter<Double>(this,"peakFC", "The fold-change threshold for sparse peak algorithm (not used if -autoparam is set)", false, new DoubleParameterType(), 4., true);
    public GediParameter<File> outXRN1 = new GediParameter<File>(this, "${prefix}.sparsePeak.tsv", "The final output of called TiSS of the sparse_peak analysis", false, new FileParameterType());

    // equal transcription parameters
    public GediParameter<File> outEquTrans = new GediParameter<File>(this, "${prefix}.equalTranscription.tsv", "The final output of called TiSS of the equal_transcription analysis", false, new FileParameterType());

    // currently unused
    public GediParameter<Integer> peakProximityNumber = new GediParameter<Integer>(this,"proxPeak", "The window size for peaks to be seen in close proximity to the current peak", false, new IntParameterType(), 3, true);
    public GediParameter<Integer> peakCompareNumber = new GediParameter<Integer>(this,"compPeak", "The x'th largest peak in the window the largest one is compared to", false, new IntParameterType(), 4, true);
    public GediParameter<String> waterXrn1Pos = new GediParameter<String>(this,"wxPos", "A string to identifying water and xrn1 samples. w for water x for xrn1. underscore character (_) for skip", false, new StringParameterType(), true);
    public GediParameter<String> waterXrn1Pairs = new GediParameter<String>(this,"wxPair", "A string to identifying water and xrn1 pairs. Identical characters form a pair (abab means pos 0/2 and 1/3 are a pair). underscore character (_) for skip. Needs to be in check with wxPos", false, new StringParameterType(), true);
    public GediParameter<Boolean> dumpAllValues = new GediParameter<Boolean>(this, "dumpMM", "dump all values used for machine learning the thresholds", false, new BooleanParameterType(), false, true);
    public GediParameter<AnalyzeType> analyzeType = new GediParameter<AnalyzeType>(this, "type", "The type of analyzis [CRNA, DRNA]", false, new EnumParameterType<>(AnalyzeType.class), true);
    public GediParameter<Double> water2xrn1FC = new GediParameter<Double>(this,"waterFC", "A fold change that XRN1 needs to exceed water on the same position", false, new DoubleParameterType(), 2., true);

    public GediParameter<Integer> upstream = new GediParameter<Integer>(this,"upstream", "The size of the upstream window", false, new IntParameterType(), 5, true);
    public GediParameter<Integer> downstream = new GediParameter<Integer>(this,"downstream", "The size of the downstream window", false, new IntParameterType(), 5, true);
    public GediParameter<Double> countThreshold = new GediParameter<Double>(this,"countThresh", "The threshold of read-start sites inside the upstream+pos+downstream window", false, new DoubleParameterType(), 5., true);
    public GediParameter<File> outCount = new GediParameter<File>(this, "${prefix}.count.tsv", "The final output of called peaks for Count", false, new FileParameterType());
    public GediParameter<File> outTsr = new GediParameter<File>(this, "${prefix}.tsr.tsv", "The final output of called peaks for TSR", false, new FileParameterType());
}
